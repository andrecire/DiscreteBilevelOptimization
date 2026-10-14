# ---------------------------------------------
# Create DD file collection
# ---------------------------------------------

import argparse
import os


def parse_aux_file(filepath):
    aux_info = {
        "num_lower_vars": None,
        "num_lower_constraints": None,
        "lower_level_variables": [],
        "lower_level_constraints": [],
        "instance_name": None,
        "mps_file": None,
        "lp_file": None
    }

    def matches(tag, variants):
        return tag.upper() in [v.upper() for v in variants]

    found_vars_section = False
    found_constraints_section = False

    with open(filepath, 'r') as file:
        lines = [line.strip() for line in file if line.strip() and not line.strip().startswith("#")]

    i = 0
    while i < len(lines):
        line = lines[i]

        if matches(line, ["@NUMVARS"]):
            i += 1
            aux_info["num_lower_vars"] = int(lines[i])
        elif matches(line, ["@NUMCONSTR", "@NUMCONSTRS"]):
            i += 1
            aux_info["num_lower_constraints"] = int(lines[i])
        elif matches(line, ["@NAME"]):
            i += 1
            aux_info["instance_name"] = lines[i]
        elif matches(line, ["@MPS"]):
            i += 1
            aux_info["mps_file"] = lines[i]
        elif matches(line, ["@LP"]):
            i += 1
            aux_info["lp_file"] = lines[i]
        elif matches(line, ["@VARSBEGIN"]):
            found_vars_section = True
            i += 1
            while i < len(lines) and not matches(lines[i], ["@VARSEND"]):
                parts = lines[i].split()
                if len(parts) != 2:
                    raise ValueError(f"Invalid variable line: '{lines[i]}'")
                var_name, obj_coeff = parts
                aux_info["lower_level_variables"].append((var_name, float(obj_coeff)))
                i += 1
        elif matches(line, ["@CONSTRBEGIN", "@CONSTRSBEGIN"]):
            found_constraints_section = True
            i += 1
            while i < len(lines) and not matches(lines[i], ["@CONSTREND", "@CONSTRSEND"]):
                aux_info["lower_level_constraints"].append(lines[i])
                i += 1
        i += 1

    # Post-processing validation
    if not found_vars_section:
        raise ValueError("Missing @VARSBEGIN ... @VARSEND section in auxiliary file.")
    if not found_constraints_section:
        raise ValueError("Missing @CONSTRSBEGIN ... @CONSTRSEND section in auxiliary file.")
    if aux_info["mps_file"] is None:
        raise ValueError("Missing @MPS definition in auxiliary file.")

    return aux_info


def infer_instance_name(mps_path):
    filename = os.path.basename(mps_path)
    if filename.endswith(".mps.gz"):
        return filename[:-7]
    if filename.endswith(".mps"):
        return filename[:-4]
    return os.path.splitext(filename)[0]


def normalize_mps_name(path):
    return infer_instance_name(os.path.basename(path))


def write_dd_file(aux_info, mps_path, output_dir):
    import gurobipy as gp

    model = gp.read(mps_path)
    model.setParam("OutputFlag", 0)

    # Check if maximization problem
    if model.ModelSense == gp.GRB.MAXIMIZE:
        raise ValueError("This script only supports minimization problems.")

    lower_vars = {v for v, _ in aux_info["lower_level_variables"]}
    lower_cons = set(aux_info["lower_level_constraints"])

    follower_vars = [v for v in model.getVars() if v.VarName in lower_vars]
    leader_vars = [v for v in model.getVars() if v.VarName not in lower_vars]

    follower_cons = [c for c in model.getConstrs() if c.ConstrName in lower_cons]
    leader_cons = [c for c in model.getConstrs() if c.ConstrName not in lower_cons]

    def build_normalized_matrix(constraints, variables):
        """
        Constructs the constraint coefficient matrix with all constraints normalized to <=.
        Equality constraints are duplicated as two inequalities.
        """
        matrix = []
        for c in constraints:
            row = [model.getCoeff(c, v) for v in variables]
            if c.Sense == gp.GRB.GREATER_EQUAL:
                row = [-coeff if coeff != 0 else 0 for coeff in row]
            elif c.Sense == gp.GRB.EQUAL:
                matrix.append(row[:])  # A * x <= b
                matrix.append([-coeff if coeff != 0 else 0 for coeff in row])  # -A * x <= -b
                continue
            matrix.append(row)
        return matrix

    def build_normalized_rhs(constraints):
        """
        Build a normalized RHS vector corresponding to the <= constraints.

        '=' constraints are split into two entries (b and -b),
        and '>=' constraints are negated.

        Returns:
            rhs: list of RHS values corresponding to normalized constraints
        """
        rhs = []

        for c in constraints:
            if c.Sense == gp.GRB.LESS_EQUAL:
                rhs.append(c.RHS)

            elif c.Sense == gp.GRB.GREATER_EQUAL:
                rhs.append(-c.RHS)

            elif c.Sense == gp.GRB.EQUAL:
                rhs.append(c.RHS)
                rhs.append(-c.RHS)

            else:
                raise ValueError(f"Unsupported constraint sense: {c.Sense}")

        return rhs

    Aleader = build_normalized_matrix(leader_cons, leader_vars)
    Bleader = build_normalized_matrix(leader_cons, follower_vars)
    Afol    = build_normalized_matrix(follower_cons, leader_vars)
    Bfol    = build_normalized_matrix(follower_cons, follower_vars)

    # b      = [c.RHS for c in leader_cons]
    # bfol   = [c.RHS for c in follower_cons]
    b = build_normalized_rhs(leader_cons)
    bfol = build_normalized_rhs(follower_cons)

    # Objective from model (leader problem)
    obj = model.getObjective()
    var_to_coeff = {obj.getVar(i).VarName: obj.getCoeff(i) for i in range(obj.size())}
    c1leader = [var_to_coeff.get(v.VarName, 0) for v in leader_vars]
    c2leader = [var_to_coeff.get(v.VarName, 0) for v in follower_vars]

    # Follower objective from AUX file
    cfol_map = dict(aux_info["lower_level_variables"])
    cfol = [cfol_map.get(v.VarName, 0) for v in follower_vars]

    nVarL = len(leader_vars)
    nVarF = len(follower_vars)
    nConsL = len(Aleader)
    nConsF = len(Bfol)

    # Function to obtain variable types as strings
    def var_type_str(v):
        if v.VType == gp.GRB.BINARY:
            return "bin"
        elif v.VType in [gp.GRB.INTEGER, gp.GRB.SEMIINT]:
            return "int"
        else:
            return "cont"
        
    # Function to obtain lower and upper bounds as strings
    def var_bounds_str(v):
        if v.LB > -gp.GRB.INFINITY and v.UB < gp.GRB.INFINITY:
            return f"{v.LB} {v.UB}"
        elif v.LB > -gp.GRB.INFINITY:
            return f"{v.LB} inf"
        elif v.UB < gp.GRB.INFINITY:
            return f"-inf {v.UB}"
        else:
            return "-inf inf"
        
    dd_path = os.path.join(output_dir, f"{infer_instance_name(mps_path)}.dd")

    os.makedirs(os.path.dirname(dd_path), exist_ok=True)
    with open(dd_path, "w") as f:
        f.write(f"{nVarL} {nVarF} {nConsL} {nConsF}\n")
        for row in Aleader: f.write(" ".join(map(str, row)) + "\n")
        for row in Bleader: f.write(" ".join(map(str, row)) + "\n")
        for row in Afol:    f.write(" ".join(map(str, row)) + "\n")
        for row in Bfol:    f.write(" ".join(map(str, row)) + "\n")
        f.write(" ".join(map(str, c1leader)) + "\n")
        f.write(" ".join(map(str, c2leader)) + "\n")
        f.write(" ".join(map(str, cfol)) + "\n")
        f.write(" ".join(map(str, b)) + "\n")
        f.write(" ".join(map(str, bfol)) + "\n")

        # Create variable bounds
        bounds_leader_str = [var_bounds_str(v) for v in leader_vars]
        bounds_follower_str = [var_bounds_str(v) for v in follower_vars]
        f.write(" ".join(bounds_leader_str) + "\n")
        f.write(" ".join(bounds_follower_str) + "\n")

        # Create types 
        types_leader_str = [var_type_str(v) for v in leader_vars]
        types_follower_str = [var_type_str(v) for v in follower_vars]

        f.write(" ".join(types_leader_str) + "\n")
        f.write(" ".join(types_follower_str) + "\n")

        # Extract original Gurobi variable names
        leader_varnames = [v.VarName for v in leader_vars]
        follower_varnames = [v.VarName for v in follower_vars]

        # Write to file
        f.write(" ".join(leader_varnames) + "\n")
        f.write(" ".join(follower_varnames) + "\n")

    print(f"Written instance to: {dd_path}")

    # # Print LP structure for verification
    # print("\n=== LEADER PROBLEM ===")
    # print(f"Objective: min c1leader^T x + c2leader^T y")
    # print("c1leader:", c1leader)
    # print("c2leader:", c2leader)
    # print("Aleader (A x):")
    # for row in Aleader:
    #     print("  ", row)
    # print("Bleader (B y):")
    # for row in Bleader:
    #     print("  ", row)
    # print("RHS (b):", b)

    # print("\n=== FOLLOWER PROBLEM ===")
    # print(f"Objective: min cfol^T y")
    # print("cfol:", cfol)
    # print("Afol (A x):")
    # for row in Afol:
    #     print("  ", row)
    # print("Bfol (B y):")
    # for row in Bfol:
    #     print("  ", row)
    # print("RHS (bfol):", bfol)

    # print("Leader variable bounds:", " ".join(bounds_leader_str))
    # print("Follower variable bounds:", " ".join(bounds_follower_str))

    # print("Leader variable types:", " ".join(types_leader_str))
    # print("Follower variable types:", " ".join(types_follower_str))

def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert one bilevel instance from .mps/.mps.gz + .aux into .dd format."
    )
    parser.add_argument(
        "--mps_file",
        required=True,
        help="Path to the input .mps or .mps.gz file.",
    )
    parser.add_argument(
        "--aux_file",
        required=True,
        help="Path to the input .aux file.",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory where the output .dd file should be created.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    aux_info = parse_aux_file(args.aux_file)

    aux_mps_name = normalize_mps_name(aux_info["mps_file"])
    provided_mps_name = normalize_mps_name(args.mps_file)
    if aux_mps_name != provided_mps_name:
        raise ValueError(
            f"AUX file expects MPS '{aux_info['mps_file']}', but received '{args.mps_file}'."
        )

    write_dd_file(aux_info, args.mps_file, args.output_dir)


if __name__ == "__main__":
    main()
