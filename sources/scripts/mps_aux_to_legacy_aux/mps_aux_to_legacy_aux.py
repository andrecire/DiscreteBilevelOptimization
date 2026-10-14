# ---------------------------------------------
# Create legacy AUX files
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
            aux_info["mps_file"] = lines[i].split("/")[-1]  # Store just the filename, not the path
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


def write_mibs_legacy_aux_file(name, aux_info, mps_file_path, output_folder):
    """
    Write a legacy-style .laux file in MIBS format, using Gurobi to get correct index order
    for lower-level variables and constraints as per the MPS file.
    """
    import gurobipy as gp

    output_dir = output_folder
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, f"{name}.laux")

    # Load MPS to get the ordering of variables and constraints
    print(f"Processing {name}: loading MPS from {mps_file_path}")
    model = gp.read(mps_file_path)

    # Create lookup for variable and constraint indices
    var_name_to_index = {v.VarName: i for i, v in enumerate(model.getVars())}
    constr_name_to_index = {c.ConstrName: i for i, c in enumerate(model.getConstrs())}

    # Index lists for lower-level variables and constraints
    lower_var_indices = [var_name_to_index[vname] for vname, _ in aux_info["lower_level_variables"]]
    lower_constr_indices = [constr_name_to_index[cname] for cname in aux_info["lower_level_constraints"]]

    # Sort indices in MPS order
    sorted_var_indices = sorted(lower_var_indices)
    sorted_constr_indices = sorted(lower_constr_indices)

    # Extract objective coefficients in the same order
    varname_to_coeff = dict(aux_info["lower_level_variables"])
    sorted_coeffs = [varname_to_coeff[model.getVars()[i].VarName] for i in sorted_var_indices]

    # Write the .laux file
    with open(output_path, "w") as f:
        f.write(f"N {aux_info['num_lower_vars']}\n")
        f.write(f"M {aux_info['num_lower_constraints']}\n")

        for idx in sorted_var_indices:
            f.write(f"LC {idx}\n")
        for idx in sorted_constr_indices:
            f.write(f"LR {idx}\n")
        for coeff in sorted_coeffs:
            # coeff_int = int(coeff)
            # f.write(f"LO {coeff_int}\n")
            f.write(f"LO {coeff}\n")

        f.write("OS 1\n")  # Always minimization by convention

    print(f"Wrote Gurobi-indexed .laux to: {output_path}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert one bilevel instance from .mps/.mps.gz + .aux into legacy .laux format."
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
        help="Directory where the output .laux file should be created.",
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

    write_mibs_legacy_aux_file(
        infer_instance_name(args.mps_file),
        aux_info,
        args.mps_file,
        args.output_dir,
    )


if __name__ == "__main__":
    main()
