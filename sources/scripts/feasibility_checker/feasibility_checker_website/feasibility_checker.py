"""
input:
    - path_to_instance mps file
    - path to instance aux file
    - path to solution file (None if not existing)
    - path to directory in which we store the sol format and logfile
      of the MIPLIB solution checker
    - bool if solutions are from mibs/filmosi or in website format (--sol_format_website)
methods:
    - hpr_infeasibility_checker -> solves high point relaxation
    - validate_cons_feasibility -> checks that for given point
      of solution file the constraints (upper and lower level)
      are satisfied using the MIPLIB2017 solution checker
    - validate_bilevel_feasibility
      - first solve for fixed upper
        level w.r.t. given solution of solution file
        the lower-level problem to optimality
        and compare the given and computed lower-level objective value
      - second solve the high point relaxation with the additional
        constraint lower-level obj <= optimal lower-level obj value of given solution.
        Then, check that the leader objective value is equivalent to the
        objective value of the solution. This guarantees that we have an
        optimistic bilevel feasible point.
"""

import gurobipy as gp
import json
import os
import subprocess
from name_based_bilevel_parser.name_based_bilevel_parser import (
    Name_Based_Bilevel_Parser,
)


class Feasibility_Checker:
    def __init__(
        self,
        path_mps_file,
        path_aux_file,
        path_solution_file,
        store_dir_sol,
        bool_sol_format_website,
    ):
        self._path_mps_file = path_mps_file
        self._path_aux_file = path_aux_file
        self._path_sol_file = path_solution_file
        self._instance_name = path_mps_file.rsplit("/", 1)[-1].split(".mps")[0]
        self._store_dir_sol_MIPLIB_Checker = store_dir_sol
        # parser contains bilevele info, e.g., lower-level var_names etc.
        self._name_based_bilevel_parser = Name_Based_Bilevel_Parser(
            path_mps_file, path_aux_file
        )
        # string "True" if solution is given in website format,
        # "False" if we use our internal format
        self._sol_format_website = bool_sol_format_website

        # tolerance for constraint violations in gurobi (default 1e-6)
        self._cons_tolerance = 1e-6
        # tolerance for lower-level optimality (used in own verification; not in gurobi)
        self._optimality_tolerance = 1e-4

        # MIPLIB feasibility_checker tolerances default is 1e-4
        # same applies also to objective tolerance -> see output.log
        self._MIPLIB_SOLChecker_linear_cons_tolerance = 1e-4
        self._MIPLIB_SOLChecker_linear_integer_tolerance = 1e-4

        # number of threads used by gurobi
        self._nr_threads = 4

        # lower-level objective_value of sol file,
        # will be computed in validate_bilevel_feasibility()
        self._lower_level_obj_sol = None

    def hpr_infeasibility_checker(self):
        """
        Solve the high point relaxation.
        If feasible point exists return True else False
        """
        # create empty gurobi model without unncessary log output
        with gp.Env(empty=True) as env:
            env.setParam("OutputFlag", 0)
            env.start()
            with gp.Model(env=env) as tmp_model:
                tmp_model = gp.read(self._path_mps_file, env)

        # solve model and stop after first feasible point
        tmp_model.Params.SOLUTION_LIMIT = 1
        tmp_model.Params.Threads = self._nr_threads
        tmp_model.optimize()

        # check if we have a feasible point
        status = tmp_model.status
        if (
            status == gp.GRB.Status.OPTIMAL
            or status == gp.GRB.Status.SOLUTION_LIMIT
        ):
            return True
        elif status == gp.GRB.Status.INFEASIBLE:
            print(
                "Instance {}: High Point Relaxation is infeasible".format(
                    self._instance_name
                )
            )
            return False
        else:
            raise Exception(
                "No Suitable Gurobi Status {} when solving high point relaxation".format(
                    status
                )
            )

        return

    def validate_cons_feasibility_miplib2017_checker(self):
        """
        We check that the solution satisfies the constraints and the stored
        objective value equals the value of the stored solution using the
        MIPLIB2017 solution checker.
        Benefit: check feasibility of solution without solver
        1) take our solution and write .sol MIPLIB file see README in
           MIPLIB2017 solution checker
           we write this (temporary file to an own folder
           in the feasibility checker directory
        2) Check feasibility with MIPLIB2017 feasibility checker
        3) analyse result of MIPLIB2017 checker
        Note we assume that the full solution including variables with
        value zero are given in our solution dictionary
        """
        # now read solution
        with open(self._path_sol_file) as file:
            dic_sol = json.load(file)

        # check that solution exists
        if self._sol_format_website == "True":
            if dic_sol["status"] == None or dic_sol["status"] == "infeasible":
                raise Exception(
                    "Instance {}: solution file does not contain a solution".format(
                        self._instance_name
                    )
                )
        else:
            assert self._sol_format_website == "False"
            if (
                dic_sol["solver_status"] == None
                or dic_sol["solver_status"] == "infeasible"
                or dic_sol["feasible"] == -1
            ):
                raise Exception(
                    "Instance {}: solution file does not contain a solution".format(
                        self._instance_name
                    )
                )

        # create temporary directory to store the solution in MIPLIB sol format
        if not os.path.exists(self._store_dir_sol_MIPLIB_Checker):
            os.makedirs(self._store_dir_sol_MIPLIB_Checker)
        # now create the MIPLIB sol file of our bilevel solution
        # create new sol file
        with open(
            self._store_dir_sol_MIPLIB_Checker
            + "/"
            + self._instance_name
            + ".sol",
            "w",
        ) as file:
            # write first the objective value to the file
            file.write("=obj=" + " " + str(dic_sol["objective_value"]) + "\n")
            # now write line by line upper level varname space value
            for varName in dic_sol["upper_level_decisions"].keys():
                file.write(
                    varName
                    + " "
                    + str(dic_sol["upper_level_decisions"][varName])
                    + "\n"
                )
            # now write line by line lower level varname space value
            for varName in dic_sol["lower_level_decisions"].keys():
                file.write(
                    varName
                    + " "
                    + str(dic_sol["lower_level_decisions"][varName])
                    + "\n"
                )
        # now we have a MIBLIB solution file
        # call the MIBLIB solution checker for the high point relaxation and our solution
        # use subprocess and the command ./solchecker modelfile solutionfile [lin_tol int_tol]
        # store output of MIPLIB solution checker to verify the result in the same dir as .sol
        with open(
            self._store_dir_sol_MIPLIB_Checker
            + "/"
            + self._instance_name
            + "_MIPLIB_SOLChecker.log",
            "w",
        ) as output:
            command = [
                "/miplib2017-testscript-v1.0.3/checker/bin/solchecker",
                self._path_mps_file,
                self._store_dir_sol_MIPLIB_Checker
                + "/"
                + self._instance_name
                + ".sol",
                "["
                + str(self._MIPLIB_SOLChecker_linear_cons_tolerance)
                + " "
                + str(self._MIPLIB_SOLChecker_linear_integer_tolerance)
                + "]",
            ]
            try:
                subprocess.run(command, stdout=output, timeout=3800)
                # run waits until the process is finished
            except subprocess.TimeoutExpired:
                print(
                    "Set timelimit was not respected by solver. Killed process on my own!"
                )
        # now we evaluate the output of the MIPLIB solution checker
        # read output
        check_sol_in_file = False
        with open(output.name, "r", encoding="UTF-8") as file:
            for line in file:
                if "Check SOL" in line:
                    check_sol_in_file = True
                    # we have the line: Check SOL: Integrality 1 Constraints 1 Objective 1
                    # where 1 indicates we pass the test and 0 note
                    # check that we pass all test
                    if "Integrality 1 Constraints 1 Objective 1" not in line:
                        # some feasibility tests failed
                        return False

        if check_sol_in_file:
            return True
        else:
            print(
                "MIPLIB Feasibility Checker Output does not include Check Sol for "
                + self._instance_name
            )
            return False

    def validate_bilevel_feasibility(self):
        """
        validate bilevel feasibility by fixing upper-level solution of
        self._path_sol_file and then solve lower-level problem to optimality.
        Compare objective value of this solved lower level with the
        objective value of the lower level if the solution of
        self._path_sol_file is fixed.
        """
        # now read solution
        with open(self._path_sol_file) as file:
            dic_sol = json.load(file)
        # check that solution exists
        if self._sol_format_website == "True":
            if dic_sol["status"] == None or dic_sol["status"] == "infeasible":
                raise Exception(
                    "Instance {}: solution file does not contain a solution".format(
                        self._instance_name
                    )
                )
        else:
            assert self._sol_format_website == "False"
            if (
                dic_sol["solver_status"] == None
                or dic_sol["solver_status"] == "infeasible"
                or dic_sol["feasible"] == -1
            ):
                raise Exception(
                    "Instance {}: solution file does not contain a solution".format(
                        self._instance_name
                    )
                )

        # read high point relaxation
        # create empty gurobi model without unncessary log output
        with gp.Env(empty=True) as env:
            env.setParam("OutputFlag", 0)
            env.start()
            with gp.Model(env=env) as model:
                model = gp.read(self._path_mps_file, env)

        # create dictionary key: var_name, value 0 only for upper level variables
        dic_sol_values = {}
        for var_name in self._name_based_bilevel_parser._upper_level_var_names:
            dic_sol_values[var_name] = 0.0
        # now we update the dictionary with the upper-level variables values
        dic_sol_values.update(dic_sol["upper_level_decisions"])
        # note that if a variable is not in upper/lower_level_decisions then it will be set to zero
        # this is fine since in the solution variables with zero value are not explicitly reported

        # fix all upper level variables
        for var_name in self._name_based_bilevel_parser._upper_level_var_names:
            var = model.getVarByName(var_name)
            if var_name in dic_sol["upper_level_decisions"].keys():
                var.lb = dic_sol["upper_level_decisions"][var_name]
                var.ub = dic_sol["upper_level_decisions"][var_name]
            else:
                # if not mentioned value is zero
                var.lb = 0.0
                var.ub = 0.0

        # delete upper level constraints
        for (
            cons_name
        ) in self._name_based_bilevel_parser._upper_level_cons_names:
            cons = model.getConstrByName(cons_name)
            model.remove(cons)
        # replace objective by lower-level objective function
        obj_expr = gp.LinExpr()
        for (
            var_name,
            coef,
        ) in self._name_based_bilevel_parser._lower_level_objective_coefs:
            var = model.getVarByName(var_name)
            obj_expr += var * coef
        model.setObjective(obj_expr, sense=gp.GRB.MINIMIZE)
        # model update
        model.update()

        # set constraint feasibility tolerance
        model.Params.FeasibilityTol = self._cons_tolerance

        # now model equals lower-level problem for fixed upper level decisions
        # solve model
        model.optimize()

        # check that problem is solved to optimality
        if model.status == gp.GRB.OPTIMAL:
            # store objective value
            obj_value_ll = model.objVal
        elif model.status == gp.GRB.INFEASIBLE:
            print(
                "Instance {}: Lower-level is infeasible w.r.t. fixed upper-level decisions".format(
                    self._instance_name
                )
            )
            return False
        else:
            raise Exception(
                "Gurobi Solver Status {} not known of instance {}".format(
                    model.status, self._instance_name
                )
            )

        # now we compute the lower-level objective of the given solution by hand
        obj_value_ll_sol = 0.0
        for (
            var_name,
            coef,
        ) in self._name_based_bilevel_parser._lower_level_objective_coefs:
            if var_name in dic_sol["lower_level_decisions"].keys():
                var_value = dic_sol["lower_level_decisions"][var_name]
            else:
                var_value = 0.0
            # add value to objective_value
            obj_value_ll_sol += var_value * coef
        # store the lower-level objective value for the
        # optimistic bilevel feasibility checker later
        self._lower_level_obj_sol = obj_value_ll_sol

        # now compare objective values of solving lower-level problem of the given point and
        # computed by gurobi with fixed upper-level
        # we adapt cplex gap |UB-LB|/(10**-10 + |LB|). Thus, we assume min-min problems, which
        # leads to |obj_value_ll_sol - obj_value_ll|/(10^-10 + obj_value_ll)
        if (
            abs((obj_value_ll_sol - obj_value_ll))
            / (10 ** (-10) + abs(obj_value_ll))
            > self._optimality_tolerance
        ):
            print(
                "Solution of Instance {} not bilevel feasible".format(
                    self._instance_name
                )
                + "LL Obj sol {} versus LL Obj fixed upper-level decisions {}".format(
                    obj_value_ll_sol, obj_value_ll
                )
            )
            return False
        else:
            return True

        return

    def check_optimistic_bilevel_feasible_point(self):
        """
        Requires that lower-level objective_value of the solution file
        is computed before, which is done in
        validate_bilevel_feasibility() -> self._lower_level_obj_sol
        Add to the high point relaxation the constraints
        lower-level objective <= lower_level_obj_value_sol
        Note assume that we have min-min problems as in the bilevel library
        Fix upper level decision and compute best follower's decision
        for the given upper-level decision
        """
        # read the hpr model
        # create empty gurobi model without unncessary log output
        with gp.Env(empty=True) as env:
            env.setParam("OutputFlag", 0)
            env.start()
            with gp.Model(env=env) as model:
                model = gp.read(self._path_mps_file, env)

        # now add the constraints ll_obj_func <= ll_obj_value_of_solution
        # check that the lower-level objective of the solution file was already computed
        # this takes place in self._validate_bilevel_feasibility()
        assert self._lower_level_obj_sol != None
        cons_name = "lower_level_obj_value"
        # build lower-level objective
        # replace objective by lower-level objective function
        obj_expr = gp.LinExpr()
        for (
            var_name,
            coef,
        ) in self._name_based_bilevel_parser._lower_level_objective_coefs:
            var = model.getVarByName(var_name)
            obj_expr += var * coef
        # set constraints
        cons = obj_expr <= self._lower_level_obj_sol
        model.addConstr(cons, name=cons_name)

        # now fix upper level decisions
        # read solution
        with open(self._path_sol_file) as file:
            dic_sol = json.load(file)

        # check that solution exists
        if self._sol_format_website == "True":
            if dic_sol["status"] == None or dic_sol["status"] == "infeasible":
                raise Exception(
                    "Instance {}: solution file does not contain a solution".format(
                        self._instance_name
                    )
                )
            obj_leader_sol = dic_sol["objective_value"]
        else:
            assert self._sol_format_website == "False"
            if (
                dic_sol["solver_status"] == None
                or dic_sol["solver_status"] == "infeasible"
                or dic_sol["feasible"] == -1
            ):
                raise Exception(
                    "Instance {}: solution file does not contain a solution".format(
                        self._instance_name
                    )
                )
            # store objective value leader for later
            obj_leader_sol = dic_sol["objective_value"]

        # create dictionary key: var_name, value 0 only for upper level variables
        dic_sol_values = {}
        for var_name in self._name_based_bilevel_parser._upper_level_var_names:
            dic_sol_values[var_name] = 0.0
        # now we update the dictionary with the upper-level variables values
        dic_sol_values.update(dic_sol["upper_level_decisions"])
        # note that if a variable is not in upper/lower_level_decisions then it will be set to zero
        # this is fine since in the solution variables with zero value are not explicitly reported

        # fix all upper level variables
        for var_name in self._name_based_bilevel_parser._upper_level_var_names:
            var = model.getVarByName(var_name)
            if var_name in dic_sol["upper_level_decisions"].keys():
                var.lb = dic_sol["upper_level_decisions"][var_name]
                var.ub = dic_sol["upper_level_decisions"][var_name]
            else:
                # if not mentioned value is zero
                var.lb = 0.0
                var.ub = 0.0

        model.update()

        # set constraint feasibility tolerance
        model.Params.FeasibilityTol = self._cons_tolerance
        # set dual reduction zero
        model.Params.DualReductions = 0
        # set number of threads
        model.Params.Threads = self._nr_threads

        model.optimize()
        # check that problem is solved to optimality
        if model.status == gp.GRB.OPTIMAL:
            # store objective value
            obj_leader = model.objVal
        else:
            print(
                "Gurobi Solver Status {} not optimal of instance {}".format(
                    model.status, self._instance_name
                )
            )
            return False
        # compare objective leader to objective of solution file
        # we adapt cplex gap |UB-LB|/(10**-10 + |LB|). Thus, we assume min-min problems, which
        # leads to |obj_leader_sol - obj_leader|/(10^-10 + obj_leader)
        if (
            abs((obj_leader_sol - obj_leader))
            / (10 ** (-10) + abs(obj_leader))
            > self._optimality_tolerance
        ):
            print(
                "Solution of Instance {} not optimistic bilevel feasible".format(
                    self._instance_name
                )
                + "Obj sol {} versus Obj HPR with ll_obj_cons {}".format(
                    obj_leader_sol, obj_leader
                )
            )
            return False
        else:
            return True

        return
