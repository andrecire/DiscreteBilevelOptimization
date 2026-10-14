"""
You have to install the MIPLIB Solution Checker
-> see readme-install-miplib-checker.md
Exemplary use of the feasibility checker for the instance moore90_2,
which is Example 2 of Moore and Bard 1990.
The test a solution given in website format
"""
from feasibility_checker import Feasibility_Checker
import json
import sys

mps_file = "./test_instances/moore90_2.mps.gz"
aux_file = "./test_instances/moore90_2.aux"
# solution is given in website format; see report
solution_file = "./test_solutions/moore90_2.res.json"
# load the solution
with open(solution_file, "r") as file:
    solution = json.load(file)
instance_name = solution["name"]

# check that the solution has a feasible point
if (
    solution["status"] == "infeasible"
    or solution["status"] == None
    or solution["status"] == "open"
):
    print("No feasible point in solution file")
else:
    # we now have a a feasible point
    solution_in_website_format = "True"
    feasibility_checker = Feasibility_Checker(
        mps_file,
        aux_file,
        solution_file,
        "./output_miplib_solution_checker",
        solution_in_website_format,
    )
    # check that high point relaxation is feasible
    if not feasibility_checker.hpr_infeasibility_checker():
        print(instance_name + " MILP relaxation is infeasible")
        sys.exit()  # we stop feasibility check
    else:
        print("MILP relaxation is feasible")

    # check that solutions satisfies constraints with MIPLIB solution checker
    if not feasibility_checker.validate_cons_feasibility_miplib2017_checker():
        print(
            instance_name
            + "solutions does not satisfies constraints according to MIPLIB solution checker"
        )
        sys.exit()  # we stop
    else:
        print("Solution passes MIPLIB solution checker")

    # we check that the point is bilevel feasible
    if not feasibility_checker.validate_bilevel_feasibility():
        print(instance_name + " solution is not bilevel feasible")
        sys.exit()  # we stop
    else:
        print("Lower-level solution is optimal")

    # check that the point is optimistic bilevel feasible
    if not feasibility_checker.check_optimistic_bilevel_feasible_point():
        print(instance_name + " solution is not optimistic bilevel feasible")
        sys.exit()
    else:
        print(
            "Solution of " + instance_name + " is optimistic bilevel feasible!"
        )
