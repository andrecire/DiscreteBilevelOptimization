**Install MIPLIB Solution Checker**

Download the MIPLIB Solution Checker: miplib2017-testscript-v1.0.3.zip into the directory of the feasibility checker.
see https://miplib.zib.de/download.html

1. unzip the install sources (already included in the git) by
   `unzip miplib2017-testscript-v1.0.3.zip -d ./miplib2017-testscript-v1.0.3`

2. Install the MIBLIB solution checker following the readme in miplib2017-testscript-v1.0.3
   only do the first step, i.e.,
   `make checker`

**Bilevel Feasibility Checker**

The single steps of the feasibility checker are described in the report.
The feasibility checker should be used with solutions in the website format; see the report.

(Note that the feasibility_checker.py also supports internal solution formats for MibS and Filmosi)

**Input**:
    - path_to_instance mps file
    - path to instance aux file
    - path to solution file (None if not existing)
    - path to directory in which we store the sol format and logfile
      of the MIPLIB solution checker

**Tolerances**:
    - self._cons_tolerance = 1e-6 (in gurobi)
    - self._optimality_tolerance = 1e-4
    - self._MIPLIB_SOLChecker_linear_cons_tolerance = 1e-4 (precision of checking linear constraints with MIPLIB checker)
    - self._MIPLIB_SOLChecker_linear_integer_tolerance = 1e-4 (precision of checking integer conditions with MIPLIB checker)

**Parameters**
    - self._nr_threads (solved by gurobi; default 4)

All of the following methods have to be successfully conducted to verify bilevel feasibility.
If one fails, the instance does not pass the feasibility check and we terminate the programm.
methods:
    - hpr_infeasibility_checker -> solves high point relaxation
    - validate_cons_feasibility_miplib2017_checker -> checks that for given point
      of solution file the constraints (upper and lower level)
      are satisfied by MIPLIB 2017 checker
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

**running the feasiblity checker**
`python run_feasibility_checker_website.py`
