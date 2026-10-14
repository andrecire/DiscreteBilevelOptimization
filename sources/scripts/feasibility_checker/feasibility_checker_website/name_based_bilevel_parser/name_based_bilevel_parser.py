"""
Input:
    - mps file
    - aux file
Creates a class object that contains the following attributes:
    - list of upper-level variable names
    - list of upper-level constraint names
    - list of lower-level variable names
    - list of lower-level constraint names
"""
import gurobipy as gp


class Name_Based_Bilevel_Parser:
    def __init__(self, path_mps_file, path_aux_file):
        self._path_mps_file = path_mps_file
        self._path_aux_file = path_aux_file
        self._instance_name = path_mps_file.rsplit("/", 1)[-1].split(".mps")[0]

        self._upper_level_var_names = []
        # list of tuples (var_name, lower_bound)
        self._upper_level_var_lb = []
        # list of tuples (var_name, upper_bound)
        self._upper_level_var_ub = []
        self._upper_level_cons_names = []
        self._lower_level_var_names = []
        # list of tuples (var_name, lower_bound)
        self._lower_level_var_lb = []
        # list of tuples (var_name, upper_bound)
        self._lower_level_var_ub = []
        self._lower_level_cons_names = []
        # list of tuples (var_name, objective_coefficient)
        self._lower_level_objective_coefs = []

        self._upper_level_nr_vars = None
        self._upper_level_nr_cons = None
        self._lower_level_nr_vars = None
        self._lower_level_nr_cons = None

        self.parse_auxfile()
        self.determine_upper_level_cons_vars()
        self.determine_variable_bounds()

    def parse_auxfile(self):
        """
        parse auxfile and determine lower-level variables and constraints
        """
        # bool indicating that variable section has began
        var_section = False
        # bool indicating that constraint section has began
        cons_section = False

        with open(self._path_aux_file, "r") as file:
            lines = file.readlines()
            # check that first line contains @NUMVARS to ensure that we have name_based format
            assert "@NUMVARS" in lines[0]
            # now we process each line
            for index, line in enumerate(lines):
                if "@NUMVARS" in line:
                    # next line contains number of lower-level variables
                    self._lower_level_nr_vars = float(lines[index + 1].strip())
                elif "@NUMCONSTRS" in line:
                    # next line contains number of lower-level constraints
                    self._lower_level_nr_cons = float(lines[index + 1].strip())
                elif "@VARSBEGIN" in line:
                    # the following lines contain the lower-level variable names
                    var_section = True
                elif "@VARSEND" in line:
                    # the end of the variable section is reached
                    var_section = False
                elif var_section:
                    # we are in the variable section, i.e., line contains a lower-level variable name
                    list_varname_objcoef = line.strip().split()
                    varname = list_varname_objcoef[0]
                    objcoef = float(list_varname_objcoef[1])
                    self._lower_level_var_names.append(varname)
                    self._lower_level_objective_coefs.append(
                        (varname, objcoef)
                    )
                elif "@CONSTRSBEGIN" in line:
                    # the following lines contain the lower-level constraint names
                    cons_section = True
                elif "@CONSTRSEND" in line:
                    # the end of the constraint section
                    cons_section = False
                elif cons_section:
                    # we are in the constraint section, i.e., the line constains a lower-level constraint name
                    self._lower_level_cons_names.append(line.strip())

    def determine_upper_level_cons_vars(self):
        """
        determine upper-level variables and constraints by reading MILP relaxation
        excluding/skipping previously determined lower-level variables and constraints
        """

        # read mps file, i.e., high point relaxation
        # create empty gurobi model without unncessary log output
        with gp.Env(empty=True) as env:
            env.setParam("OutputFlag", 0)
            env.start()
            with gp.Model(env=env) as hpr_model:
                hpr_model = gp.read(self._path_mps_file, env)
        # check that we have determined lower level variables
        assert self._lower_level_var_names != []
        # now iterate over all variables
        for var in hpr_model.getVars():
            varname = var.VarName
            if varname not in self._lower_level_var_names:
                self._upper_level_var_names.append(var.VarName)

        # now iterate over all constraints
        for cons in hpr_model.getConstrs():
            consname = cons.ConstrName
            # check that we have determined lower level constraints
            # note this set can be empty if we have only bounds
            if consname not in self._lower_level_cons_names:
                self._upper_level_cons_names.append(consname)

        # determine number of upper-level variables and constraints
        self._upper_level_nr_vars = len(self._upper_level_var_names)
        assert (
            self._upper_level_nr_vars
            == len(hpr_model.getVars()) - self._lower_level_nr_vars
        )
        self._upper_level_nr_cons = len(self._upper_level_cons_names)
        assert (
            self._upper_level_nr_cons
            == len(hpr_model.getConstrs()) - self._lower_level_nr_cons
        )

    def determine_variable_bounds(self):
        """
        Store variable lower and upper bounds in list of the format
        [(var_name, lb), ...], ....
        """

        # read the high point relaxation
        # read mps file, i.e., high point relaxation
        # create empty gurobi model without unncessary log output
        with gp.Env(empty=True) as env:
            env.setParam("OutputFlag", 0)
            env.start()
            with gp.Model(env=env) as hpr_model:
                hpr_model = gp.read(self._path_mps_file, env)

        # start with lower-level variables
        assert self._lower_level_var_names != []
        for var_name in self._lower_level_var_names:
            # get the variable by name
            var = hpr_model.getVarByName(var_name)
            self._lower_level_var_lb.append((var_name, var.lb))
            self._lower_level_var_ub.append((var_name, var.ub))

        # upper level variables
        for var_name in self._upper_level_var_names:
            self._upper_level_var_lb.append((var_name, var.lb))
            self._upper_level_var_ub.append((var_name, var.ub))
