# Network Relaxations for Combinatorial Bilevel Optimization

Code, computational results, and solution files accompanying the paper:

**Leonardo Lozano, David Bergman, and Andre A. Cire (2026). [Network Relaxations for Combinatorial Bilevel Optimization Under Linear Interactions](https://pubsonline.informs.org/doi/10.1287/opre.2025.2457). *Operations Research*. DOI: [10.1287/opre.2025.2457](https://doi.org/10.1287/opre.2025.2457).**

The implementation uses decision diagrams (DDs) to construct network relaxations for bilevel optimization. This repository includes the Java solvers, instance-conversion and feasibility-checking utilities, data for the computational analysis, and solutions obtained on the [BOBILib](https://bobilib.org/) benchmark.

## Repository contents

```text
.
├── README.md
├── solution_files/
│   ├── improved_boblib.csv
│   ├── improved_BOBLIB/          # 804 improved-solution JSON files
│   ├── optimal.csv
│   └── optimal/                  # 344 optimal-solution JSON files
└── sources/
    ├── code/
    │   ├── BilevelDDOR_DISCRETE/
    │   └── BilevelDDOR_MIXED/
    ├── data_examples/
    ├── results_tables_plots/
    │   ├── analysis.ipynb
    │   ├── benchmark_info.csv
    │   ├── collection-benchmark.csv
    │   ├── results_Oct19_2025_cont.csv
    │   └── results_solvers_plus_dd_Mar24_2026_mibs_bc_augmented.csv
    └── scripts/
        ├── feasibility_checker/
        ├── mps_aux_to_dd/
        └── mps_aux_to_legacy_aux/
```

All command examples below assume the repository root is the working directory, unless a `cd` command is shown.

## Solution files

The [solution_files](solution_files/) directory contains solutions from the October 19, 2025 experiments, using the DD solver with network width **8**. Each CSV indexes the JSON files in its corresponding subfolder.

| Folder | CSV index | Instances | Selection criterion |
| --- | --- | ---: | --- |
| [improved_BOBLIB](solution_files/improved_BOBLIB/) | [improved_boblib.csv](solution_files/improved_boblib.csv) | 804 | DD found a feasible objective at least 0.1 below the recorded BOBILib objective, or found a feasible solution when no objective was recorded. |
| [optimal](solution_files/optimal/) | [optimal.csv](solution_files/optimal.csv) | 344 | DD found a feasible solution and reported an absolute relative optimality gap below `1e-10`. |

Note: the two collections overlap, with 330 instances occurring in both. The optimal collection also includes instances for which DD proved optimality without improving the recorded objective. Comparisons use the BOBILib snapshot dated from October 19th, 2025.

Three additional entries were excluded from this public solution collection: `mine-166-5_50_50`, `ran14x18-disj-8_0_100`, and `weish23.dat-0` due to numerical issues. The original aggregate experimental CSV under `sources/` is retained for the paper's analysis.

### CSV columns

| Column | Meaning |
| --- | --- |
| `name` | BOBILib instance identifier. |
| `width` | DD network width, fixed at 8 for these collections. |
| `known_objective_value` | Objective recorded in the benchmark snapshot; blank when unavailable. |
| `LB` | Lower bound on the objective. When absent from the processed DD results, reconstructed as `UB - optgap * abs(UB)`. |
| `UB` | Objective of the feasible solution, matching `objective_value` in the JSON within numerical tolerance. |
| `optgap` | Relative optimality gap as a fraction; multiply by 100 to express it as a percentage. |

The bounds follow the minimization convention used in the processed results and can be negative. They are solver-reported numerical results; the JSON files under each subfolder contains the solution assignments.

### JSON format

Files are named `dd_<name>.dd_8-res.json`, where `<name>` is the CSV instance identifier. For example, `K5040W05.KNP` maps to `dd_K5040W05.KNP.dd_8-res.json`.

Each file contains:

- `name`: the input instance name, including the `.dd` suffix.
- `bilevel_type`: the bilevel convention used by the instance.
- `objective_value`: the objective of the exported solution.
- `upper_level_decisions` and `lower_level_decisions`: dictionaries mapping variable names to their values.
- `status` and `difficulty`: metadata fields that may be empty; use the CSV and collection membership to identify optimal results.

The original `.mps`/`.mps.gz` and `.aux` instances are available from [BOBILib](https://bobilib.org/). Representative examples are included in [sources/data_examples](sources/data_examples/).

Solution files were checked using BOBILib feasibility checker available at https://bobilib.org/misc.html.

## Java solvers

[sources/code](sources/code/) contains two Eclipse-style Java projects, with source files under `src/`:

- **`BilevelDDOR_DISCRETE`**: for binary problems with integer-valued data. This variant uses integer data structures and is recommended for instances satisfying these assumptions.
- **`BilevelDDOR_MIXED`**: accommodates general variable types and noninteger data using floating-point data structures.

Both variants use the same algorithm parameters described in Section 6 of the paper.

### Dependencies and compilation

The solvers require a Java development environment and IBM CPLEX/Concert. The experiments used **CPLEX 22.1.2**. Obtain CPLEX separately with an appropriate license.

Import the appropriate project into Eclipse, or compile it using your preferred Java build setup. Update the CPLEX JAR path in the project's `.classpath` to match your installation, and configure the CPLEX native-library path for your operating system. The included project settings contain paths from the development environment.

The main classes are:

- `BilevelBOBILIB_DISCRETE.Main`
- `BilevelBOBILIB_MIXED.Main`

The repository provides source code, not prebuilt solver JARs. If you export runnable JARs named `dd_solver_discrete.jar` and `dd_solver_mixed.jar`, configure their main classes and runtime dependencies before using the commands below.

### Running the solver

```bash
java -jar dd_solver_discrete.jar <fileName> <timeLimit> <numCores> <outputFolder> <algSwitch> <setWidth> <logicalCuts>
java -jar dd_solver_mixed.jar <fileName> <timeLimit> <numCores> <outputFolder> <algSwitch> <setWidth> <logicalCuts>
```

| Argument | Meaning |
| --- | --- |
| `fileName` | Path to an instance in the custom `.dd` format. |
| `timeLimit` | Time limit in seconds. |
| `numCores` | Number of solver threads. |
| `outputFolder` | Existing directory for the solution JSON. |
| `algSwitch` | Algorithm mode, described below. |
| `setWidth` | Network budget/width; `0` selects a heuristic default. |
| `logicalCuts` | Whether to enable the cuts from Section 5.5 (`true` or `false`). |

Algorithm modes:

- `0`: exact approach.
- `1`: root network relaxation; the discrete version also exports an MPS/AUX reformulation. Create `output_mps/` before requesting this export.
- `2`: high-point relaxation only.
- `3`: high-point relaxation augmented with the cuts from Section 5.5.

For example, after compiling and configuring the runtime:

```bash
mkdir -p output
java -jar dd_solver_discrete.jar path/to/instance.dd 7200 1 output 0 8 true
```

## Conversion and feasibility checking

The Python conversion utilities require **Gurobi's Python interface (`gurobipy`)** and a suitable Gurobi installation/license. CPLEX and Gurobi are separate dependencies and are not bundled with this repository.

### Convert MPS/AUX to DD

```bash
python3 sources/scripts/mps_aux_to_dd/mps_aux_to_dd.py \
  --mps_file path/to/instance.mps.gz \
  --aux_file path/to/instance.aux \
  --output_dir path/to/output_dir
```

This produces one `.dd` file for the Java solver. Uncompressed `.mps` input is also supported.

### Convert AUX to the legacy format

```bash
python3 sources/scripts/mps_aux_to_legacy_aux/mps_aux_to_legacy_aux.py \
  --mps_file path/to/instance.mps.gz \
  --aux_file path/to/instance.aux \
  --output_dir path/to/output_dir
```

This produces a `.laux` file.

### Check a solution

[sources/scripts/feasibility_checker/feasibility_checker_website](sources/scripts/feasibility_checker/feasibility_checker_website/) includes the BOBILib feasibility checker and an example driver for the JSON solution format. It requires Gurobi and the external MIPLIB 2017 solution checker.

Follow the [checker setup notes](sources/scripts/feasibility_checker/feasibility_checker_website/readme_feasibility_checker_website.md), configure the instance, solution, and output paths in the example driver, and run:

```bash
cd sources/scripts/feasibility_checker/feasibility_checker_website
python3 run_feasibility_checker_website.py
```

The checker tests constraint feasibility, integrality, and follower optimality for the supplied leader decisions.

## Computational results

[sources/results_tables_plots](sources/results_tables_plots/) contains the analysis notebook and supporting CSVs:

- `analysis.ipynb`: tables and figures for the computational study.
- `benchmark_info.csv`: the open BOBILib instances tested in the paper.
- `collection-benchmark.csv`: instance statistics and benchmark information.
- `results_Oct19_2025_cont.csv`: results for the main experiments in Section 6.
- `results_solvers_plus_dd_Mar24_2026_mibs_bc_augmented.csv`: experiments combining DD relaxations with other solvers, discussed in Appendix B.

To analyze the recorded results, use a Python environment with Jupyter, pandas, NumPy, Matplotlib, and seaborn. Run the notebook with `sources/results_tables_plots/` as the working directory, and create its figure-output directory first:

```bash
cd sources/results_tables_plots
mkdir -p figures
jupyter notebook analysis.ipynb
```

The CSVs support analysis of the recorded experiments without rerunning the optimization solvers. Repeating the solver experiments requires the original benchmark instances and the solver dependencies described above.

Please feel free to reach out to us with any questions, concerns, or errors you notice in the code or results; we would be happy to hear from you, help, and make corrections.

## Citation

Please cite the paper when using this code or these results:

```bibtex
@article{lozano2026network,
  author  = {Lozano, Leonardo and Bergman, David and Cire, Andre A.},
  title   = {Network Relaxations for Combinatorial Bilevel Optimization Under Linear Interactions},
  journal = {Operations Research},
  year    = {2026},
  doi     = {10.1287/opre.2025.2457},
  url     = {https://doi.org/10.1287/opre.2025.2457}
}
```
