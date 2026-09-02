# gtfsOutliers

> **Generalized Trimmed Functional Score (GTFS) for functional shape outlier detection**

This repository provides the official R package implementation and simulation toolkits for the **GTFS** algorithm, designed to detect shape outliers in functional data.

---

## Repository Structure

The project is structured as a standard R package environment coupled with workflow and reproducibility scripts:

### Core Package Functions (`R/`)
*   **`GTFS_main.R`**: Core implementation of the proposed GTFS algorithms and C-statistic calculations.
*   **`LTFS_calculation.R`**: Wrapper for the benchmark LTFS pipeline.
*   **`add_intercept_add_centering.R`**: Preprocessing utilities for handling data centering and adaptive random intercepts.
*   **`plotting_utils.R`**: Visualization suites for functional curves contaminated with outliers.
*   **`simulation_generators.R`**: Synthetic functional shape outlier data generation engines (`simul_1` to `simul_6`).

### Workflow & Execution Scripts (`scripts/`)
*   **`GTFS_analysis_example.R`**: Step-by-step tutorial demonstrating GTFS analysis pipelines.
*   **`data_generation_plotting_centering_smoke_test.R`**: Comprehensive pipeline validation and graphing test runner.
*   **`main_simulation.R`**: Full-scale multi-iteration simulation script benchmarking empirical power.
*   **`Real_data_human_activity_analysis.R`**: Real data application on human activity recognition data.
*   **`Real_data_tecator_analysis.R`**: Real data application on Tecator data.

---

## Installation & Quick Start

You can install this development version directly from GitHub using the `devtools` package:

```R
# In your R console
if (!require("devtools")) install.packages("devtools")
devtools::install_github("HyungJun94/gtfs-functional-outliers")
```

## 📜 References

*   **LTFS**: Ren, H., Chen, N., & Zou, C. (2017). Projection-based outlier detection in functional data. *Biometrika*, 104(2), 411-423.
*   **MDP/ReMDP**: Ro, K., Zou, C., Wang, Z., & Yin, G. (2015). Outlier detection for high-dimensional data. *Biometrika*, 102(3), 589–599.
