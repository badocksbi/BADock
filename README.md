# BADock

Prediction of protein interactions binding affinity from unbound tertiary structures

## Getting Started

### Summary

* **db**: Contains the two datasets used in the analysis, Affinity Benchmark 2 and Docking Benchmark 5
* **src**: Contains the source code
* **workspace**: Contains the results of the analysis
* **revision**: Contains the analysis of the revision

### src

#### Binding Affinity Docking

Prediction of binding affinities using all the docking poses from the unbound protein structures.

* **setup_submitter.py**: Script that permits to perform the docking of the structures and scoring of the obtained poses in the whole dataset.<br />
  * Using the option 'docking', it executes ***setup_docking.py***, which runs Patchdock to obtain the docking poses (decoys), and runs Fiberdock to refine and rescore the decoys. Example:
    ```
    python /path/to/BADock/src/Binding_Affinity_Docking/setup_submitter.py -s docking -m
    ```
  * Using the option 'scoring', it executes ***setup_score.py***, which scores the decoys using the statistical potentials EPAIR, ES3DC and E3D. Example:
    ```
    python /path/to/BADock/src/Binding_Affinity_Docking/setup_submitter.py -s docking -m
    ```
* **setup_docking.py**: Runs Patchdock to obtain the docking poses and Fiberdock to refine and rescore them
* **setup_score.py**: Scores the decoys using statistical potentials

