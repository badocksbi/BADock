# BADock

Prediction of protein interactions binding affinity from unbound tertiary structures

## Getting Started

### Summary

* **BADock/db**: Contains the two datasets used in the analysis, Affinity Benchmark 2 and Docking Benchmark 5
* **BADock/src**: Contains the source code
* **BADock/workspace**: Contains the results of the analysis
* **BADock/revision**: Contains the analysis of the revision

### BADock/src

#### Binding Affinity Docking

Prediction of binding affinities using all the docking poses from the unbound protein structures <span style="color:blue">*(not necessary, the results are already in BADock/workspace)*</span>.

* **setup_submitter.py**: Script that permits to perform the docking of the structures and scoring of the obtained poses in the whole dataset.
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
* **score_nativecomplex.py**: Scores the native complexes of the dataset
* **score_nativecomplex.py**: Scores the native complexes of the dataset
* **affinity_analysis.py**: Performs the prediction of binding affinities using the statistical potentials
* **CCharPPI_affinityparse.py**: Performs the prediction of binding affinities using the CCharPPI server results


#### Binding Mechanism Docking

Analysis of the conformational space of the encounter complexes in four classes (Near-Native, Face-Face, Face-Back and Back-Back).

