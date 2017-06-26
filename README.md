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

Prediction of binding affinities using all the docking poses from the unbound protein structures.

* **setup_submitter.py**: Script that permits to perform the docking of the structures and scoring of the obtained poses in the whole dataset.
  * Using the option 'docking', it executes ***setup_docking.py***, which runs Patchdock to obtain the docking poses (decoys), and runs Fiberdock to refine and rescore the decoys. Example:
    ```
    python /path/to/BADock/src/Binding_Affinity_Docking/setup_submitter.py -s docking
    ```
  * Using the option 'scoring', it executes ***setup_score.py***, which scores the decoys using the statistical potentials EPAIR, ES3DC and E3D. Example:
    ```
    python /path/to/BADock/src/Binding_Affinity_Docking/setup_submitter.py -s docking
    ```
* **setup_docking.py**: Runs Patchdock to obtain the docking poses and Fiberdock to refine and rescore them.

* **setup_score.py**: Scores the decoys using statistical potentials.

  ***Until here, it is not necessary to run the scripts because the results are already in BADock/workspace)***

* **score_nativecomplex.py**: Scores the native complexes of the dataset.

* **affinity_analysis.py**: Performs the prediction of binding affinities using the statistical potentials.

* **CCharPPI_affinityparse.py**: Performs the prediction of binding affinities using the CCharPPI server results.


#### Binding Mechanisms Docking

Analysis of the conformational space of the encounter complexes in four classes (Near-Native, Face-Face, Face-Back and Back-Back).

* **setup_submitter.py**: Script that permits to submit the different parts of the Binding Mechanism Docking analysis in a computer cluster.

* **setup_decoys.py**: Runs Patchdock to obtain the docking poses of the Docking Benchmark 5 and Fiberdock to refine and rescore them.

* **setup_decoys.py**: Runs Patchdock to obtain the docking poses of the Docking Benchmark 5 and Fiberdock to refine and rescore them.

* **setup_location.py**: Calculates the location (Near-Native, Face-Face, Face-Back and Back-Back) for each decoy of each protein in the Docking Benchmark 5.

* **setup_score.py**: Scores the decoys using statistical potentials.

  ***Until here, it is not necessary to run the scripts because the results are already in BADock/workspace)***

* **analyze_scores.py**: Get the scores of each decoy for each protein-protein interaction and analyze the different locations (Near-Native, Face-Face, Face-Back and Back-Back)
.
* **print_score_results.py**: Plot the boxplots of scores depending on the location (Near-Native, Face-Face, Face-Back and Back-Back) (Figure 1, Supplementary Figures S2, S3, S4).

* **plot_method_example.py**: Plot the figures of the examples of decoy classification by location (Supplementary Figure S1).

* **show_byprotein.py**: Print the score distributions of each location for each PPI.


#### BioLib

Python library used for the analyses of BADock.



### BADock/workspace

#### Binding Affinity Docking

* **ba_jobs**: Folder containing the Affinity results of all the protein-protein interactions in Affinity Benchmark 2.
* **img**: Folder containing the resulting images.
* **CCharPPI_AffinityBenchmark.csv**: File containing the scores from the CCharPPI server
* **scores_native.txt**: File containing the scores of the native decoys

#### Binding Mechanisms Docking

* **benchmark5_jobs**: Folder containing the Mechanisms results of all the protein-protein interactions in Docking Benchmark 5.
* **figures**: Folder containing the resulting images.
* **figure_method.csv**: File containing the images of the examples of location (from *plot_method_example.py*).
* **scores_locations_E+OX.txt**: File containing the scores of the protein-protein interactions



### BADock/revision

#### BADock/revision/workspace

This folder contains all the results asked in the revision and includes a report of these changes explaining them and how to obtain them.

* **revision_questions_report.pdf**: Report explaining the major changes in the review
* **Figures**: Folder containing the figures included in the review
* **Tables**: Folder containing the tables included in the review

#### BADock/revision/src

This folder contains all the source code to obtain the changes of the review and includes a README explaining how to execute the scripts. 

* **README**: Report explaining how to execute the scripts of the changes in the revision

