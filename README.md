# Using MDAnalysis for Efficient Simulation Pre- and Post-Processing [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/MDAnalysis/WorkshopPrace2021/HEAD)

## Contents

This repository contains the materials for the 3-day [MDAnalysis 2021 PRACE/SURF
workshop](https://www.mdanalysis.org/2021/04/09/prace-workshop/).


The workshop is structured in the following manner:

### Day 1:

1. `Day1-Session1-Lecture`
  * `Welcome.ipynb`

    Workshop welcome presentation.

  * `Lecture1_Molecules.ipynb`

    The first lecture covering the general structure, objects and workflows
    of MDAnalysis.


2. `Day1-Session1-Practical`
  * `session1_practical.ipynb`

    Practical complementing `Day1-Session1-Lecture/Lecture1_Molecules.ipynb`
    introducing the fundamental objects of MDAnalysis (e.g. `Universe`,
    `AtomGroup`, `Atom`, `Residues`, `Segments`).


3. `Day1-Session2-Lecture`
  * `session2-dynamics.ipynb`

    Lecture covering the main concepts of coordinate I/O in MDAnalysis.


4. `Day1-Session2-Practical`
  * `session2.ipynb`

    Practical complementing `Day1-Session2-Lecture/session2-dynamics.ipynb`.
    This practical demonstrates the main objects linked to coordinate I/O,
    and demonstrates how they can be used for basic trajectory analysis.


### Day 2:

1. `Day2-Session1-Lecture`
  * `lib.distances.pdf`

    Lecture covering both the `MDAnalysis.lib.distances` module and the
    MDAnalysis' analysis framework (e.g. `AnalysisBase`).


2. `Day2-Session1-Practical`
  * `day2session1-lib.distances.ipynb`

    Practical complementing `Day2-Session1-Lecture/lib.distances.pdf`.
    Demonstrates how to calculate hydrogen bonds using both
    `MDAnalysis.lib.distances` and `MDAnalysis.analysis.hydrogenbonds`.


3. `Day2-Session2-Lecture`
  * `universe_creation_and_manipulation-lecture.ipynb`

    Lecture covering how to;
      * Create, copy, and merge `Universe` objects
      * Manipulating positions (e.g. translate, rotate, PBC wrap/unwrap)
      * Talk to other libraries through `MDAnalysis.converters`


4. `Day2-Session2-Practical`
  * `day2session2-Analysis.ipynb`

    Follow-up to `Day2-Session1-Lecture/lib.distances.pdf` and
    `Day2-Session1-Practical/day2session1-lib.distances.ipynb`. This practical
    demonstrates how one could write their own `AnalysisBase`-derived analysis
    class.


### Day 3:

1. `Day3-Session1-Lecture`
  * `MDAParallelization.pdf`

    Lecture convering various strategies for parallelizing MDAnalysis-based
    workflows.


2. `Day3-Session1-Practical`
  * `parallelism.ipynb`

    Practical complementing `Day3-Session1-Lecture/MDAParallelization.pdf`.
    Here we demonstrate how both multiprocessing and MPI-based parallelism
    can be employed to analyze lipid headgroup orientation.
    **Note: the MPI portion of this practical assumes access to SURF's LISA HPC platform, see `Day3-Session1-Practical/README.md` for more information.**


3. `Day3-Session2-Lecture`
  * `more-features.ipynb`

    Final lecture demonstrating some of the other features available in
    MDAnalysis such as auxiliary data reading, and on-the-fly transformations.


4. `Day3-Session2-Practical`
  * `creating_and_manipulating_universes.ipynb`

    Practical complementing `Day2-Session2-Lecture/universe_creation_and_manipulation-lecture.ipynb`
    This practical demonstrates how to seemlessly interface [RDKit](https://github.com/rdkit/rdkit)
    and [OpenMM](https://github.com/openmm/openmm) with MDANalysis.


5. `Day3-Closing-Lecture`
  * `MDAnalysis_Workshop_final_no_logo.pdf`

    SURF / PRACE workshop closing remarks.



## Setting up your python environment

In order to access this material, several python packages must be installed. A full
list can be seen inside `environment.yml`.

To set up the environment (using conda), the following can be done:

```bash
conda env create --name envname --file=environments.yml
jupyter contrib nbextension install --user
jupyter nbextension enable splitcell/splitcell
jupyter nbextension enable rubberband/main
jupyter nbextension enable exercise2/main
jupyter nbextension enable autosavetime/main
jupyter nbextension enable collapsible_headings/main
jupyter nbextension enable codefolding/main
jupyter nbextension enable limit_output/main
jupyter nbextension enable toc2/main
```

**Note: this workshop uses the beta release of MDAnalysis 2.0.0.**


## Binder

The tutorial materials can be accessed online via binder.
To launch the binder instance, click [here](https://mybinder.org/v2/gh/MDAnalysis/WorkshopPrace2021/HEAD).


## License

The code examples in this repository are licensed under the GPL version 2.0 license. Non-code content is licensed under the Creative Commons Attribution-ShareAlike 3.0 International license.
