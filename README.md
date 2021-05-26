# Using MDAnalysis for Efficient Simulation Pre- and Post-Processing [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/MDAnalysis/WorkshopPrace2021/HEAD)

## Contents

This repository contains the materials for the MDAnalysis 2021 PRACE/SURF workshop.

The directories contain lectures and practical materials for seassion over the 3 days.


## Setting up python environment

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


## Binder

The tutorial materials can be accessed online via binder.
To launch the binder instance, click [here](https://mybinder.org/v2/gh/MDAnalysis/WorkshopPrace2021/HEAD).


## License

The code examples in this repository are licensed under the GPL version 2.0 license. Non-code content is licensed under the Creative Commons Attribution-ShareAlike 3.0 International license.
