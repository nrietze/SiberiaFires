# Code and documentation
This repository contains the code and output for Rietze et al. (in prep.): Topography and pre-fire conditions shape the underestimated burn heterogeneity of Siberian tundra fires.

> [!WARNING]
> Update when officially done.
 > The necessary data is publicly available under [![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.12650945-blue)](https://doi.org/10.5281/zenodo.12650945).

## Repository structure:
> [!WARNING]
> Update when officially done.
Here is the structure of this repo, files have been excluded from this tree for readability.

```bash
├───code
│   │
│   │  ZOIB_model.R
│   │
│   ├───classification
│   │   │  prep.R
│   │   │  predict_burned_area.R
│   │   │  predict_water.py
│   │
│   └───figures_and_tables
│          Table_1.R
│          figure_1.R
│          figure_2.R
│          figure_3.R
│    
├───data
│   └───geodata
│       ├───feature_layers
│       └───raster
│
├───figures  
│
└───tables
```

- The scripts in `classification` are used to prepare the PlanetScope imagery and run the random forest classification.
  - `prep.R` is used to crop, rename and prepare raster files for the image classification.
  - `predict_burned_area.R` is used to execute the classification of burned areas, performing validation and predicting the burned area maps.
- The scripts in `figures_and_tables` are used to generate the main and supplementary figures as well as the supporting tables and are named appropriately.
- The script `ZOIB_model.R` contains the code for the zero-one inflated beta regression. 
- The folder `data` is empty and should contain the data that can be downloaded from Zenodo (see link on top).
- The folder `figures` contains the figures that are produced in the correspoinding scripts.
- The folder `tables` contains the tables that are produced in the correspoinding scripts.

[to top](https://github.com/nrietze/SiberiaFires/main/README.md)

## Software requirements
The data pre-processing and data analysis was using R 4.2.2 (2022-10-31 ucrt). Newer versions of these software packages will likely work, but have not been tested.

Code development and processing were carried out in Windows 10 (64 bit), but execution should (in theory) be platform independent.

[to top](https://github.com/nrietze/SiberiaFires/main/README.md)

## Contact
Code development and maintenance: Nils Rietze ([nils.rietze@uzh.ch](nils.rietze@uzh.ch))

[to top](https://github.com/nrietze/SiberiaFires/main/README.md)

## Acknowledgements
> [!WARNING]
> Update when officially done.

From the manuscript:
*N.R. was supported through the TRISHNA Science and Electronics Contribution (T-SEC), an ESA PRODEX project (contract no. 4000133711).*

[to top](https://github.com/nrietze/SiberiaFires/main/README.md)
<!--- ## License
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
-->

## Citation
When citing elements in this repository, please cite as:

Rietze, N., Assmann, J. J., Heim, R., Troeva, E., Schaepman-Strub, G. (in prep.). Topography and pre-fire conditions shape the underestimated burn heterogeneity of Siberian tundra fires. 

[to top](https://github.com/nrietze/SiberiaFires/main/README.md)

## License

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

[to top](https://github.com/nrietze/SiberiaFires/main/README.md)
