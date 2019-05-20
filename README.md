# Curating metal-organic frameworks to compose robust gas sensor arrays in dilute conditions

Jupyter Notebook to produce data for:

> A. Sturluson, R. Sousa, Y. Zhang, M. T. Huynh, C. Laird, A. H. P. York, C. Silsby, C. Chang, C. M. Simon. Curating metal-organic frameworks to compose robust gas sensor arrays in dilute conditions.

The notebooks can be viewed without having the necessary programming languages (Julia 1.1 and Python 3), but in order to run the cells within the notebooks, you will need Julia 1.1/Python 3</br>
[Click me for Julia 1.1](https://julialang.org/)</br>
[Click me for Python 3](https://www.python.org/)</br>
# Structure of repository
```
.
├── ...
├── src/
│   ├── Sensing.ipynb                     # A Julia v1.1 notebook used to find metal-organic frameworks to compose robust sensor arrays
│   ├── scatter_plot_matrix.ipynb         # A Python 3 notebook used to create a scatter plot matrix of the Henry constants of the MOFs
│   └── Extract_Henry.ipynb               # A Julia v1.1 notebook used to extract Henry constants from experimental CO2 and SO2 adsorption isotherms
└── data/
    ├── expt_data/                        # A directory containing `.csv` files of all CO2 and SO2 adsorption isotherms used in this study, alongside the figures that were used to extract the data from (see references below)
    ├── cif_viz/                          # A directory containing generated images of the MOFs (using the VisIt visualization software)
    └── henry_constants.csv               # A file containing the extracted Henry constants from all adsorption isotherms (See `src/Extract_Henry.ipynb`)
```

# MOFs used in the study
The following MOFs were used in the paper. The choice of MOFs was determined by the availability of CO<sub>2</sub> and SO<sub>2</sub> adsorption isotherms.

| MOF | CO<sub>2</sub> Reference | SO<sub>2</sub> Reference |
| :--- |    :----:      | :----:         |
|MFM-600| 10.1021/jacs.8b08433 | 10.1021/jacs.8b08433 |
|MFM-601| 10.1021/jacs.8b08433 | 10.1021/jacs.8b08433 |
|Mg-MOF-74| 10.1021/cm401270b | 10.1039/c1ee01720a |
|Ni(bdc)(ted)<sub>0.5</sub>| 10.1021/cm401270b | 10.1016/j.micromeso.2009.11.026 |
|Zn(bdc)(ted)<sub>0.5</sub>| 10.1021/cm401270b | 10.1016/j.micromeso.2009.11.026 |
|MFM-300(In)| 10.1002/adma.201602338 | 10.1002/adma.201602338 |
|NOTT-202a | 10.1021/ja401061m | 10.1021/ja401061m |
|NOTT-300 | 10.1038/nchem.1457 | 10.1038/nchem.1457 |
|Co<sub>3</sub>[Co(CN)<sub>6</sub>]<sub>2</sub> | 10.1021/ic902397w | 10.1021/ic902397w |
|Zn<sub>3</sub>[Co(CN)<sub>6</sub>]<sub>2</sub> | 10.1021/ic902397w | 10.1021/ic902397w |
|KAUST-7 | 10.1038/s41467-019-09157-2 | 10.1038/s41467-019-09157-2 |
|KAUST-8 | 10.1038/s41467-019-09157-2 | 10.1038/s41467-019-09157-2 |

 
