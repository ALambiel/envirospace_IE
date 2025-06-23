# Identifying Ecological Infrastructure
This repository contains R scripts, configuration files and documentation for the identification and update of Ecological Infrastructure (EI), according to the method proposed by [GE21](https://www.ge21.ch/). 
This work has been designed to be used via a shared environment on JupyterHub and sripts have therefore been written to being run in this kind of environment. 

## 📌 Overview 
The aim of these scripts is to produce raster layers for each of the indicators (see figure below) used to build the pillars of the EI, and then perform a hierarchical prioritization of all this information in order to identify EI. The users decide on the input data to be used as well as on certain methodological parameters. Then, scripts run appropriate models/softwares to process data and produce documented outputs.
The methodology is implemented through an automated **R-based processing pipeline** which is described in this document.

Our objectives:
- Automate the computation of ecological indicators from heterogeneous datasets
- Centralize data and scripts
- Facilitate regular updates of the EI with minimal user intervention
- Ensure full traceability and reproducibility of results.
This project is part of a broader effort to enhance reproducibility, transparency, and scalability in environmental monitoring and conservation planning.

## 📖 Methodology and references 
The methodology presented here is based on an assessment of the ecological quality of the territory, hereafter named biodiversity diagnosis. The biodiversity diagnosis relies on four pillars, each representing a set of indicators that collectively aim to characterize biodiversity from a specific perspective. To facilitate the integration of these various dimensions of biodiversity, the spatial prioritization tool [Zonation](https://zonationteam.github.io/Zonation5/) is employed to rank each spatial unit. Using an iterative algorithm, each cell across the entire territory is assigned a unique value reflecting its relative importance compared to others. From the biodiversity diagnosis, 30% of the territory with the highest priority score can be determined, forming the Geneva Ecological Infrastructure according to target 3 of the CBD. This represents an area where biodiversity covers the most representative proportion of the selected dimensions to consider in the systematic conservation planning

<img src="https://github.com/ALambiel/envirospace_IE/blob/main/images/methodology_summary.png" data-canonical-src="https://github.com/ALambiel/envirospace_IE/blob/main/images/methodology_summary.png" width="1000" />

For more informations, please refers to: 
- Honeck, Erica, Arthur Sanguet, Martin A. Schlaepfer, Nicolas Wyler, and Anthony Lehmann. 2020. “Methods for Identifying Green Infrastructure.” *SN Applied Sciences* 2 (11): 1916. [https://doi.org/10.1007/s42452-020-03575-4](https://doi.org/10.1007/s42452-020-03575-4).
- Sanguet, Arthur, Nicolas Wyler, Benjamin Guinaudeau, Noé Waller, Loreto Urbina, Laurent Huber, Claude Fischer, and Anthony Lehmann. 2023. “Mapping Ecological Infrastructure in a Cross-Border Regional Context.” *Land* 12 (11): 2010. [https://doi.org/10.3390/land12112010](https://doi.org/10.3390/land12112010). 

## 📂 Folder structure

```
envirospace_IE/
│── config/                            # contains files for deploying required environments
│   ├── rspatial.yml                   # 1st (main) environment with R packages for spatial analysis 
│   ├── invest3141.yml                 # 2nd environment with R and Python packges/modules required to run InVEST 3.14.1
│── indicators/                        # contains scripts for indicator calculation                    
│   ├── pillar_indicator.r             # one script per indicator
│   ├── ...
│── prioritization/                    # contains scripts for prioritization 
│── README.md                          # this documentation
│── tools/                             # additional scripts for optional pre-processing
```

## 🚀 Getting started
From your `/home` workspace on JupyterLab, open a **terminal** and:

1. Clone this repository
```bash
git clone https://github.com/ALambiel/envirospace_IE
cd envirospace_IE
```
Steps 2 and 3 have to be done both for `rspatial` and `invest3141`.

2. Deploy `rspatial` environment and create associated kernel
```bash
conda env create -f config/rspatial
conda activate rspatial

R
IRkernel::installspec(name = 'rspatial', displayname = 'Spatial analysis with R')
q()
```

3. Deploy `invest3141` environment and create associated kernels
```bash
conda env create -f config/invest3141
conda activate invest3141
pip install natcap.invest==3.14.1

R
IRkernel::installspec(name = 'invest3141r', displayname = 'InVEST 3.14.1 with R')
q()

python -m ipykernel install --user --name invest3141 --display-name "InVEST 3.14.1"
```

4. Run R script

## 🤝 Reproducibility & FAIR Principles
This project adheres to the FAIR principles (Findable, Accessible, Interoperable, Reusable) and promotes open, transparent, and reproducible science. All scripts are annotated, and metadata is generated automatically to ensure traceability.

## 🛡 Licence
To be defined

## Aknowledgment
The authors would like to thank the Swiss Federal Office of the Environment (FOEN) for their financial support and express their gratitude to the RPT4 pilot committee. The authors gratefully acknowledge the members of the GE21 expert group from the University of Geneva (UNIGE), the Geneva School of Engineering, Architecture and Landscape (HEPIA), the Conservatory and Botanical Garden of Geneva (CJB) and the Cantonal Office for Agriculture and Nature (OCAN) for their previous work on establishing indicators and, more broadly, on developing the method for identifying ecological infrastructure, from which the present work is a continuation.

## Related publication
Lambiel, Audrey, Gregory Giuliani, Nathan Külling, and Anthony Lehmann (in preparation). A Digital Twin approach for the identification and update of Ecological Infrastructure.
