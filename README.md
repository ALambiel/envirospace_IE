# Identifying Ecological Infrastructure
This repository contains R scripts, configuration files, and documentation for the identification and update of Ecological Infrastructure (EI), according to the method proposed by GE21. 
This work has been designed to be used via a shared environment on JupyterHub. The scripts have therefore been written to being run in this kind of environment. 

## 📌 Overview 
The aim of these scripts is to produce raster layers for each of the indicators used to build the pillars of the ecological infrastructure. The users decide on the data to be used as well as on certain methodological parameters. Then, the scripts link this data to the appropriate tools/software (e.g., InVEST, MaxEnt, etc.) and produce aligned, standardized layers. Metadata associated with the outputs is generated automatically.
The project is part of a broader effort to enhance reproducibility, transparency, and scalability in environmental monitoring and conservation planning.
- Automate the computation of ecological indicators from heterogeneous datasets.
- Centralize data and scripts.
- Facilitate regular updates of the EI with minimal user intervention.
- Ensure full traceability and reproducibility of results.


## 📖 Methodology and references 

Figure (how to add .png?)

Honeck, Erica, Arthur Sanguet, Martin A. Schlaepfer, Nicolas Wyler, and Anthony Lehmann. 2020. “Methods for Identifying Green Infrastructure.” SN Applied Sciences 2 (11): 1916. https://doi.org/10.1007/s42452-020-03575-4.

Lambiel, Audrey, Gregory Giuliani, Nathan Külling, and Anthony Lehmann (in preparation) A Digital Twin approach for the identification and update of Ecological Infrastructure.
https://www.ge21.ch/

Sanguet, Arthur, Nicolas Wyler, Benjamin Guinaudeau, Noé Waller, Loreto Urbina, Laurent Huber, Claude Fischer, and Anthony Lehmann. 2023. “Mapping Ecological Infrastructure in a Cross-Border Regional Context.” Land 12 (11): 2010. https://doi.org/10.3390/land12112010. 

The methodology is implemented through an automated **R-based processing pipeline** which is described in this document.

## 📂 Folder structure

```
repo_root/
│── config/
│   ├── mainenv.yml
│   ├── invest3149.yml
│── indicators/                            
│   ├── pillar_indicator_additional.r
│   ├── ...
│── prioritization/ 
│── README.md                            # This documentation
```

## Getting started
From your personnal worspace on JupyterLab:

1. Clone this repository

```bash
git clone <repo_url>
cd <repo_name>
```

2. Deploy environments and install dependencies

```bash
conda env create -f <env_path>
conda activate <env_name>

# R kernel
R
IRkernel::installspec(name = 'kernel_name', displayname = 'Kernel name')
q()

# python kernel
python -m ipykernel install --user --name kernel_name --display-name "Kernel name"
```

3. Run script

## 🚀 Reproducibility & FAIR Principles
This project adheres to the FAIR principles (Findable, Accessible, Interoperable, Reusable) and promotes open, transparent, and reproducible science. All scripts are annotated, and metadata is generated automatically to ensure traceability.

## 🛡 Licence
To be defined
