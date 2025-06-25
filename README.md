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

This repository follows the above structure:
```
envirospace_IE/
│── config/                            # contains files for deploying required environments
│   ├── rspatial.yml                   # 1st (main) environment with R packages for spatial analysis 
│   └── invest3141.yml                 # 2nd environment with R and Python packges/modules required to run InVEST 3.14.1
│── run/
│   ├── <ind_pillar_indicator.r>       # one script for each indicator
│   ├── ... 
│   |   └── ...
│   └── tools/                         # additional scripts for optional pre-processing
└── README.md                          # this documentation
```
Execution of the scripts (indicators or prioritization) creates new `OUTPUTS` folder, structured as follows:

```
working_directory/                     # depending on the folder you enter in the script
└── OUTPUTS/
    ├── INDICATORS/
    |   └── <version>/                 # sub-folder to store all results from a single version
    |      ├── <pillar>/               # one for each pillar, either “structure”, “ecosystem service”, ‘function’ or “composition”
    |      |   └── <indicator>/        # one for each indicator that forms a pillar
    |      |      ├── <result.tif>     # one or more rasters
    |      |      └── METADATA.txt     # an automatically generated file with information on input data/parameters/etc.  
    |      └── scratch/                # contains any intermediate results if retained 
    └── PRIORITIZATION/
        └── <version>/                # One subfolder per prioritization version
            ├── zonation_output/
            │   ├── rankmap.tif       # Final prioritization output raster
            │   └── METADATA.txt      # Summary of processing and parameters
            └── zonation_settings/
                ├── settings.z5       # Zonation settings file
                ├── all_files.txt     # List of input rasters with weights/groups
                └── weight_group.txt  # Optional: group weights file (if used)
```

## 🚀 Getting started : environments and indicators
From your `/home` workspace on JupyterLab, open a **terminal** and:

1. Clone this repository
```bash
git clone https://github.com/ALambiel/envirospace_IE
cd envirospace_IE
```

2. Deploy `rspatial` environment and create associated kernel
```bash
conda env create -f config/rspatial.yml
conda activate rspatial

R
IRkernel::installspec(name = 'rspatial', displayname = 'Spatial analysis with R')
q()
```

3. Deploy `invest3141` environment and create associated kernels
```bash
conda env create -f config/invest3141.yml
conda activate invest3141
pip install natcap.invest==3.14.1

R
IRkernel::installspec(name = 'invest3141r', displayname = 'InVEST 3.14.1 with R')
q()

python -m ipykernel install --user --name invest3141 --display-name "InVEST 3.14.1"
```

4. Run script

Open any script from `run` folder, read it and modify: 
- data/folder path
- parameters

Then, execute the code by launching a console and selecting the correct kernel.

## 🚀 Getting started : Prioritization with Zonation 5
To run the **prioritization script** (`prio_zonation.R`), you’ll need to install Zonation 5 using its AppImage.

1. Download and extract Zonation AppImage 
 [Zonation5_Linux.zip (latest release)](https://github.com/zonationteam/Zonation5/releases/latest/download/Zonation5_Linux.zip)

2. Unzip the archive and load  `Zonation5` file in your JupyterLab environment
3. Make it executable and extract it, by running this in your terminal, with the correct path
```bash
chmod +x correctpath/Zonation5
correctpath/Zonation5 --appimage-extract
```

This creates a folder called `squashfs-root/`, which contains the executable.

4. Update Your R Script

In the `prio_zonation.R` script, make sure to set the executable path:
```r
z5_exe <- "/your/path/to/squashfs-root/usr/bin/z5"
```
and modify all the input parameters to fit your prioritization purposes

5. Launch the Prioritization

Now run the script from JupyterLab using the `rspatial` environment. Zonation will read the indicator layers, generate a prioritization based on weights/groups/masks, and export the results.

## 📈 Future Improvements
- Production of a tutorial/practical guide 
- Provision of a test dataset

## 🤝 Reproducibility & FAIR Principles
This project adheres to the FAIR principles (Findable, Accessible, Interoperable, Reusable) and promotes open, transparent, and reproducible science. All scripts are annotated, and metadata is generated automatically to ensure traceability.

## 🛡 Licence
©UNIGE/GE21 2025. This work is licensed under Creative commons Attribution-NonCommercial-ShareAlike 4.0 International terms.

<img src="https://mirrors.creativecommons.org/presskit/buttons/88x31/svg/by-nc-sa.svg" data-canonical-src="https://mirrors.creativecommons.org/presskit/buttons/88x31/svg/by-nc-sa.svg"/> 

You are free to use, modify, and distribute this software under the terms of the CC-BY-NC-SA 4.0. For more details, see the full license text: [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/).

## Aknowledgment
The authors would like to thank the Swiss Federal Office of the Environment (FOEN) for their financial support and express their gratitude to the RPT4 pilot committee. The authors gratefully acknowledge the members of the GE21 expert group from the University of Geneva (UNIGE), the Geneva School of Engineering, Architecture and Landscape (HEPIA), the Conservatory and Botanical Garden of Geneva (CJB) and the Cantonal Office for Agriculture and Nature (OCAN) for their previous work on establishing indicators and, more broadly, on developing the method for identifying ecological infrastructure, from which the present work is a continuation.

## Related publication
Lambiel, Audrey, Gregory Giuliani, Nathan Külling, and Anthony Lehmann (in preparation). A Digital Twin approach for the identification and update of Ecological Infrastructure.
