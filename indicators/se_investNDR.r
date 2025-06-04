####################################################################################
################## Nutrient Delivery Ratio modelling with InVEST ###################
####################################################################################
#                                                                                  #
#                        PLEASE READ BEFORE RUNNING                                #
#                                                                                  #
#   This script allows you to map nutrient sources from watersheds and their       # 
#   transport to the stream, as a way to to assess the service of nutrient         #
#   retention by natural vegetation. It using the InVEST Nutrient Delivery Ratio   #                         
#   model, which uses data on the landscape (topography, land use types,           #
#   precipitation) to estimate how water and nutrients move through the soil to    #
#   water bodies.                                                                  # 
#                                                                                  #
#   The first part of the script sets up the working environment and loads and     # 
#   prepares the data required for the process. In particular, it ensures that     #
#   all the data is correctly aligned with the reference raster.                   #
#   The second part of the script is used firstly to write the python script       #
#   used to launch the InVEST model, and then to read it via a system command.     #
#   The final raster is obtained by reading the output(s) of the InVEST model      #
#   and aggregating them. The final raster is clipped to the study area raster.    # 
#                                                                                  #
#   The results of this script are automatically saved in the "SE" subfolder of    #
#   the project folder. Metadata file is automatically written in output folder.   #
#   For the script to run correctly, you only need to modify sections 'Paths' and  #
#   'Parameters' below.                                                            #
#   You don't need to make any other changes.                                      #
#                                                                                  # 
#   Make sure you select the kernel associated with the invest3141.yml environment #
#   when you run the script.                                                       # 
#                                                                                  #
####################################################################################
####################################################################################

# Paths ----------------------------------------------------------------------------
    #  !! Make sure that the data you enter meets this requirements : 
    #  ==>  Respect the specified data extension
    #  ==>  Ensure compliance with the points required for input data to the InVEST model:
    #       https://invest.readthedocs.io/en/latest/models.html#nutrient-delivery-ratio
    #       http://releases.naturalcapitalproject.org/invest-userguide/latest/en/ndr.html

    # Study area raster (.asc or .tif)
    study_area_path <- "envirospace/projects/GE21/IE/..." 

    # Datas for InVEST model
    # Path to a digital elevation model raster (in meters ; .asc or .tif)
    dem_path <- "envirospace/projects/GE21/IE/..."
    # Path to a landcover raster (.asc or .tif)
    lulc_path <- "envirospace/projects/GE21/IE/..."
    # Path to a runoff proxy raster (unitless ; .asc or .tif)
    runoff_proxy_path <- "envirospace/projects/GE21/IE/..."
    # Path to a watersheds vector (.shp ; polygon/multipolygon)
    watersheds_path <- "envirospace/projects/GE21/IE/..."
    # Path to a biophysical table with biophysical properties related to nutrient load and retention for each LULC class (.csv ; sep="," ; UTF-8 encoded)
    biophysical_table_path <- "envirospace/projects/GE21/IE/..."
    

# Parameters ------------------------------------------------------------------------
    # Project
    # Name of the main shared project folder
    shared_directory <- "envirospace/projects/GE21/IE"
    # Specify the name of an existing project or choose your new project's name
    # Please note that if you enter an existing project name, previously calculated results for this indicator may be overwritten.
    project_name <- "test"  
    # Name of the pillar
    pillar_name <- "SE" 
    # Name of the indicator
    indicator_name <- "NDR"
    # Give a short descrition of the indicator
    description <- paste0("Rétention des nutriments (modèle InVEST)", "\n",
                          "http://releases.naturalcapitalproject.org/invest-userguide/latest/en/ndr.html")
    
    # Datas and computing parameters 
    # CRS of your projet, e.g. to which your data are
    mycrs <- "epsg:2056" 
    # Pattern used to select which of the outputs of the InVEST model should be aggregated to obtain the final result.
    # See here: http://releases.naturalcapitalproject.org/invest-userguide/latest/en/ndr.html#interpreting-results
    pattern_invest_output <- "effective_retention"
    # Function to be applied to InVEST outputs to aggregate them. Not used if only one output of the InVEST model is retained. 
    aggregation_function <- "sum"

    # InVEST model
    # Threshold flow accumulation (in number of pixel) = number of upslope pixels that must flow into a pixel before it is classified as a stream
    threshold_flow_accumulation <- 140
    # Borselli K parameter (unitless) = calibration parameter that determines the shape of the relationship between hydrologic connectivity
    k_param <- 2
    # Subsurface critical length for nitrogen (meters) = distance traveled (subsurface and downslope) after which it is assumed that soil retains nitrogen at its maximum capacity
    subsurface_critical_length_n <- 200
    # Subsurface maximum retention efficiency for nitrogen (ratio) = maximum nitrogen retention efficiency that can be reached through subsurface flow
    subsurface_eff_n <- 0.8

    # Clean up options
    # Do you want to delete the contents of the "scratch" folder at the end of the calculation? 'YES' or 'NO'
    scratch_to_trash <- "NO"
    # Do you want to delete the progress and error files generated while the script is running when it finishes? 'YES' or 'NO'
    # n.b.: If "YES" but an error occurs, the two files will not be deleted in all cases to allow debugging.  
    track_to_trash <- "YES"

####################################################################################
############################### DO NOT MODIFY ######################################

# Record the start time of the script
script_start_time <- Sys.time()

# Create a file for error management in the 'RUN' folder
# Get the path of the session
session_path <- Sys.getenv()[["JPY_SESSION_NAME"]]
# Separate all elements in the path
session_path <- strsplit(session_path, '/')[[1]]
# Isolate the name of the script
script_name <- tools::file_path_sans_ext(tail(session_path, n=1))
# Set the name of error file
error_file <- paste0(project_name, "_", script_name, "_err.out")
# Create it
err <- file(error_file, open = "a")

# Encapsulate the entire script in a tryCatch block
tryCatch({
    # -----------------------------------------------------------------------------------
    # 1) Preparation --------------------------------------------------------------------
    # 1.1) Environment and project ------------------------------------------------------
        # Get user name 
        user <- as.character(Sys.info()["user"])

        # Get the date
        date <- Sys.Date() 

        # Get the name of the current conda environment
        # Get the R_HOME environment variable, which gives the path to the R installation director
        r_home_path <- Sys.getenv()[["R_HOME"]]
        # Split it 
        r_home_path <- strsplit(r_home_path, '/')[[1]]
        # Find the index of the "envs" directory in the split path
        index_envs <- which(r_home_path == "envs")
        # Get the name of the environment, which is the directory following "envs" in the path
        env <- r_home_path[index_envs + 1]
    
        Sys.setenv('GTIFF_SRS_SOURCE'='EPSG')
        Sys.setenv('PROJ_LIB' = paste0("/home/",user,"/.conda/envs/", env, "/share/proj"))
    
        # Path to the python binary in the 'invest3141' environment where InVEST is installed
        env_path <- paste0("/home/", user, "/.conda/envs/", env, "/bin/python")
    # -----------------------------------------------------------------------------------
    # 1.2) Packages and functions -------------------------------------------------------
        # Packages
        # Custom 'install & load' function
        install.load.package <- function(x) {
            if (!require(x, character.only = TRUE)) {
                install.packages(x, repos = "http://cran.us.r-project.org")
            }
            require(x, character.only = TRUE)
        }
        # List of packages required  
        required_packages <- c("terra") # spatial data analysis
        # Executing install & load for each package
        sapply(required_packages, install.load.package)

        # Functions
        # align_spatial_object(input_path, reference_raster, output_raster) = aligns a spatraster or spatvector whose access path is given by 'input_path' with a 'reference_raster'. The spatial object thus aligned is then saved under 'output_path' with the same original name, and this updated path is returned. 
        align_spatial_object <- function(input_path, reference_raster, output_path) {
            # Load the input spatial object as a raster or vector based on file extension
            spatial_object <- if (grepl("\\.tif", input_path)) rast(input_path) else vect(input_path)
            # If spatial object is a SpatRaster
            if (inherits(spatial_object, "SpatRaster")) {
                # Check if all values in the raster are integers to select the resampling method
                is_integer <- all(as.integer(values(spatial_object)) == values(spatial_object), na.rm = TRUE)
                # Choose the resampling method based on whether the raster has integer values
                method <- if (is_integer) "near" else "bilinear"
                # Resample the raster to match the reference raster using the chosen method
                spatial_object <- resample(spatial_object, reference_raster, method = method)
                # Construct the output path using the input filename
                output_path <- file.path(output_path, basename(input_path))
                # Save the resampled raster to the output path
                writeRaster(spatial_object, output_path, overwrite = TRUE)
            # If spatial object is a spatvector    
            } else if (inherits(spatial_object, "SpatVector")) {
                # Project the vector to the same CRS as the reference raster
                spatial_object <- project(spatial_object, crs(reference_raster))
                # Crop the vector to the extent of the reference raster
                spatial_object <- crop(spatial_object, ext(reference_raster))
                # Construct the output path using the input filename
                output_path <- file.path(output_path, basename(input_path))
                # Save the processed vector to the output path
                writeVector(spatial_object, output_path, overwrite = TRUE)
                }
            return(output_path)
            }
    # -----------------------------------------------------------------------------------
    # 1.3) Working directories ----------------------------------------------------------
        # Path to the working directory
        work_directory <- file.path(shared_directory, "OUTPUTS","INDICATEURS", project_name) 
        # Folder for the many intermediate results, which can be deleted at the end
        scratch_folder <- file.path(work_directory, "scratch", paste(pillar_name, indicator_name, date, sep="_"))
        # Directory for final outputs
        output_folder <- file.path(work_directory, pillar_name, indicator_name)
    
        # Create directories if they don't exist
        dir.create(work_directory, showWarnings = FALSE, recursive = TRUE)
        dir.create(scratch_folder, showWarnings = FALSE, recursive = TRUE)
        dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
    
        # Initialising the tracking file 
        tracking_file <- paste0(project_name, "_", script_name, ".txt")
        writeLines(paste(Sys.time(), "RUNNING ..."), con=tracking_file)
        # Initialising the metadata text file
        info <- c(paste0("Projet : ", project_name, "\n",
                         "Date : ", date, "\n",
                         "User : ", user, "\n\n",
                         pillar_name, " - ", indicator_name, "\n",
                         description, "\n\n"))
    # -----------------------------------------------------------------------------------
    # 1.4) Load and prepare data --------------------------------------------------------
        write(paste(Sys.time(),"LOAD AND PREPARE DATA"), tracking_file, append=TRUE)
        start_time <- Sys.time()
    
        # Study area
        study_area_raster <- rast(study_area_path)
        crs(study_area_raster) <- mycrs

        # InVEST inputs
        # Get the list of the land use/land cover type 
        biophysical_table <- read.csv(biophysical_table_path, sep=",", header=TRUE)
        lulc_list <- unique(biophysical_table$LULC_name)
        # Save information from biophysical table
        info_biophysical_table <- capture.output(print(biophysical_table))
        
        # Read input rasters for InVEST model and make sure to allign them with study area raster ; update the path to the corrected raster
        dem_aligned_path <- align_spatial_object(dem_path, study_area_raster, scratch_folder)
        lulc_aligned_path <- align_spatial_object(lulc_path, study_area_raster, scratch_folder)
        runoff_proxy_aligned_path <- align_spatial_object(runoff_proxy_path, study_area_raster, scratch_folder)
        watersheds_aligned_path <- align_spatial_object(watersheds_path, study_area_raster, scratch_folder)
        
        end_time <- Sys.time()
        time_loading <- difftime(end_time, start_time, units="mins")
        # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(), "done"), tracking_file, append=TRUE)
        # Complete metadata file 
        info <- c(info, paste0("INPUTS : ", "\n",
                               " * Zone d'étude : ", study_area_path, "\n",
                               " * Raster DEM : ", dem_path, "\n",
                               " * Raster LULC  : ", lulc_path, "\n",
                               " * Raster 'runoff proxy' : ", runoff_proxy_path, "\n",
                               " * Vector watersheds : ", watersheds_path, "\n",
                               " * Table biophysique : ", biophysical_table_path, "\n"))
        info <- c(info, info_biophysical_table)
        info <- c(info, paste0(" * Temps de lecture : ", time_loading, " min", "\n\n"))
    # 2) Modelling ----------------------------------------------------------------------
    # 2.1) Prepare the Python script ----------------------------------------------------
        write(paste(Sys.time(), "WRITE PYTHON SCRIPT"), tracking_file, append=TRUE)
        
        # Set arguments for the InVEST model
        invest_arguments <- paste0("args ={", "\n",
                                   "'workspace_dir':", "'",scratch_folder,"'", "," , "\n",
                                   "'results_suffix':", "'",project_name, "'","," , "\n",
                                   "'dem_path':", "'",dem_aligned_path, "'","," , "\n",
                                   "'lulc_path':", "'",lulc_aligned_path, "'","," , "\n",
                                   "'runoff_proxy_path':","'",runoff_proxy_aligned_path, "'", "," , "\n",
                                   "'watersheds_path':","'",watersheds_aligned_path, "'", "," , "\n",
                                   "'biophysical_table_path':","'",biophysical_table_path, "'", "," , "\n",
                                   "'calc_p' : True,","\n",
                                   "'calc_n' : True,","\n",
                                   "'threshold_flow_accumulation':",threshold_flow_accumulation,"," , "\n",
                                   "'k_param':",k_param, "," , "\n",
                                   "'subsurface_critical_length_n':",subsurface_critical_length_n,"," , "\n",
                                   "'subsurface_eff_n':",subsurface_eff_n,"," , "\n",
                                   "}")
        # Initiate the .py file
        py_script_path <- paste(scratch_folder, "investNDR.py", sep="/")
        file(py_script_path, open="a")
        # Prepare all the lines
        py_script <- paste0("import natcap.invest.ndr.ndr", "\n", # import the nutrient delivery ratio invest model from the natcap module
                            "import csv", "\n", # import csv module to read table
                            invest_arguments, "\n", # read the arguments 
                            "natcap.invest.ndr.ndr.execute(args)") # execute the model 
        # Save lines in the .py file
        write(py_script, py_script_path)
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(), "done"), tracking_file, append=TRUE)
    # -----------------------------------------------------------------------------------
    # 2.2) Running InVEST model ---------------------------------------------------------
        write(paste(Sys.time(), "RUN INVEST MODEL"), tracking_file, append=TRUE)
        start_time <- Sys.time()
    
        # Call python script with command system
        cmd <- paste(env_path, py_script_path)
        system(cmd, intern=TRUE)
    
        end_time <- Sys.time()
        time_prediction <- difftime(end_time, start_time, units="mins")
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(), "done"), tracking_file, append=TRUE)
        # Complete metadata file with modelling parameters and running time 
        info <- c(info, paste0("SCRIPT PYTHON POUR LE MODELE INVEST", "\n",
                               "##", "\n",
                               py_script, "\n",
                               "##", "\n",
                               "Temps de calcul : ", time_prediction, " min", "\n\n"))
    # -----------------------------------------------------------------------------------
    # 2.3) Prepare final result ---------------------------------------------------------
        write(paste(Sys.time(), "PREPARE FINAL RESULT"), tracking_file, append=TRUE)
        start_time <- Sys.time()
    
        # Read InVEST output raster(s) and ensure projection
        predictions_list <- list.files(file.path(scratch_folder, "intermediate_outputs"), pattern=pattern_invest_output, full.names = TRUE)
        predictions_stack <- rast(predictions_list)
        crs(predictions_stack) <- mycrs
        
        # Obtain the final raster by aggregating the various InVEST outputs (if required)
        final_raster <- app(predictions_stack, fun = match.fun(aggregation_function))
        
        # Mask final raster with the study area 
        final_raster <- mask(final_raster, study_area_raster)
    
        # Save the result  
        writeRaster(final_raster, file.path(output_folder, "nutrients_effective_retention_n_p_pred.tif"), overwrite=TRUE)

        end_time <- Sys.time()
        time_finalising <- difftime(end_time, start_time, units="mins")
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(), "done"), tracking_file, append=TRUE)
        # Complete metadata file with outputs' informations 
        info <- c(info, paste0("OUTPUTS", "\n",
                               " * Prédiction finale sauvegardée dans le dossier : ", output_folder, "\n",
                               " * Résolution : ", res(final_raster)[1], " m", "\n",
                               " * CRS : ", crs(final_raster, describe=TRUE)[1], "\n",
                               " * Nombre de rasters compilés pour obtenir la prédiction finale : ", length(predictions_list)))
        info <- c(info, paste0("  ", predictions_list))
        info <- c(info, paste0(" * Fonction d'aggrégation : ", aggregation_function, "\n",
                               " * Temps de finalisation : ", time_finalising, " min"))
    # -----------------------------------------------------------------------------------
        # Save the metadata file in the output folder
        writeLines(info, file.path(output_folder, "METADATA.txt")) 
        # Update the progress tracking file
        write(paste(Sys.time(), "COMPLETED."), tracking_file, append=TRUE)
    # -----------------------------------------------------------------------------------
}, error = function(e) {
            # If an error occured, save the error message in the error file
            writeLines(paste(Sys.time(), ": An error occurred - ", conditionMessage(e)), err)
            write(paste(Sys.time(), "!! FAILED !!"), tracking_file, append=TRUE)
        }
    )

# Record the end time of the script
script_end_time <- Sys.time()
total_run_time <- difftime(script_end_time, script_start_time, units="mins")

# Update the metadata file with the total run time 
write(paste0("##### Durée totale : ", total_run_time, " min #####"), file.path(output_folder, "METADATA.txt"), append =TRUE) 

# Close error file
close(err)
# Manage scratch folder content
if(scratch_to_trash=="YES"){ # if the user has chosen "YES" to allow the contents of the sracth folder to be deleted 
    unlink(scratch_folder, recursive =TRUE) # deletion of the scratch folder
    }
# Manage files progress and error tracking 
if(track_to_trash=="YES" & file.size(error_file)==0){ # if the user has chosen "YES" to allow deletion AND that no errors have occurred 
    unlink(tracking_file) # deletion of the progress tracking file
    unlink(error_file) # deletion of the error tracking file 
    }