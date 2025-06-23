####################################################################################
################## Sediment Delivery Ratio modelling with InVEST ###################
####################################################################################
#                                                                                  #
#                        PLEASE READ BEFORE RUNNING                                #
#                                                                                  #
#   This script allows you to quantify and map overland sediment generation and    #
#   delivery to the stream by using the InVEST Sediment Delivery Ratio model,      #
#   which rely on climate (in particular rain intensity), soil properties,         #
#   topography, and vegetation, as well as anthropogenic factors.                  #
#                                                                                  #
#   The first part of the script sets up the working environment and loads and     # 
#   prepares the data required for the process. In particular, it ensures that     #
#   all the data is correctly aligned with the reference raster.                   #
#   The second part of the script is used firstly to write the python script       #
#   used to launch the InVEST model, and then to read it via a system command.     #
#   The final raster is obtained by reading the output(s) of the InVEST model      #
#   and aggregating them. The final raster is clipped to the study area raster.    # 
#                                                                                  #
#   The results of this script are automatically saved in the "ES" subfolder of    #
#   the project folder. Metadata file is automatically written in output folder.   #
#   For the script to run correctly, you only need to modify sections 'Paths' and  #
#   'Parameters' below.                                                            #
#   You don't need to make any other changes.                                      #
#                                                                                  # 
#   Make sure you select the kernel associated with the invest3141 environment     #
#   when you run the script.                                                       # 
#                                                                                  #
####################################################################################
####################################################################################

# Paths ----------------------------------------------------------------------------
    #  !! Make sure that the data you enter meets this requirements : 
    #  ==>  Respect the specified data extension
    #  ==>  Ensure compliance with the points required for input data to the InVEST model:
    #       https://invest.readthedocs.io/en/latest/models.html#sediment-delivery-ratio
    #       http://releases.naturalcapitalproject.org/invest-userguide/latest/en/sdr.html

    # Study area raster (.asc or .tif)
    study_area_path <- "path/to/data/..." 

    # Datas for InVEST model
    # Path to a digital elevation model raster (in meters ; .asc or .tif)
    dem_path <- "path/to/data/..."
    # Path to a map of rainfall erosivity, reflecting the intensity and duration of rainfall (MJ*mm/(h*ha*year) ; .asc or .tif)
    erosivity_path <- "path/to/data/..."
    # Path ot a map of soil erodibility, the susceptibility of soil particles to detachment and transport by rainfall and runoff (t*h*ha/(ha*MJ*mm) ; .tif ou .asc)
    erodibility_path <- "path/to/data/..."
    # Path to a landcover raster (.asc or .tif)
    lulc_path <- "path/to/data/..."
    # Path to a watersheds vector (.shp ; polygon/multipolygon)
    watersheds_path <- "path/to/data/..."
    # Path to a biophysical table with biophysical properties related to nutrient load and retention for each LULC class (.csv ; sep="," ; UTF-8 encoded)
    biophysical_table_path <- "path/to/data/..."

# Parameters ------------------------------------------------------------------------
    # Project
    # Name of the main shared project folder
    shared_directory <- "path/to/the/root/of/the/shared/folder"
    # Specify the name of an existing project or choose your new project's name
    # Please note that if you enter an existing project name, previously calculated results for this indicator may be overwritten.
    project_name <- "version name"  
    # Name of the pillar
    pillar_name <- "ES" 
    # Name of the indicator
    indicator_name <- "SDR"
    # Give a short descrition of the indicator
    description <- paste0("Capacity of natural elements to limit erosion and therefore retain sediment (InVEST model).", "\n", 
                          "http://releases.naturalcapitalproject.org/invest-userguide/latest/en/sdr.html#")

    # Datas and computing parameters 
    # CRS of your projet, e.g. to which your data are
    mycrs <- "epsg:2056" 
    # Pattern used to select which of the outputs of the InVEST model should be aggregated to obtain the final result.
    # See here: http://releases.naturalcapitalproject.org/invest-userguide/latest/en/sdr.html#interpreting-results
    pattern_invest_output <- "avoided_export"
    # Function to be applied to InVEST outputs to aggregate them. Not used if only one output of the InVEST model is retained. 
    aggregation_function <- "sum"

    # InVEST model
    # Threshold flow accumulation (in number of pixel) = number of upslope pixels that must flow into a pixel before it is classified as a stream
    threshold_flow_accumulation <- 140
    # Borselli K parameter (unitless) = calibration parameter that determines the shape of the relationship between hydrologic connectivity
    k_param <- 2
    # Maximum SDR value (ratio) = the maximum SDR value that a pixel can have. It is a function of soil texture. More precisely, it is defined as the fraction of topsoil particles finer than coarse sand.
    sdr_max <- 0.8
    # Borselli IC0 Parameter (unitless)
    ic_0_param <- 0.5
    # Maximul L value (unitless) = the maximum allowed value of the slope length parameter (L) in the LS factor. Values of that exceed this are thresholded to this value.
    l_max <- 122

    # Clean up options
    # Do you want to delete the contents of the "scratch" folder at the end of the calculation? 'YES' or 'NO'
    # Temporary results are saved here. 
    scratch_to_trash <- "YES"
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
            # Remove the spatial object from memory
            rm(spatial_object)
            return(output_path)
            }
    # -----------------------------------------------------------------------------------
    # 1.3) Working directories ----------------------------------------------------------
        # Path to the working directory
        work_directory <- file.path(shared_directory, "OUTPUTS","INDICATORS", project_name) 
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
        info <- c(paste0("Version : ", project_name, "\n",
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
        # save information from biophysical table
        info_biophysical_table <- capture.output(print(biophysical_table))
    
        # Read input raster for InVEST model and make sure to align them with study area raster ; update the path to the corrected raster
        dem_aligned_path <- align_spatial_object(dem_path, study_area_raster, scratch_folder)
        erosivity_aligned_path <- align_spatial_object(erosivity_path, study_area_raster, scratch_folder)
        erodibility_aligned_path <- align_spatial_object(erodibility_path, study_area_raster, scratch_folder)
        lulc_aligned_path <- align_spatial_object(lulc_path, study_area_raster, scratch_folder)
        watersheds_aligned_path <- align_spatial_object(watersheds_path, study_area_raster, scratch_folder)
    
        end_time <- Sys.time()
        time_loading <- difftime(end_time, start_time, units="mins")
        # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(), "done"), tracking_file, append=TRUE)
        # Complete metadata file with running time 
        info <- c(info, paste0("INPUTS : ", "\n",
                               " * Study area: ", study_area_path, "\n",
                               " * DEM raster: ", dem_path, "\n",
                               " * Precipitation erosivity raster: ", erosivity_path, "\n", 
                               " * Soil erodibility raster: ", erodibility_path, "\n", 
                               " * LULC raster: ", lulc_path, "\n",
                               " * Vector watersheds : ", watersheds_path, "\n",
                               " * Biophysical table : ", biophysical_table_path))
        info <- c(info, info_biophysical_table)
        info <- c(info, paste0(" * Reading time: ", time_loading, " min", "\n\n"))
    # 2) Modelling ----------------------------------------------------------------------
    # 2.1) Prepare the Python script ----------------------------------------------------
        write(paste(Sys.time(), "WRITE PYTHON SCRIPT"), tracking_file, append=TRUE)
    
        # Set arguments for the InVEST model
        invest_arguments <- paste0("args ={", "\n",
                                   "'workspace_dir':", "'",scratch_folder,"'", "," , "\n",
                                   "'results_suffix':", "'",project_name, "'","," , "\n",
                                   "'dem_path':", "'",dem_aligned_path, "'","," , "\n",
                                   "'erosivity_path':", "'",erosivity_aligned_path, "'","," , "\n",
                                   "'erodibility_path':", "'",erodibility_aligned_path, "'","," , "\n",
                                   "'lulc_path':", "'",lulc_aligned_path, "'","," , "\n",
                                   "'watersheds_path':","'",watersheds_aligned_path, "'", "," , "\n",
                                   "'biophysical_table_path':","'",biophysical_table_path, "'", "," , "\n",
                                   "'threshold_flow_accumulation':",threshold_flow_accumulation,"," , "\n",
                                   "'k_param':",k_param, "," , "\n",
                                   "'sdr_max':",sdr_max, "," , "\n",
                                   "'ic_0_param':",ic_0_param, "," , "\n",
                                   "'l_max':",l_max, "," , "\n","}")
        # Initiate the .py file
        py_script_path <- paste(scratch_folder, "investSDR.py", sep="/")
        file(py_script_path, open="a")
        # Prepare all the lines
        py_script <- paste0("import natcap.invest.sdr.sdr", "\n", # Import the sediment delivery ratio InVEST model from the natcap module
                            "import csv", "\n", # Import csv module to read table
                            invest_arguments, "\n", # Read the arguments
                            "natcap.invest.sdr.sdr.execute(args)") # Execute the model 
        # Save lines in the .py file
        write(py_script, py_script_path, append=TRUE)
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
        info <- c(info, paste0("PYTHON SCRIPT FOR INVEST MODEL", "\n",
                               "##", "\n",
                               py_script, "\n",
                               "##", "\n",
                               "Computation time: ", time_prediction, " min", "\n\n"))
    # -----------------------------------------------------------------------------------
    # 2.3) Final prediction ------------------------------------------------------------
        write(paste(Sys.time(), "FINAL PREDICTION"), tracking_file, append=TRUE)
        start_time <- Sys.time()
    
        # Read InVEST output raster(s) and ensure projection
        predictions_list <- list.files(scratch_folder, pattern=pattern_invest_output, full.names=TRUE)
        predictions_stack <- rast(predictions_list)
        crs(predictions_stack) <- mycrs
    
        # Obtain the final raster by aggregating the various InVEST outputs (if required)
        final_raster <- app(predictions_stack, fun = match.fun(aggregation_function))
    
        # Mask final raster with the study area
        final_raster <- mask(final_raster, study_area_raster)
    
        # Save the result  
        writeRaster(final_raster, file.path(output_folder, "sediments_retention_pred.tif"), overwrite=TRUE)

        end_time <- Sys.time()
        time_finalising <- difftime(end_time, start_time, units="mins")
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(), "done"), tracking_file, append=TRUE)
        # Complete metadata file with outputs' informations 
        info <- c(info, paste0("OUTPUTS", "\n",
                               " * Final prediction saved in the: ", output_folder, "\n",
                               " * Resolution : ", res(final_raster)[1], " m", "\n",
                               " * CRS : ", crs(final_raster, describe=TRUE)[1], "\n",
                               " * Number of rasters compiled to obtain the final prediction: ", length(predictions_list)))
        info <- c(info, paste0("  ", predictions_list))
        info <- c(info, paste0(" * Aggregation function: ", aggregation_function, "\n",
                               " * Completion time: ", time_finalising, " min"))
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
write(paste0("##### Total duration: ", total_run_time, " min #####"), file.path(output_folder, "METADATA.txt"), append =TRUE) 

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
