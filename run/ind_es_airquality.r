#####################################################################################
######################## Air quality regulation based on LAI ########################
#####################################################################################
#                                                                                   #
#                        PLEASE READ BEFORE RUNNING                                 #
#                                                                                   #
#   This script allows you to estimate air quality by computing a raster of the     #
#   Leaf Area Index (LAI) based on an aggregate of several rasters of Normalized    #
#   Difference Vegetation Index (NDVI).                                             #
#                                                                                   #
#   The first part of the script sets up the working environment and loads and      # 
#   prepares the data required for the process. Next, the NDVI rasters are aggre-   #
#   gated according to the formula specified, to form a single NDVI raster. In the  #
#   last part of the script, this NDVI raster is converted into an LAI raster       #
#   according to the formula entered by the user. The final raster is then cut from #
#   the study area raster.                                                          # 
#                                                                                   #
#   The results of this script are automatically saved in the "SE" subfolder of the #
#   project folder.Metadata file is automatically written in output folder.         #
#   For the script to run correctly, you only need to modify sections 'Paths' and   #
#   'Parameters' below.                                                             #
#   You don't need to make any other changes.                                       #
#                                                                                   # 
#   Make sure you select the kernel associated with the rspatial environment when   #
#   you run the script.                                                             # 
#                                                                                   #
#####################################################################################
#####################################################################################

# Paths ----------------------------------------------------------------------------
#  !! Make sure that the data you enter meets the following requirements : 
    #  ==>   Spatial data must share the same CRS (and for rasters the same resolution)
    #  ==>   Respect the specified data extension

    # Study area raster (.asc or .tif). Will be used as a reference for the CRS, the resolution and the extent of the outputs of this script
    study_area_path <- "path/to/data/..." 
    # Folder path containing NDVI rasters (which are in .asc or .tif)
    ndvi_folder_path <- "path/to/data/..."
    
# Parameters ------------------------------------------------------------------------
    # Project
    # name of the main shared project folder
    shared_directory <- "path/to/the/root/of/the/shared/folder"
    # specify the name of an existing project or choose your new project's name
    # Please note that if you enter an existing project name, previously calculated results for this indicator may be overwritten.  
    project_name <- "version name" 
    # Name of the pillar
    pillar_name <- "ES" 
    # Name of the indicator
    indicator_name <- "AIRQUALITY"
    # Give a short descrition of the indicator
    description <- "Filtration of particles by the leaf surface, estimated using the Leaf Area Index (LAI)."

    # Datas and computing parameters 
    # CRS of your projet, e.g. to which your data are 
    mycrs <- "epsg:2056"        
    # Function to be used when aggregating NDVI rasters located in the 'ndvi_folder_path' folder 
    # (either 'max', 'min', 'median', 'sum', 'mean')
    aggregation_function <- "max"
    # Enter here (in string form) the conversion formula to be applied to the NDVI to obtain the LAI.
    # !! Make sur to use the letter 'r' to designate NDVI raster
    lai_conversion_formula <- "0.57*exp(2.33*r)"

    # Clean up options
    # Do you want to delete the contents of the "scratch" folder at the end of the calculation? 'YES' or 'NO'
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
        # raster_calculator(f, r) = apply a formula 'f' (string) to the values of a raster 'r'. The formula 'f' must have 'r' as its unique unknown. 
        raster_calculator <- function(f, r) {
            eval(parse(text = f))
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
        write(paste(Sys.time(), "LOAD AND PREPARE DATA"), tracking_file, append=TRUE)
        start_time <- Sys.time()
    
        # Coordinates reference systeme 
        mycrs <- crs(mycrs)
    
        # Study area
        study_area_raster <- rast(study_area_path)
        crs(study_area_raster) <- mycrs
    
        # Get a list of NDVI rasters' path
        ndvi_list <- list.files(ndvi_folder_path, full.names=TRUE)

        end_time <- Sys.time()
        time_loading <- difftime(end_time, start_time, units="mins")
    # -----------------------------------------------------------------------------------
        # fill the progress tracking file
        write(paste(Sys.time(), "done"), tracking_file, append=TRUE)
        # complete metadata file with inputs' informations 
        info <- c(info, paste0("INPUTS : ", "\n",
                               " * Study area: ", study_area_path, "\n",
                               " * NDVI rasters: ", ndvi_folder_path, "\n"))
        info <- c(info, paste0(ndvi_list, "\n"))
        info <- c(info, paste0(" * Reading time: ", time_loading, " min", "\n"))
    # -----------------------------------------------------------------------------------
    # 2) Computing Leaf Area Index ------------------------------------------------------
    # 2.1) Process NDVI rasters ---------------------------------------------------------
        write(paste(Sys.time(), "PROCESSING NDVI RASTERS"), tracking_file, append=TRUE)
        start_time <- Sys.time()
    
        # Stack all NDVI rasters 
        ndvi_stack <- rast(ndvi_list)
        
        # Initialise the synthaxis to aggregate the NDVI rasters according to the 'aggregation_function' function determined by the user. 
        # We want the NA present in certain rasters not to be taken into account in the aggregation.
        aggreg_function <- paste0(aggregation_function, "(r, na.rm=TRUE)")
    
        # Use 'raster_calculator()' to apply the 'aggreg_function' to the NDVI rasters
        final_ndvi <- raster_calculator(f=aggreg_function, r=ndvi_stack)
        
        # Initialize a matrix to reclass the remaining NA
        m <- c(NA, 0)
        mx <- matrix(m, byrow=TRUE, ncol=2)
        # Make sure there is no NA left in the final NDVI raster
        final_ndvi <- classify(final_ndvi, mx)

        # Save result
        writeRaster(final_ndvi, file.path(scratch_folder, paste0(aggregation_function, "_NDVI.tif")), overwrite=TRUE)
    
        end_time <- Sys.time()
        time_process_ndvi <- difftime(end_time, start_time, units="mins")
        # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(), "done"), tracking_file, append=TRUE)
        # Complete metadata file with computing and output information
        info <- c(info, paste0("NDVI RASTERS PROCESSING: ", "\n",
                               " * Number of rasters compiled: ", length(ndvi_list), "\n",
                               " * Function applied to raster aggregation: ", aggregation_function, "\n",
                               " * Computation time: ", time_process_ndvi, " min", "\n\n"))
    # 2.2) Convert NDVI in LAI raster ---------------------------------------------------
        write(paste(Sys.time(), "COMPUTING LAI"), tracking_file, append=TRUE)
        start_time <- Sys.time()
    
        # Apply 'lai_conversion_formula' to the NDVI raster with 'raster_calculator()' function
        lai_raster <- raster_calculator(f=lai_conversion_formula, r=final_ndvi)
        
        # Mask the final raster with the study area
        lai_raster <- mask(lai_raster, study_area_raster)

        # Save final raster
        writeRaster(lai_raster,file.path(output_folder,"LAI.tif"), overwrite=TRUE)
    
        end_time <- Sys.time()
        time_computing_lai <- difftime(end_time, start_time, units="mins")
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(), "done"), tracking_file, append=TRUE)
        # Complete metadata file with computing and output information
        info <- c(info, paste0("LAI COMPUTATION: ", "\n",
                               " * LAI = ", lai_conversion_formula, ", where 'r' represent NDVI raster", "\n",
                               " * Computation time: ", time_computing_lai, " min", "\n\n"))
        info <- c(info, paste0("OUTPUTS", "\n",
                               " * Results saved in the: ", output_folder, "\n",
                               " * Resolution : ", res(lai_raster)[1], " m", "\n",
                               " * CRS : ", crs(lai_raster, describe=TRUE)[1], "\n"))
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
write(paste0("##### Total duration : ", total_run_time, " min #####"), file.path(output_folder, "METADATA.txt"), append =TRUE) 

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
