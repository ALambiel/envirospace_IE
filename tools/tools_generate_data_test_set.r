####################################################################################
########################## Generate a data test set ################################
####################################################################################
#                                                                                  #
#                        PLEASE READ BEFORE RUNNING                                #
#                                                                                  #
#   This script allows you to crop raster and vector data based on a specified     #
#   study area. It uses the 'terra' packages to perform spatial operations.        #
#                                                                                  #
#                                                                                  # 
#   Make sure you select the kernel associated with the mainenv.yml environment    #
#   when you run the script.                                                       # 
#                                                                                  #
####################################################################################
####################################################################################


# Paths ----------------------------------------------------------------------------
    # list to store paths of input raster data
    input_raster_data_path <- list("/home/path/to/raster1.tif",
                                  "/home/path/to/raster2.tif")
    # List to store paths of input vector data
    input_vector_data_path <- list("/home/path/to/vector1.shp",
                                  "/home/path/to/vector2.shp")
    # Path to the study area file (either raster or vector)
    study_area_path <- "/home/path/to/studyarea.tif"
    # specify if the study area object is either a 'raster' or a 'vector' type object
    study_area_type <- "raster"
    # Specify either a 'reduction_factor' (as numeric) or the desired extent of the 'bounding_box' directly
    # The other parameter must be set to NULL
    reduction_factor <- 50 
    bounding_box <- NULL # (xmin, xmax, ymin, ymax)
    output_data_path <- "/home/path/to/output/folder"

####################################################################################
############################### DO NOT MODIFY ######################################

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
        library(terra) # spatial data analysis

        # Functions
        # Function to crop raster data
        crop_raster <- function(raster_path, extent_obj, output_path) {
            # extract the file name from path 
            file_name <- basename(raster_path)
            # load raster data 
            raster <- rast(raster_path)
            # crop raster data to extent 
            cropped_raster <- crop(raster, extent_obj)
            # save it in the outout folder
            writeRaster(cropped_raster, file.path(output_path, file_name), overwrite=TRUE)
            return(cropped_raster)
        }

        # Function to crop vector data
        crop_vector <- function(vector_path, extent_obj, output_path) {
            # extract file name from path
            file_name <- basename(vector_path)
            # load vector data 
            vector <- vect(vector_path)
            # crop vector data  to extent
            cropped_vector <- crop(vector, extent_obj)
            # save it 
            writeVector(cropped_vector, file.path(output_path, file_name), overwrite =TRUE)
            return(cropped_vector)
        }
        
    # -----------------------------------------------------------------------------------
    # 1.2) Load data and prepare data ---------------------------------------------------
        # create output directory if it doesn't exist
        dir.create(output_data_path, showWarnings=FALSE, recursive=TRUE)
        # load study area and append its path to the corresponding list, depending on the type of data (either raster or vector)
        if(study_area_type == "vector"){ 
            study_area <- vect(study_area_path)
            input_vector_data_path <- append(input_vector_data_path, study_area_path)
            } else {
            study_area <- rast(study_area_path)
            input_raster_data_path <- append(input_raster_data_path, study_area_path)
            }
    # -----------------------------------------------------------------------------------
    # 2) Generate data test set ---------------------------------------------------------
    # 2.1) Get test area extent ---------------------------------------------------------
        # determine the extent of the test area
        if(!is.null(reduction_factor)){ # by reducing the extent of study area with the reduction factor 
            extent_obj <- ext(study_area)/reduction_factor
            }
        if(!is.null(bounding_box)){ #or by using the specified bounding box
            extent_obj <- ext(bounding_box)
            } 
    # 2.2) Crop data  ------------------------------------------------------------------
        # crop all vector datas stored in the list
        vector_output <- lapply(input_vector_data_path, function(path){crop_vector(path, extent_obj, output_data_path)})
        # crop all raster datas stored in the list
        raster_output <- lapply(input_raster_data_path, function(path){crop_raster(path, extent_obj, output_data_path)})
    # -----------------------------------------------------------------------------------
}, error = function(e) {
            # If an error occured, save the error message in the error file
            writeLines(paste(Sys.time(), ": An error occurred - ", conditionMessage(e)), err)
        }
    )
# Close error file
close(err)
# Manage files progress and error tracking 
if(file.size(error_file)==0){ # if no errors have occurred 
    unlink(error_file) # deletion of the error tracking file 
    }
