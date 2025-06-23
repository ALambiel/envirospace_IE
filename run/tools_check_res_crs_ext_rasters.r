####################################################################################
##################### Check and correct CRS/EXT and resolution #####################
####################################################################################
#                                                                                  #
#                        PLEASE READ BEFORE RUNNING                                #
#                                                                                  #
#   This script allows you to check that the crs, extent and resolution of the     #
#   data contained in the 'infold' folder are identical to those of 'ref', correct #
#   them if necessary and save the aligned versions in the 'outfold' folder.       #
#                                                                                  # 
#   Make sure you select the kernel associated with the rspatial environment       #
#   when you run the script.                                                       # 
#                                                                                  #
####################################################################################
####################################################################################


# Paths ----------------------------------------------------------------------------
    infold <- "/home/path/to/input/folder"
    outfold <- "/home/path/to/output/folder"
    ref <- "/home/path/to/reference/raster.tif"

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
        # crs_check(ref, r)= compares CRS of a raster 'r' and a reference raster 'ref' and projects 'r' into the CRS of 'ref' if necessary
        crs_check <- function(ref, r) {
            if (!crs(ref) == crs(r)) {
                r <- project(r, crs(ref))
                }
            return(r)
            }
        # ext_res_check(ref, r)= compares the resolution and extent of a raster 'r' with those of a reference raster 'ref' and resamples 'r' to 'ref' if necessary
        ext_res_check <- function(ref, r) {
            if (!all(res(ref) == res(r)) || !ext(ref) == ext(r)) {
                r <- resample(r, ref, method="bilinear")
                }
            return(r)
            }
        # check_and_correct_raster(ref, r)= corrects 'r' by using the functions crs_check() and ext_res_check()
        # The function returns a list, with each sub-element containing the raster ‘r’ (corrected if necessary) and the verification tags. 
        check_and_correct_raster <- function(ref, r) {
            # get the name of the raster
            nam <- tools::file_path_sans_ext(basename(sources(r)))
            # check and correct CRS
            r <- crs_check(ref, r)
            # check and correct extent and resolution
            r <- ext_res_check(ref, r)
            # restore the correct name to the raster, potentially lost during resampling 
            varnames(r) <- nam
            # check raster properties to check that it has been corrected
            crs_match <- crs(ref) == crs(r)
            ext_match <- ext(ref) == ext(r)
            res_match <- all(res(ref) == res(r))
            all_match <- crs_match && ext_match && res_match
            return(list(
                r = r,
                crs_match = crs_match,
                ext_match = ext_match,
                res_match = res_match,
                all_match = all_match
                ))
            }
    # -----------------------------------------------------------------------------------
    # 1.2) Load data --------------------------------------------------------------------
        # read reference raster 
        rref <- rast(ref)
        # read path to each raster to be corrected
        rpath <- list.files(infold, full.names=TRUE)
        # if needed, create 'outfold'
        dir.create(outfold)
    # -----------------------------------------------------------------------------------
    # 2) Check and correct rasters ------------------------------------------------------
        # apply 'check_and_correct_raster()' function to each raster
        check <- lapply(rpath, function(path){
            r <- rast(path)
            r_clean <- check_and_correct_raster(ref=rref, r=r)
        })
        # extrat and save the corrected raster 
        res <- lapply(check, function(results){
            if(results$all_match==TRUE){
                rfinal <- results$r
                nam <- varnames(rfinal)
                writeRaster(rfinal, file.path(outfold, paste0(nam, ".tif")), overwrite=TRUE)
                return(rfinal)
                } else {
                return(cat("La couche n'est pas corrigée."))
                  }
        })
    print(res)
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
