####################################################################################
#################################### Title #########################################
####################################################################################
#                                                                                  #
#                        PLEASE READ BEFORE RUNNING                                #
#                                                                                  #
#   This script allows you XXXX                                                    #
#                                                                                  #
#   Describe the main steps of the script                                          # 
#                                                                                  #
#   The results of this script are automatically saved in the "XXX"                #
#   subfolder of the project folder. Metadata file is automatically written in the #
#   output folder.                                                                 #
#   For the script to run correctly, you only need to modify sections 'Paths' and  #
#   'Parameters' below.                                                            #
#   You don't need to make any other changes.                                      #
#                                                                                  # 
#   Make sure you select the kernel associated with the XXXXXXX environment        #
#   when you run the script.                                                       # 
#                                                                                  #
####################################################################################
####################################################################################

# Paths ----------------------------------------------------------------------------
    #  !! Make sure that the data you enter meets the following requirements : 
    #  ==>   requirement 1

    study_area_path <- "path/to/data/..."
# Parameters ------------------------------------------------------------------------
    # Project
    # Path to the main shared project folder
    shared_directory <- "path/to/the/root/of/the/shared/folder"
    # Specify the name of an existing project or choose your new project's name
    # Please note that if you enter an existing project name, previously calculated results for this indicator may be overwritten.  
    version <- ""
    # Name of the pillar (either "COMPOSITION", "STRUCTURE", "FONCTION", "SE")
    pillar_name <- "" 
    # Name of the indicator (In uppercase and without spaces/special characters !!)
    indicator_name <- ""
    # Give a short descrition of the indicator
    description <- ""

    # Datas and computing parameters
    # CRS of your projet, e.g. to which your data are
    mycrs <- "EPSG:XXXX"

    # Add potential options (e.g. aggregation of final rasters to a lower resolution, etc.).
    
    # Clean up options
    # Do you want to delete the contents of the "scratch" folder at the end of the calculation? 'YES' or 'NO'
    # Temporary -non compiled- files were saved here. 
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
error_file <- paste0(version, "_", script_name, "_err.out")
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
        # list of packages required  
        required_packages <- c() # please detail the main goal of each packages 
        # Executing install & load for each package
        sapply(required_packages, install.load.package)

        # Functions
        # add here custom functions used in the script
    # -----------------------------------------------------------------------------------
    # 1.3) Working directories ----------------------------------------------------------
       # Path to the working directory
        work_directory <- file.path(shared_directory, "OUTPUTS","INDICATORS") 
        # Folder for the many intermediate results, which can be deleted at the end
        scratch_folder <- file.path(work_directory, pillar_name, indicator_name, version, "scratch")
        # Directory for final outputs
        output_folder <- file.path(work_directory, pillar_name, indicator_name, version)
    
        # Create directories if they don't exist
        dir.create(work_directory, showWarnings = FALSE, recursive = TRUE)
        dir.create(scratch_folder, showWarnings = FALSE, recursive = TRUE)
        dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
    
        # Initialising the tracking file 
        tracking_file <- paste0(version, "_", script_name, ".txt")
        write(paste(Sys.time(), "RUNNING ..."), tracking_file, append = TRUE)
        # Initialising the metadata text file
        info <- c(paste0("Version name : ", version, "\n",
                         "Date : ", date, "\n",
                         "User : ", user, "\n\n",
                         pillar_name, " - ", indicator_name, "\n",
                         description, "\n\n"))
    # -----------------------------------------------------------------------------------
    # 1.4) Load and prepare data --------------------------------------------------------
        write(paste(Sys.time(),"LOAD AND PREPARE DATA"), tracking_file, append=TRUE)
        start_time <- Sys.time()
    
        # Coordinates reference systeme 
        mycrs <- crs(mycrs)
    
        end_time <- Sys.time()
        time_loading <- difftime(end_time, start_time, units="mins")
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(),"done"), tracking_file, append=TRUE)
        # Complete metadata file with inputs' informations 
        info <- c(info, paste0("INPUTS", "\n",
                               # list here all the input data
                               "  Loading time : ", time_loading, " min", "\n\n"))
    # -----------------------------------------------------------------------------------
    # 2) Title --------------------------------------------------------------------------
    # 2.X) Subtitle X -------------------------------------------------------------------
        write(paste(Sys.time(), "SUBTITLE"), tracking_file, append=TRUE)
        start_time <- Sys.time()
    
        end_time <- Sys_time()
        time_XXXX <- difftime(end_time, start_time, units="mins")
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(),"done"), tracking_file, append=TRUE)
        # Complete metadata file
        info <- c(info, paste0("SUBTITLE", "\n",
                               # all useful informations about process (function, parameters, etc.)
                               # add computing time
                              ))
    # 2.X) Subtitle X -------------------------------------------------------------------
        write(paste(Sys.time(), "SUBTITLE"), tracking_file, append=TRUE)
        start_time <- Sys.time()
    
        end_time <- Sys_time()
        time_XXXX <- difftime(end_time, start_time, units="mins")
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(),"done"), tracking_file, append=TRUE)
        # Complete metadata file
        info <- c(info, paste0("SUBTITLE", "\n",
                               # all useful informations about process (function, parameters, etc.)
                               # add computing time
                              ))
     # 2.X) Prepare final result --------------------------------------------------------
        write(paste(Sys.time(), "PREPARE FINAL RESULT"), tracking_file, append = TRUE)
        start_time <- Sys.time()
    
        end_time <- Sys_time()
        time_finalising <- difftime(end_time, start_time, units="mins")
        # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(), "done"), tracking_file, append = TRUE)
        # Complete metadata file with ouptuts informations
        info <- c(info, paste0("OUTPUTS", "\n",
                               " * Final layers saved in: ", output_folder, "\n",
                               " * Resolutions : ", , " m", "\n",
                               " * CRS : ", , 
                               " * Completion time:", time_finalising, " min", "\n\n"))
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
