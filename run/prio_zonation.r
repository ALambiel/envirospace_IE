####################################################################################
######################## Running Zonation prioritization project  ##################
####################################################################################
#                                                                                  #
#                        PLEASE READ BEFORE RUNNING                                #  
# Description:                                                                     #
# This script combines environmental or planning indicators                        #
# (referred to as "pillars") into a hierarchical prioritization                    #
# using the Zonation 5 software.                                                   #
#                                                                                  #
# Requirements:                                                                    #
# - Zonation 5 must be installed (AppImage format)                                 #
# - Installation info: https://github.com/Alambiel/envirospace_IE                  #
#                                                                                  #
# Script structure:                                                                #
# 1. Define input indicators (pillars)                                             #
# 2. Set Zonation parameters                                                       #
#    (More info: https://zonationteam.github.io/Zonation5/)                        #
# 3. Choose weighting options for prioritization                                   #
# 4. Generate settings file for the Zonation project                               #
# 5. Launch Zonation and generate outputs                                          #
#                                                                                  #
# Notes:                                                                           #
# - Make sure the Zonation 5 AppImage is accessible and executable                 #
# - Outputs include prioritization maps and report files                           #
#                                                                                  #
#   The results of this script are automatically saved in the "zonation_output"    #
#   subfolder of the project folder. Metadata file is automatically written in the #
#   output folder.                                                                 #
#   For the script to run correctly, you only need to modify sections 'Paths' and  #
#   'Parameters' below.                                                            #
#   You don't need to make any other changes.                                      #
#                                                                                  # 
#   Make sure you select the kernel associated with the XXXXXXX.yml environment    #
#   when you run the script.                                                       # 
#                                                                                  #
####################################################################################
####################################################################################

# Paths ----------------------------------------------------------------------------
    # Path to Zonation5 AppImage (z5 executable)
    z5_exe <- "/your/path/to/squashfs-root/usr/bin/z5"

    # Path to the main shared project folder
    shared_directory <- "/home/..."

    # Subdirectories containing raster files 
    dir_structure <- "path/to/folder/STRUCTURE"
    dir_es <- "path/to/folder/ES"
    dir_comp_lr <- "path/to/COMPOSITION/.../LR"
    dir_comp_nlr <- "path/to/COMPOSITION/.../NLR"

    # Mask file path
    mask_file<-""

# Parameters ------------------------------------------------------------------------
    # Specify the name of an existing project or choose your new project's name
    # Please note that if you enter an existing project name, previously calculated results for this indicator may be overwritten.  
    version <- ""

    # Datas and computing parameters
    # CRS of your projet, e.g. to which your data are
    mycrs <- "EPSG:2056"

    # Clean up options
    # Do you want to delete the contents of the "scratch" folder at the end of the calculation? 'YES' or 'NO'
    # Temporary -non compiled- files were saved here. 
    scratch_to_trash <- "YES"
    # Do you want to delete the progress and error files generated while the script is running when it finishes? 'YES' or 'NO'
    # n.b.: If "YES" but an error occurs, the two files will not be deleted in all cases to allow debugging.  
    track_to_trash <- "YES"

# Zonation specific parameters: 

    # Zonation marginal loss rule: one of "ABF", "CAZ1", "CAZ2", "CAZMAX"
    MLR <- "CAZ2"

    # Use Weights? (TRUE/FALSE)
    use_weights <- TRUE

    # Use Weight Groups? (TRUE/FALSE)
    use_weight_groups <- TRUE

    # Use habitat mask? (TRUE/FALSE)
    use_habitat_mask <- FALSE

    # Weights per pillar (structure, ES, RL, NRL) 
    weights <- c(1, 1, 2, 1)

    # Weight group IDs (for higher order aggregation) 
    weight_groups <- c(1, 2, 3, 3)

    # Weight group values (e.g. for group 1,2,3, a value of 20,30,50, total 100)
    weight_groups_values<- c(20, 30, 50)


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
        # Functions
        get_files <- function(path, pattern) {
          list.files(path, recursive = TRUE, full.names = TRUE, pattern = pattern)
        }
    
        # add here custom functions used in the script
    # -----------------------------------------------------------------------------------
    # 1.3) Working directories ----------------------------------------------------------
       # Path to the working directory
        work_directory <- file.path(shared_directory, "OUTPUTS", "PRIORITIZATION", version) 
        # Folder for the many intermediate results, which can be deleted at the end
        scratch_folder <- file.path(work_directory, "scratch")
        # Directory for final outputs
        output_folder <- file.path(work_directory,"zonation_output")
        # Directory for zonation settings files
        settings_folder<-file.path(work_directory,"zonation_settings")
    
        # Create directories if they don't exist
        dir.create(work_directory, showWarnings = FALSE, recursive = TRUE)
        dir.create(scratch_folder, showWarnings = FALSE, recursive = TRUE)
        dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
        dir.create(settings_folder, showWarnings = FALSE, recursive = TRUE)
    
        # Initialising the tracking file 
        tracking_file <- paste0(version, "_", script_name, ".txt")
        write(paste(Sys.time(), "RUNNING ..."), tracking_file, append = TRUE)
        # Initialising the metadata text file
        info <- c(paste0("Projet : ", version, "\n",
                         "Date : ", date, "\n",
                         "User : ", user, "\n\n"))
    # -----------------------------------------------------------------------------------
    # 1.4) Load and prepare data --------------------------------------------------------
        write(paste(Sys.time(),"LOAD AND PREPARE DATA"), tracking_file, append=TRUE)
        start_time <- Sys.time()
    
    # Output setting files paths
        settings_file <- file.path(settings_folder, "settings.z5")  
        all_files_txt <- file.path(settings_folder, "all_files.txt")
        weight_groups_txt <- file.path(settings_folder, "weight_group.txt")
    
    # reading the rasters in the subdirectories
    # Structure, ES, Composition Red list, Composition not red list. only tif
        Str <- get_files(dir_structure, "_25m.tif$")
        ES <- get_files(dir_es, ".tif$")
        CompLR <- get_files(dir_comp_lr, ".tif$")
        CompNLR <- get_files(dir_comp_nlr, ".tif$")
    
    # Generating a list of the pillar files
        piliers <- list(Str, ES, CompLR, CompNLR)
    
        end_time <- Sys.time()
        time_loading <- difftime(end_time, start_time, units="mins")
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(),"done"), tracking_file, append=TRUE)
        # Complete metadata file with inputs' informations 
        info <- c(info, paste0("INPUTS", "\n",
                               # list here all the input data
                               "  Temps de lecture : ", time_loading, " min", "\n\n"))
    # -----------------------------------------------------------------------------------
    # 2) GENERATE setting files  --------------------------------------------------------
    # 2.1) Settings.z5  -----------------------------------------------------------------
       write(paste(Sys.time(), "Settings.z5"), tracking_file, append=TRUE)
        start_time <- Sys.time()

    settings_lines <- c(
      paste("feature list file =", basename(all_files_txt))
    )

    if (use_weight_groups) {
      settings_lines <- c(settings_lines, paste("weight groups file =", basename(weight_groups_txt)))
    }

    if (use_habitat_mask) {
      settings_lines <- c(settings_lines, paste("analysis area mask layer =", habitat_mask_path))
    }

    # Export settings.z5 file
    writeLines(settings_lines, settings_file)

        end_time <- Sys.time()
        time_XXXX <- difftime(end_time, start_time, units="mins")
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(),"done"), tracking_file, append=TRUE)
        # Complete metadata file
         info <- c(info, paste0("ZONATION SETTINGS FILE", "\n",
                       " - Settings file: ", settings_file, "\n",
                       " - Marginal loss rule: ", MLR, "\n",
                       " - Habitat mask: ", ifelse(use_habitat_mask, habitat_mask_path, "None"), "\n",
                       " - Generation time: ", time_XXXX, " min\n\n"))
    
    # 2.2) Generate list of all input files with weights  -------------------------------
        write(paste(Sys.time(), "Generated settings.z5 file"), tracking_file, append=TRUE)
        start_time <- Sys.time()

        df_list <- mapply(function(files, w, g) {
          data.frame(weight = w, wgrp = g, filename = files)
        }, piliers, weights, weight_groups, SIMPLIFY = FALSE)

        final_df <- do.call(rbind, df_list)

        write.table(final_df, file = all_files_txt, row.names = FALSE, quote = FALSE)
    
        end_time <- Sys.time()
        time_XXXX <- difftime(end_time, start_time, units="mins")
    
        info <- c(info, paste0("FEATURE LIST FILE", "\n",
                       " - Total raster inputs: ", nrow(final_df), "\n",
                       " - Weights used: ", paste(weights, collapse = ", "), "\n",
                       " - Generation time: ", time_XXXX, " min\n\n"))
 # 2.3) Generate weight group file  -------------------------------
        write(paste(Sys.time(), "Generated all_files.txt"), tracking_file, append=TRUE)
        start_time <- Sys.time()

        
        if (use_weight_groups) {
          wg_data <- data.frame(
            col1 = unique(weight_groups),
            col2 = weight_groups_values[seq_along(unique(weight_groups))] 
            )
          write.table(wg_data, file = weight_groups_txt, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
        }
    
        end_time <- Sys.time()
        time_XXXX <- difftime(end_time, start_time, units="mins")
    
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(),"done"), tracking_file, append=TRUE)
        # Complete metadata file
        info <- c(info, paste0("WEIGHT GROUP FILE", "\n",
                       " - Groups: ", paste(unique(weight_groups), collapse = ", "), "\n",
                       " - Group values: ", paste(weight_groups_values, collapse = ", "), "\n",
                       " - Generation time: ", time_XXXX, " min\n\n"))
     # 3) Run Zonation analysis --------------------------------------------------------
        write(paste(Sys.time(), "Generated weight_groups.txt"), tracking_file, append = TRUE)
        start_time <- Sys.time()

        # Determine composite flag based on user entry
        weight_flag <- if (use_weight_groups) {
          if (use_habitat_mask) "-Wa" else "-W"
        } else if (use_weights) {
          if (use_habitat_mask) "-wa" else "-w"
        } else if (use_habitat_mask) {
          "-a"
        } else {
          NULL
        }

        args <- c(
             weight_flag,
          paste0("--mode=", MLR),
          settings_file,
          paste(work_directory,"zonation_output",sep="/")
        )
    # Lauchning zonation call
    system2(command = z5_exe, args = args)
    
        end_time <- Sys.time()
        time_finalising <- difftime(end_time, start_time, units="mins")
        # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(), "done"), tracking_file, append = TRUE)
        # Complete metadata file with ouptuts informations
        info <- c(info, paste0("OUTPUTS", "\n",
                       " - Final outputs saved in: ", output_folder, "\n",
                       " - Coordinate Reference System: ", mycrs, "\n",
                       " - Time to finalize: ", time_finalising, " min\n\n"))
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
