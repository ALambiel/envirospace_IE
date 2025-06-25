####################################################################################
########################### Natural Habitat selection ##############################
####################################################################################
#                                                                                  #
#                          PLEASE READ BEFORE RUNNING                              #
#                                                                                  #
#   This script uses a habitat map to produce a raster for each habitat categories # 
#   and its zone of influence.                                                     # 
#                                                                                  #
#   The first part of the script sets up the working environment and loads and     # 
#   prepares the data required for the process. The second part of the script      #
#   is used to reclassify the habitat map categories into groups of natural        #
#   habitat of interest. The script then creates a new raster for each of the      #
#   habitats of interest, firstly by extracting it from the rasterised map, and    #
#   then by applying a zonal statistic to establish a zone of influence for the    #
#   habitat. Finally, each layer is saved according to the priority status         #
#   accorded to the habitat. 
#                                                                                  #
#   The results of this script are automatically saved in the "COMPOSITION"        #
#   subfolder of the project folder. Metadata file is automatically written        #
#   in output folder.                                                              #
#   For the script to run correctly, you only need to modify sections 'Paths' and  #
#   'Parameters' below.                                                            #
#   You don't need to make any other changes.                                      #
#                                                                                  # 
#   Make sure you select the kernel associated with the rspatial environment when  #
#   you run the script.                                                            # 
#                                                                                  # 
####################################################################################
####################################################################################


# Paths ----------------------------------------------------------------------------
    #  !! Make sure that the data you enter meets the following requirements : 
    #  ==>   Spatial data must share the same CRS (and for rasters the same resolution)
    #  ==>   The table of correspondence between the habitat categories defined in the 
    #        map of natural habitat and the categories you wish to retain in the analysis 
    #        must contain at least three columns: 
    #         1) Reference categories, THE SAME as those found in the habitat map and which you wish to reclassify 
    #         2) New categories, which allows you to group together certain reference categories under a single habitat 
    #         3) Protection status for each new categories: either 'LR' (red list) or 'NLR' (non-red list) 
    #  ==>   Respect the specified data extension

    # Study area raster (.asc or .tif)
    study_area_path <- "path/to/data/..."
    # Map of natural habitats for the above-mentioned study area, where each polygon is assigned a habitat category (.shp) 
    habitat_map_path <- "path/to/data/..."
    # Correspondence table between the categories of the habitat map and the new ones + their conservation status (.csv, sep=",") 
    corr_table_path <- "path/to/data/..."  

# Parameters ------------------------------------------------------------------------
    # Project
    # path to the main shared project folder
    shared_directory <- "path/to/the/root/of/the/shared/folder"
    # Specify the name of an existing project or choose your new project's name
    # Please note that if you enter an existing project name, previously calculated results for this indicator may be overwritten.  
    version <- "version name" 
    # Name of the pillar 
    pillar_name <- "COMPOSITION" 
    # Name of the indicator
    indicator_name <- "HABITAT"
    # Give a short descrition of the indicator
    description <- "Natural and semi-natural habitats suitable for the development of plant and animal species."

    # Datas and computing parameters
    # CRS of your projet, e.g. to which your data are
    mycrs <- "EPSG:2056" 
    # Name of the variable corresponding to the categories of habitat in the 'habitat_map_path' object 
    col_map_ref <- "CODE_MN"
    # Name of the column in 'corr_table_path' corresponding to the reference categories
    col_table_ref <- "CODE_MN" 
    # Name of the column in 'corr_table_path' indicating the new categories for reclassifying the column with the reference categories
    col_table_new <- "CAT_RECLASS" 
    # Name of the column in 'corr_table_path' indicating the priority statu
    col_prio <- "PRIO" 
    # A list of categories that should not be included in the analysis 
    excluded_categories <- c("Autres surfaces dures", "Surfaces dures", "Bâtiments", "Chemins", "Chemins imperméables", "Routes", "Sols nus sans végétation", "Voies ferrées", "Routes - Bâtiments", "Sols et substrats nus", "NONMN", " ", "", 4444, 901, 1000, 903, 905, 906, NA)
    # Function to apply for zonal statistics
    zonal_stat_function <- "sum"
    # Range of influence of habitats, e.g. radius size of the circular moving window to be applied when calculating focal statistics (in metres). 
    circular_window_size <- 100 

    # Clean up options
    # Do you want to delete the contents of the "scratch" folder at the end of the calculation? 'YES' or 'NO'
    # The temporary -non-compiled- results are saved here.
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
        required_packages <- c("terra", # Spatial data analysis
                               "stringr", # Common string opperations 
                               "dplyr") # Tool for working with data frame
        # Executing install & load for each package
        sapply(required_packages, install.load.package)

        # Functions
        # clean_characters(x) = remove from 'x' all accents and special symbols, and convert it to lowercase
        clean_characters <- function(x) {
            # If 'x' is numeric, do nothing, if else ... 
            if(!is.numeric(x)){
                # Remove accents and convert to lowercase
                x <- tolower(iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT"))
                # Remove special characters using regex
                x <- gsub("[^a-z0-9 ]", "", x)
                }
            return(x)
            }
    
        # clean_obj(obj, col) = applies the clean_characters() function to variables from 'obj' (neither data frame or spatvector). Names of selected variables must be indicate in the vector 'col'. The function return a cleaned object (same class as the input).
        clean_obj <- function(obj, col) {
            # Check if obj is a dataframe
            if (inherits(obj, "data.frame")) {
                # Check if all specified columns exist
                check_cols <- setdiff(col, colnames(obj))
                if (length(check_cols) > 0) {
                    return(NULL)
                }
                # Apply cleaning function to specified columns
                obj_clean <- obj %>%
                mutate(across(all_of(col), ~ clean_characters(.)))
                }
                
            else if (inherits(obj, "SpatVector")) {
                # Check if all specified columns exist
                check_cols <- setdiff(col, names(obj))
                if (length(check_cols) > 0) {
                    return(NULL)
                    }
                # Apply cleaning function to specified columns
                obj_clean <- obj
                val <- values(obj_clean[col])
                newval <- val %>%
                mutate(across(all_of(col), ~ clean_characters(.)))
                for(i in col){
                    values(obj_clean[,i]) <- newval[i]
                    }
                } else {
                return(NULL)
                }
            return(obj_clean)
            }
    # -----------------------------------------------------------------------------------
    # 1.3) Working directories ----------------------------------------------------------
        # Path to the working directory
        work_directory <- file.path(shared_directory, "OUTPUTS","INDICATORS")
        # Folder for the many intermediate results, which can be deleted at the end
        scratch_folder <- file.path(work_directory, pillar_name, indicator_name, version, "scratch")
        # Directory for final outputs
        output_folder <- file.path(work_directory, pillar_name, indicator_name, version)
        # Subfolder for 'Red List' habitat (LR)
        LR <- file.path(output_folder, "LR") 
        # Subfolder for 'Non Red List' habitat (NLR)
        NLR <-file.path(output_folder, "NLR") 
            
        # Create the project folder and all its subfolders (if needed)
        dir.create(work_directory, showWarnings = FALSE, recursive = TRUE)
        dir.create(scratch_folder, showWarnings = FALSE, recursive = TRUE)
        dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
        dir.create(LR, showWarnings = FALSE, recursive = TRUE)
        dir.create(NLR, showWarnings = FALSE, recursive = TRUE)
            
        # Initialising the tracking file 
        tracking_file <- paste0(version, "_", script_name, ".txt")
        writeLines(paste(Sys.time(), "RUNNING ..."), con=tracking_file)
        # Initialising the metadata text file
        info <- c(paste0("Version : ", version, "\n",
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
            
        # Study area
        study_area_raster <- rast(study_area_path)
        crs(study_area_raster) <- mycrs
            
        # Habitat map 
        habitat_map_vector <- vect(habitat_map_path)
        crs(habitat_map_vector) <- mycrs
            
        # Correspondence table
        corr_table <- read.csv(corr_table_path, sep=",")
        # Make sure to consider the targeted column 'col_map_ref', 'col_table_ref' , 'col_table_new' and 'col_prio' provided by user
        # To do so, select columns as specified in parameter section and rename them
        # In 'habitat_map_vector'
        p1 <- match(col_map_ref,names(habitat_map_vector)) 
        names(habitat_map_vector)[p1] <- "REF" 
        # In 'corr_table'
        p2 <- match(c(col_table_ref, col_table_new, col_prio), names(corr_table)) #same in the correspondance table
        names(corr_table)[p2] <- c("REF", "RECLASS", "PRIO")
        
        # Clean categories in variables indicated by the user by applying 'clean_obj' function
        habitat_map_vector <- clean_obj(habitat_map_vector, c("REF"))
        corr_table <- clean_obj(corr_table, c("REF", "RECLASS"))
        # Clean also the vector of excluded categories 
        excluded_cat <- clean_characters(excluded_categories)
        
        end_time <- Sys.time()
        time_loading <- difftime(end_time, start_time, units="mins")
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(),"done"), tracking_file, append=TRUE)
        # Complete metadata file with inputs' informations 
        info <- c(info, paste0("INPUTS : ", "\n",
                                       " * Study area: ", study_area_path, "\n",
                                       " * Natural habitat map: ", habitat_map_path, "\n",
                                       " * Reference MN categories: ", col_map_ref, "\n",
                                       " * Correspondence table: ", corr_table_path, "\n",
                                       " * Categories not included: ","\n"))
        info <- c(info, paste0("  '", excluded_categories, "'","\n"))
        info <- c(info, paste0(" * Loading time : ", time_loading, " min", "\n\n"))
    # -----------------------------------------------------------------------------------
    # 2) Computing ----------------------------------------------------------------------
    # 2.1) Reclassification of habitat map categories -----------------------------------
        write(paste(Sys.time(), "RECLASSIFICATION"), tracking_file, append=TRUE)
        start_time <- Sys.time()
            
        # Reclassification of habitat map categories in habitat of interest
        habitat_map_vector$RECLASS <- corr_table$RECLASS[match(habitat_map_vector$REF, corr_table$REF)]
            
        end_time <- Sys.time()
        time_reclassify <- difftime(end_time, start_time, units="mins")
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(),"done"), tracking_file, append=TRUE)
        # Complete metadata file
        info <- c(info, paste0("RECLASSIFICATION ", "\n",
                               " * Matching categories: ", "\n",
                               "  ",  col_map_ref, " <== ", col_table_new, "\n"))
        info <- c(info, paste0("  '", corr_table[,"REF"],"' <== '", corr_table[,"RECLASS"],"'", "\n"))
        info <- c(info, paste0(" * Computation time: ", time_reclassify, " min", "\n\n"))
   # 2.2) Creating layers for each habitat ----------------------------------------------
        write(paste(Sys.time(), "CREATING LAYERS FOR EACH HABITAT"), tracking_file, append=TRUE)
        start_time <- Sys.time()
            
        # Get the list of new categories
        habitat_list <- unique(habitat_map_vector$RECLASS)

        result <- lapply(habitat_list, function(category){
            # Do not consider categories specified in the 'excluded_cat' vector
            if(!(category %in% excluded_cat)){ 
                # Select only poylgons corresponding to the 'category'
                selected_habitat <- habitat_map_vector[habitat_map_vector$RECLASS==category]
                
                # Rasterize it with the study area raster as a canva (same resolution, extent and CRS). 
                # As 'cover=TRUE', the pixel's value is now the fraction of a cell that is covered by the polygon
                selected_habitat_raster <- rasterize(selected_habitat, study_area_raster, background=0, touches=TRUE, cover=TRUE)
                # Mulitply the values of selected_habitat_raster by 100 to reach entier value 
                selected_habitat_raster <- selected_habitat_raster*100
                
                # Compute focal statistic 
                # Create a mooving circular window
                circular_window <- focalMat(selected_habitat_raster, circular_window_size, "circle") 
                # Attribute the same weight (1) for all cells in the window, those outside are not considered (NA)
                circular_window <- ifelse(circular_window == 0, NA, 1) 
                # Apply the selected function (zonal_stat_function) to the cells in the window, and assign the resulting value to the central cell.
                focal_raster <- focal(selected_habitat_raster, circular_window, fun = zonal_stat_function)

                # Mask the resulting raster with the study area 
                final_raster <- mask(focal_raster, study_area_raster)

                # Save the final raster in the correct folder accordinf to protection status
                # Get the file name and replace potential space in habitat's name to avoid future error in Zonation
                nam <- str_replace_all(category, "[[:space:]]|'|-", "_")
                # Get the priority statut associated to the habitat
                priority <- unique(corr_table[corr_table$RECLASS == category, "PRIO"])
                # Save the layer 
                if(priority == "LR"){
                    writeRaster(final_raster,filename=paste0(LR, "/", nam, "_focal.tif"), overwrite=TRUE)
                    }else{
                    writeRaster(final_raster,filename=paste0(NLR, "/", nam, "_focal.tif"), overwrite=TRUE)
                    }
                # save raw result (e.g. the habitat layer without the focal statistic) in the scratchfold
                writeRaster(selected_habitat_raster,filename=paste0(scratch_folder,"/", nam, ".tif"), overwrite=TRUE)
                
                # Update the progress tracking file 
                progress <- round(which(habitat_list==category)/length(habitat_list)*100, 2)
                write(paste0("> Progress : ", progress," % -- ",category,": ok."), tracking_file, append=TRUE)
                # Return the resolution and CRS of the layer produced to check that all the layers have been produced correctly. 
                return(list("RES"=res(final_raster)[1],"CRS"=crs(final_raster, describe=TRUE)[1]))
                 }
            })
    
        # Filters out all the null values from 'result'
        result <- Filter(Negate(is.null), result)
        # Check that all outputs have the same resolution and CRS
        check <- unique(result)
        if(length(check)==1){
            output_resolution <- result[[1]]$RES
            output_crs <- result[[1]]$CRS
            } else {
                write(paste0("Warning ! It seems that there are different resolutions or CRS in the outputs : ", result), error_file, append=TRUE)
            }

        end_time <- Sys.time()
        time_computing <- difftime(end_time, start_time, units="mins")
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(), "done"), tracking_file, append=TRUE)
        # Complete metadata file with computing and outputs' informations
            info <- c(info, paste0("EXTRACTION AND ZONAL STATISTICS", "\n",
                                       " * Applied function: ", zonal_stat_function, "\n",
                                       " * Circular window size: ", circular_window_size, " m", "\n",
                                       " * Comuptation time: ", time_computing, " min", "\n\n"))
        info <- c(info, paste0("OUTPUTS ", "\n",
                               " * Final layers saved in: ", output_folder, "\n",
                               " * Resolution : ", output_resolution, " m", "\n",
                               " * CRS : ", output_crs, "\n",
                               " * Number of layers produced: ", length(result), "\n",
                               " * RL habitats :", "\n"))
        info <- c(info, paste0("  ", list.files(LR), "\n"))
        info <- c(info, paste0(" * NRL habitats:", "\n"))
        info <- c(info, paste0("  ", list.files(NLR), "\n"))
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
if(scratch_to_trash=="YES"){ 
    unlink(scratch_folder, recursive =TRUE) # deletion of files in the scratch folder
    }
# Manage files progress and error tracking 
if(track_to_trash=="YES" & file.size(error_file)==0){ # if the user has chosen "YES" to allow deletion AND that no errors have occurred 
    unlink(tracking_file) # Delete the progress tracking file
    unlink(error_file) # Delete the error tracking file 
    }
