####################################################################################
############################# Green spaces Diversity ###############################
####################################################################################
#                                                                                  #
#                        PLEASE READ BEFORE RUNNING                                #
#                                                                                  #
#   This script allows you to map diversity (in terms of number of different       #
#   green environments) and equipartition of green spaces in a pixel's sur-        #
#   roundings. It is based on the "diversity" attribute (ranging from 1 to 8),     #
#   which separates habitat categories of the habitat map into different types of  #
#   green spaces, and uses the Shannon diversity index formula to asses the        #
#   diversity of green environments in a pixel's surroundings.                     # 
#                                                                                  #
#   The first part of the script sets up the working environment and loads and     # 
#   prepares the data required for the process. The second part of the script      #
#   assigns diversity value to each pixel in the 5m resolution habitat map, based  #
#   on the attribute table provided. The map of habitat is then rasterised on the  #
#   model of the raster of the study area, and extended outside this area  by      #
#   Euclidean allocation. The final final diversity raster is calculated  by Focal #
#   Statistics over a circle of 200 m radius, based on the formula used by  the    #
#   'lsm_l_shdi' function in the 'landscapemetrics' package, which reproduces the  #
#   landscape metrics calculated by the FragStat software. The final raster  is    #
#   then aggregated to 25m resolution and clipped to the study area raster.        #
#                                                                                  #
#   The results of this script are automatically saved in the "STRUCTURE"          #
#   subfolder of the project folder. Metadata file is automatically written        #
#   in output folder.                                                              #
#   For the script to run correctly, you only need to modify sections 'Paths' and  #
#   'Parameters' below.                                                            #
#   You don't need to make any other changes.                                      #
#                                                                                  # 
#   Make sure you select the kernel associated with the mainenv environment when   #
#   you run the script.                                                            # 
#                                                                                  #
####################################################################################
####################################################################################


# Paths ----------------------------------------------------------------------------
    #  !! make sure that the data you enter meets the following requirements : 
    #  ==>   Spatial data must be in the same CRS
    #  ==>   The raster specified in 'study_area_path' will be used as a reference, meaning that the products will have the same resolution, extent and CRS as it. An aggregated version of the final raster could be produced.

    # Study area raster (.asc or .tif) 
    study_area_path <- "path/to/data/..."
    # Map of natural habitats for the above-mentioned study area, where each polygon is assigned a habitat category (.shp) 
    habitat_map_path <- "path/to/data/..."
    # Attribute table linking all categories of the habitat map with a green environment type (.csv, sep=",") 
    mn_attribute_table_path <- "path/to/data/..."  

# Parameters ------------------------------------------------------------------------
    # Project
    # Name of the main shared project folder
    shared_directory <- "path/to/the/root/of/the/shared/folder"
    # Specify the name of an existing project or choose your new project's name
    # Please note that if you enter an existing project name, previously calculated results for this indicator may be overwritten.  
    project_name <- "version name"
    # Name of the pillar
    pillar_name <- "STRUCTURE" 
    # Name of the indicator
    indicator_name <- "DIVERSITY"
    # Give a short descrition of the indicator
    description <- paste0("Shannon diversity index for areas of selected natural habitats. It takes into account the abundance and equipartition of habitats within a given perimeter.", "\n", 
                          "It is based on the approach used by the Fragstat software to calculate diversity as a landscape metric.")

    # Datas and computing parameters 
    # CRS of your projet, e.g. to which your data are
    mycrs <- "EPSG:2056"
    # Name of the variable corresponding to the categories of habitat in the 'habitat_map_path' object
    col_map_ref <- "CODE_MN"
    # Name of the column in 'mn_attribute_table_path' corresponding to the reference categories
    col_table_ref <- "CODE_MN"
    # Name of the column in 'mn_attribute_table_path' corresponding to the naturalness index value
    col_table_index <- "DIVERSITY"
    # All the values that can be taken by the index mentioned under 'col_table_index' and that correspond to the different types of green spaces. 
    diversity_classes <- c(1:8)
    # Background value, e.g. Value representing non-green spaces 
    background_value <- 99
    # Buffer size used to extend the raster using Euclidean allocation (in meters)
    buffer_size  <- 50
    # Radius size of the circular moving window to be applied when calculating focal statistics (in metres) 
    circular_window_size <- 200
    # Maximum ratio of NA present in a window for the diversity calculation to be carried out
    na_threshold <- 0.5

    # Aggregation options
    # Do you want to compute an aggregated version of the final raster, at a lower resolution of 'study_area_path' reference raster ? 
    # If not, set 'do_aggragte' as NULL
    # If yes, assign to 'do_aggragte' the value of the factor by which the resolution of 'study_area_path' should be multiplied. 
    # For example, if 'study_area_path' has a resolution of 5m and you want to aggregate the raster at a resolution of 25m, doaggr should be set to 5.
    do_aggregate <- 5
    # Function to apply for aggregation 
    aggregation_function <- "mean"

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
        # list of packages required  
        required_packages <- c("terra", # spatial data analysis
                               "raster", # spatial data analysis
                               "FNN", # Fast Nearest Neighbor Search Algorithms and Applications
                               "doParallel", # Provides a parallel backend for the %dopar% function using the parallel package
                               "foreach") # Support for the foreach looping construct
        # Executing install & load for each package
        sapply(required_packages, install.load.package)

        # Functions
        # extend_by_euclidean_allocation(r, b) = Extend the  raster 'r' by applying a buffer 'b' (in metres) per Eulidian allocation, i.e. by assigning to the new pixels the value of the nearest pixel from 'r'.
        extend_by_euclidean_allocation <- function(r, b){
            # Extend the 'r' by a buffer size 'b' with NA values
            R1 <-  extend(r, ext(r)+b, fill=NA) 
            # Calculate the Euclidean distance for each NA cell in the extend raster
            R2 <- distance(R1)
            # Initiate the result raster with the same structure as the extended raster
            R3 <- rast(R1)
            # Fill the final raster based on the following conditions
            
            # Condition 1 : fill R3 where R2 == 0 with values from R1
            R3[R2==0] <- R1[R2==0]
            
            # Condition 2 : fill R3 where 0 < R2 <= b with the values of the nearest non-NA pixel in R1
            cells_to_fill <- which(values(R2) > 0 & values(R2) <= b)
            # Get the coordinates of all cells in 'R1'
            coords <- xyFromCell(R1, 1:ncell(R1))
            # Get the indices of non-NA cells in 'R1'
            non_na_indices <- which(!is.na(values(R1)))
            # Get their coordinates
            non_na_coords <- coords[non_na_indices, ]
            # Get their values
            non_na_values <- values(R1)[non_na_indices]
            # Find the nearest non-NA neighbor for each cell in R3 that needs filling
            knn_result <- get.knnx(data = non_na_coords, query = coords[cells_to_fill, ], k=1) # k=1 means we looking for 1 neighbor only 
            # Extract indices of the nearest neighbor found
            nearest_indices <- knn_result$nn.index
            # Fill the identified cells in 'R3' with the values of their nearest non-NA neighbord from R1
            R3[cells_to_fill] <- non_na_values[nearest_indices]
            
            # Condition 3 : fill R3 where R2 > b with NA
            R3[R2 > b] <- NA
            
            # Return the final raster
            return(R3)
            }
    

        # shannon_diversity(values, all_classes, threshold) = calculates the Shannon diversity index (based on: https://rdrr.io/github/r-spatialecology/landscapemetrics/man/lsm_l_shdi.html) for class values in 'values'. It includes all possible classes specified in 'all_classes' and checks that the proportion of NA values in 'values' does not exceed 'threshold'. This function is intended for use in the focal() function.
        shannon_diversity <- function(values, all_classes, threshold) {
            # Check if the proportion of NA values exceeds the threshold; if so, return NA
            if (sum(is.na(values)) / length(values) >= threshold) {
                return(NA)
                }
            # Count occurrences for each class
            class_counts <- table(factor(values, levels = all_classes))
            # Calculate proportions of each class
            proportions <- class_counts / sum(class_counts)
            # Calculate the Shannon diversity index
            shdi <- -sum(proportions * log(proportions), na.rm=TRUE)
            return(shdi)
            }
    
        # create_block(i, raster, crs, blocks, overlap) = defines the extent needed to extract sub-raster/block 'i' from 'raster', based on 'blocks' from the blockSize function (which suggests chunk sizes and corresponding row numbers for processing Raster* objects in chunks). 'crs' is the coordinate reference system of 'raster', and 'overlap' ensures that the blocks overlap to avoid edge effects. 
        create_block <- function(i, raster, crs, blocks, overlap) {
            # Calculate the starting row for the block, ensuring it doesn't go below the first row
            start_row <- max(blocks$row[i] - overlap, 1)
            # Calculate the ending row for the block, ensuring it doesn't go beyond the last row
            end_row <- min(blocks$row[i] + blocks$nrows[i] - 1 + overlap, nrow(raster))
            # Define the extent of the block using the start and end rows and the full width of the raster
            block_extent <- extent(raster, start_row, end_row, 1, ncol(raster))
            # Crop the raster to the defined block extent
            block_raster <- crop(raster, block_extent)
            # Assign the specified coordinate reference system to the cropped raster block
            crs(block_raster) <- crs
            # Return the sub-raster/block 
            return(block_raster)
            }
    # -----------------------------------------------------------------------------------
    # 1.3) Working directories ----------------------------------------------------------
        # Path to the working directory
        work_directory <- file.path(shared_directory, "OUTPUTS", "INDICATORS", project_name)
        # Folder for the many intermediate results, which can be deleted at the end
        scratch_folder <- file.path(work_directory, "scratch", paste(pillar_name, indicator_name, date, sep = "_"))
        # Directory for final outputs
        output_folder <- file.path(work_directory, pillar_name, indicator_name)
    
        # Create the project folder and all its subfolders (if needed)
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
    
        # Habitat map 
        habitat_map_vector <- vect(habitat_map_path)
        crs(habitat_map_vector) <- mycrs
    
        # Attribute table
        mn_attribute_table <- read.csv(mn_attribute_table_path, sep=",")
        # Make sure to consider the targeted column 'col_map_ref', 'col_table_ref' and 'col_table_index' provided by user
        # To do so, select columns as specified in parameter section and rename them
        # In 'habitat_map_vector'
        p1 <- match(col_map_ref,names(habitat_map_vector))
        names(habitat_map_vector)[p1] <- "REF"
        # In 'mn_attribute_table'
        p2 <- match(c(col_table_ref, col_table_index), names(mn_attribute_table))
        names(mn_attribute_table)[p2] <- c("REF", "INDEX")
    
        end_time <- Sys.time()
        time_loading <- difftime(end_time, start_time, units="mins")
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(), "done"), tracking_file, append=TRUE)
        # Complete metadata file with inputs' informations 
        info <- c(info, paste0("INPUTS", "\n",
                               " * Study area: ", study_area_path, "\n",
                               " * Natural habitat map: ", habitat_map_path, "\n",
                               " * Attribute table: ", mn_attribute_table_path, "\n",
                               " * Reading time: ", time_loading, " min", "\n\n"))
    # -----------------------------------------------------------------------------------
    # 2) Compute diversity -------------------------------------------------------------
    # 2.1) Assigning diversity values --------------------------------------------------
        write(paste(Sys.time(), "ASSIGNING DIVERSITY VALUES"), tracking_file, append=TRUE)
        start_time <- Sys.time()
    
        # Reclassification of habitat map categories in type of green spaces 
        habitat_map_vector$RECLASS <- mn_attribute_table$INDEX[match(habitat_map_vector$REF, mn_attribute_table$REF)]
    
        end_time <- Sys.time()
        time_reclassify <- difftime(end_time, start_time, units="mins")
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(), "done"), tracking_file, append=TRUE)
        # Complete metadata file 
        info <- c(info, paste0("ASSIGNING NATURALNESS VALUES", "\n",
                               " * Number of diversity classes: ", length(diversity_classes), "\n",
                               " * Background value: ", background_value, "\n",
                               " * Matching table: ", "\n",
                               "  ",  col_map_ref, " <== ", col_table_index, "\n"))
        info <- c(info, paste0("  '", mn_attribute_table[,"REF"],"' <== '", mn_attribute_table[,"INDEX"],"'", "\n"))
        info <- c(info, paste0(" * Computation time: ", time_reclassify, " min", "\n"))
    # 2.2) Rasterize the reclassified habitat map ---------------------------------------
        write(paste(Sys.time(), "RASTERISING"), tracking_file, append=TRUE)
        start_time <- Sys.time()
    
        # Extract only informations about diversity and save it in another vector layer
        diversity_vector <- habitat_map_vector["RECLASS"]
        # Rasterize it with the study area raster as a canva (same resolution, extent and CRS)
        diversity_raster <- rasterize(diversity_vector, study_area_raster, field = "RECLASS")
    
        # Extend the rasterized diversity map by Euclidean allocation to avoid edge effect
        extended_diversity_raster <- extend_by_euclidean_allocation(diversity_raster, buffer_size)
    
        # Save temporary layers in scratch folder
        writeRaster(diversity_raster, file.path(scratch_folder, "mn_diversity.tif"), overwrite=TRUE)
        writeRaster(extended_diversity_raster, file.path(scratch_folder, paste0("mn_diversity_buff", buffer_size, "m.tif")), overwrite=TRUE)

        # Free memory
        rm(diversity_raster)
        gc()

        end_time <- Sys.time()
        time_rasterizing <- difftime(end_time, start_time, units = "mins")
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(), "done"), tracking_file, append=TRUE)
        # Complete metadata file 
        info <- c(info, paste0("RASTERISATION AND EXTENSION BY EUCLIDEAN ALLOCATION", "\n",
                               " * Buffer size: ", buffer_size, " m", "\n",
                               " * Computation time: ", time_rasterizing, " min", "\n\n"))
    # 2.3) Compute focal statistics ---------------------------------------------------------
        write(paste(Sys.time(), "COMPUTE ZONAL STATISTICS"), tracking_file, append=TRUE)
        start_time <- Sys.time()
    
        # Create a mooving circular window
        circular_window <- focalMat(extended_diversity_raster, circular_window_size, "circle")
        # Attribute the same weight (1) for all cells in the window, those outside are not considered (NA)
        circular_window <- ifelse(circular_window == 0, NA, 1)
    
        # Get suggested chunk sizes to crop the diversity raster in sub-rasters
        # This give us the suggested 'row' numbers at which to start the blocks for reading and writing + 'nrows' the number of rows in each block and the total number 'n' of blocks 
        blocks <- blockSize(raster(extended_diversity_raster))
        # Get the overlap size according to the radius size of circular window and the resolution of the raster 
        overlap <- round(circular_window_size / res(extended_diversity_raster)[1])+1
        # Get a list of raster blocks with the 'create_block()' function 
        blocks_list <- lapply(1:blocks$n, function(i){
            create_block(i, raster(extended_diversity_raster), mycrs, blocks, overlap)
        })

        # Initialise parallel computing 
        # Detect the number of cores available and leave one aside
        num_cores <- detectCores()-1
        # Register the parallel backend with the foreach package.
        registerDoParallel(num_cores)
        # Apply the 'shannon_diversity()' funciton through the 'focal()' function, for each block in parallel
        results <- foreach(block = blocks_list, .combine = 'c', .packages = c('terra', 'raster')) %dopar% {
            focal(block, w = circular_window, fun = function(x){shannon_diversity(x, diversity_classes, na_threshold)}, na.rm=FALSE)
            }
        # Close the cluster
        stopImplicitCluster()
    
        # Combine rasters from each block into a single raster by calling 'mosaic' function and apply the 'max' function in the overlaping area
        focal_raster <- do.call(mosaic, c(results, fun = max, na.rm = TRUE))
        focal_raster <- rast(focal_raster)

        # Save temporary layers in scratch folder
        writeRaster(focal_raster, file.path(scratch_folder, "mn_diversity_focal_raw.tif"), overwrite = TRUE)

        # Free memory
        rm(extended_diversity_raster)
        gc()
    
        end_time <- Sys.time()
        time_focal_stat <- difftime(end_time, start_time, units = "mins")
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(), "done"), tracking_file, append=TRUE)
        # Complete metadata file with computing and outputs' informations
        info <- c(info, paste0("ZONAL STATISTICS    ", "\n",
                               " * aPPLIED FUNCTION: ", "\n"))
        info <- c(info, deparse(shannon_diversity))
        info <- c(info, paste0(" * Circular window size: ", circular_window_size, " m", "\n",
                               " * Maximum ratio of NA accepted in the window: ", na_threshold, "\n",
                               " * Computation time", time_focal_stat, " min", "\n\n"))
    # 2.4) Prepare final result ---------------------------------------------------------
        write(paste(Sys.time(), "PREPARE FINAL RESULT"), tracking_file, append=TRUE)
        start_time <- Sys.time()
    
        # Come back to an unbuffered raster
        final_raster <- crop(focal_raster, study_area_raster)
        final_raster <- mask(final_raster, study_area_raster)
    
        # If needed, aggregate the final raster to a lower resolution
        if (!is.null(do_aggregate)) {
          aggregated_raster <- aggregate(final_raster, fact = do_aggregate, fun = aggregation_function)
        }
        
        # Save outputs
        writeRaster(final_raster, file.path(output_folder, paste0("mn_diversity_focal_", res(final_raster)[1], "m.tif")), overwrite = TRUE)
        if (!is.null(do_aggregate)) {
          writeRaster(aggregated_raster, file.path(output_folder, paste0("mn_diversity_focal_", unique(res(final_raster)*do_aggregate), "m.tif")), overwrite = TRUE)
        }
    
        end_time <- Sys.time()
        time_finalising <- difftime(end_time, start_time, units="mins")
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(), "done"), tracking_file, append=TRUE)
        # Complete metadata file with computing and outputs' informations
        info <- c(info, paste0("OUTPUTS", "\n",
                               " * Final layers saved in: ", output_folder, "\n",
                               " * Resolutions : ", res(final_raster)[1], " et ", if (!is.null(do_aggregate)) res(aggregated_raster)[1], " m", "\n", 
                               " * Aggregation function: ", aggregation_function, "\n",
                               " * CRS : ", crs(final_raster, describe = TRUE)[1], "\n",
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
if (scratch_to_trash == "YES") {
    unlink(scratch_folder, recursive = TRUE)  # Delete the scratch folder
    }
# Manage files progress and error tracking
if (track_to_trash == "YES" & file.size(error_file) == 0) { # if the user has chosen "YES" to allow deletion AND that no errors have occurred
    unlink(tracking_file)  # Delete the progress tracking file
    unlink(error_file)  # Delete the error tracking file
}
