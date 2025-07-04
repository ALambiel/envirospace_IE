####################################################################################
############## Fragmentation by computing Effective Mesh Size ######################
####################################################################################
#                                                                                  #
#                        PLEASE READ BEFORE RUNNING                                #
#                                                                                  #
#   This script allows you to map fragmentation, literally "effective mesh size",  # 
#   which is an index of landscape fragmentation (Jaeger, 2000) used by the FOEN.  #
#   It simply takes into account the surface area of "patches", i.e. the natural   #
#   habitat entities that make up the landscape. The map of natural environments   #
#   is classified into "barrier" environments (buildings, roads, etc.) and         #
#   habitats favorable to biodiversity (patches and corridors: forests and forest  #
#   belts, for example). Within an area studied, the effective mesh size will      #
#   therefore be equal to the total area of the zone if there are no barriers,     #
#   and equal to 0 if there are only barriers. In other words, the lower the       #
#   indicator value, the more fragmented (poor) the landscape, and the higher the  #
#   indicator value, the less fragmented (good) the landscape.
#                                                                                  #
#   The first part of the script sets up the working environment and loads and     # 
#   prepares the data required for the process. The second part of this script     #
#   uses an SQL query to select the habitat considered as barriers and to extract  #
#   them from the 5m map of natural habitat and assign them a value of 1.          #
#   Natural environments are assigned a value of 0.                                #
#   This new map is then rasterised on the model of the raster of the study area,  #
#   and extended outside this area  by Euclidean allocation. A further step        #
#   identifies all natural habitat patches (0) and assigns them an ID value. This  #
#   raster is then passed to the focal function, which uses zonal statistics to    #
#   calculate effective mesh size over a circle of 200 m radius, based on the      #
#   formula used by  the  'lsm_c_mesh' function in the 'landscapemetrics' package, #
#   which reproduces the landscape metrics calculated by the FragStat software.    #
#   The final raster is then aggregated to 25m resolution and clipped to the study #
#   area raster.                                                                   # 
#                                                                                  #
#   The results of this script are automatically saved in the "STRUCTURE"          #
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
    #  ==>   Spatial data must be in the same CRS
    #  ==>   The raster specified in 'study_area_path' will be used as a reference, meaning that the products will have the same resolution, extent and CRS as it. 

    # Study area raster (.asc or .tif) 
    study_area_path <- "path/to/data/..."
    # Map of natural habitats for the above-mentioned study area, where each polygon is assigned a habitat category (.shp) 
    habitat_map_path <- "path/to/data/..."

# Parameters ------------------------------------------------------------------------
    # Project
    # Name of the main shared project folder
    shared_directory <- "path/to/the/root/of/the/shared/folder"
    # Specify the name of an existing project or choose your new project's name
    # Please note that if you enter an existing project name, previously calculated results for this indicator may be overwritten.  
    version <- "version name"
    # Name of the pillar
    pillar_name <- "STRUCTURE" 
    # Name of the indicator
    indicator_name <- "FRAGMENTATION"
    # Give a short description of the indicator
    description <- paste0("Landscape fragmentation index (Jaeger, 2000) used by the FOEN.", "\n",
                          "It takes into account the surface area of the ‘patches’, i.e. the natural habitat entities that make up the landscape.", "\n",
                          "The lower the value of the indicator, the more fragmented the landscape (equal to 0 if there are only barriers).")

    # Datas and computing parameters 
    # CRS of your projet, e.g. to which your data are
    mycrs <- "EPSG:2056"
    # SQL query to apply at 'habitat_map_path'. Selection criterion for environments considered as barriers (assigned value: 1). Such query must be entered in the list as a single element('OR' operator between the elements of the list)  
    selection_criteria <- list(
        "CODE_MN %in% c(101, 102, 107, 218, 901, 902, 906, 910, 911, 913, 924, 1007, 1008)", #OR
        "CODE_MN > 10000", #OR
        "CODE_MN %in% c(903, 909, 908) & PROV == 13")
    # You can also select elements (SELECTED) of ‘habitat_map_path’ as barriers only if they are adjacent to certain element (SELECTOR). 
    adjacency_criteria <- list(
        SELECTED="CODE_MN %in% c(912, 4444)", # selection of categories 912 and 4444 ...
        SELECTOR="CODE_MN == 901") # ... only if they are adjacent to polygons of category 901

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
        required_packages <- c("terra", # spatial data analysis
                               "raster", # spatial data analysis
                               "FNN", # Fast Nearest Neighbor Search Algorithms and Applications
                               "dplyr", # Tool for working with data frame
                               "sf", # Support for simple features, a standardized way to encode spatial vector data
                               "purrr", # set of tools for working with functions and vectors
                               "foreach", # Support for the foreach looping construct
                               "doParallel") # Provides a parallel backend for the %dopar% function using the parallel package
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

    
        # apply_sql(x, sql, adjacent) = applies selection expressions in 'sql' and 'adjacent' (if provided) to an sf object 'x'. The expressions in 'sql' are used to select elements of 'x' based on their values, while the expression in 'adjacent' adds a spatial criterion.  
        apply_sql <- function(x, sql, adjacent = NULL) {
            # Selection of x elements according to given criteria
            res_list <-  map(sql, ~ filter(x, !!rlang::parse_expr(.x)))
            # combining results
            res_combined <- bind_rows(res_list) %>%
            # Remove any duplicates
                distinct()
            
            # If an adjacency criterion is provided, apply it
            if (!is.null(adjacent)) {
                # Extract the elements corresponding to the selection elements
                selector <- x %>% filter(., !!rlang::parse_expr(adjacent$SELECTOR))
                # Extract the items to be selected on the basis of their value and their position in relation to the selector
                selected <- x %>% filter(., !!rlang::parse_expr(adjacent$SELECTED))
                # Obtain indexes of elements that are well adjacent to the selector
                index <- which(sapply(st_touches(selected, selector), length) > 0)
                # Extract them
                adj_selected <- selected[index, ] 
                
                # Combine with the items selected in the first part
                res_combined <- bind_rows(res_combined, adj_selected) %>%
                # Remove any duplicates
                    distinct()
                }
            return(res_combined)
            }
    

        # effective_mesh_size(values, threshold, resolution) = calculates the effective mesh size (see here: https://rdrr.io/github/r-spatialecology/landscapemetrics/man/lsm_c_mesh.html) for a vector of patch IDs 'values' from a patch map raster. It computes the surface area of each patch based on the 'resolution' of the raster map, ignoring zero values. The function checks that the proportion of NA values in 'values' does not exceed 'threshold'. This function is intended for use in the focal() function.
        effective_mesh_size <- function(values, threshold, resolution){
            # Check if the proportion of NA values exceeds the threshold; if so, return NA
            if (sum(is.na(values))/length(values) >= threshold) {
                return(NA)
                }
            # Remove 0 as they do not correspond to patches ID
            patches_ID <- values[values != 0] 
            # Count the number of pixels in each patch
            count_patches <- table(patches_ID)
            # Calculate the area of each patch
            area_patches <- count_patches*resolution^2
            # Calculate the total area
            total_area <- length(values)*resolution^2
            # Calculate the effective mesh size
            mesh_value <- (sum(area_patches^2) / total_area) * (1 / 10000)
            return(mesh_value)
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
        work_directory <- file.path(shared_directory, "OUTPUTS","INDICATORS")
        # Folder for the many intermediate results, which can be deleted at the end
        scratch_folder <- file.path(work_directory, pillar_name, indicator_name, version, "scratch")
        # Directory for final outputs
        output_folder <- file.path(work_directory, pillar_name, indicator_name, version)
    
        # Create the project folder and all its subfolders (if needed)
        dir.create(work_directory, showWarnings = FALSE, recursive = TRUE)
        dir.create(scratch_folder, showWarnings = FALSE, recursive = TRUE)
        dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
    
        # Initialising the tracking file 
        tracking_file <- paste0(version, "_", script_name, ".txt")
        write(paste(Sys.time(), "RUNNING ..."), tracking_file, append = TRUE)
        # Initialising the metadata text file
        info <- c(paste0("Version : ", version, "\n",
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
        habitat_map_vector <- st_read(habitat_map_path)
        st_crs(habitat_map_vector) <- mycrs
    
        end_time <- Sys.time()
        time_loading <- difftime(end_time, start_time, units="mins")
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(), "done"), tracking_file, append=TRUE)
        # Complete metadata file with inputs' informations 
        info <- c(info, paste0("INPUTS : ", "\n",
                               " * Study area: ", study_area_path, "\n",
                               " * Natural habitat map: ", habitat_map_path, "\n",
                               " * Loading time: ", time_loading, " min","\n\n"))
    # -----------------------------------------------------------------------------------
    # 2) Compute fragmentation ----------------------------------------------------------
    # 2.1) Extraction of barrier habitats -----------------------------------------------
        write(paste(Sys.time(), "EXTRACTION OF BARRIER HABITATS"), tracking_file, append=TRUE)
        start_time <- Sys.time()
    
        # Create a new object by selecting the 'habitat_map_vector' categories corresponding to barriers according to SQL query
        selected_barrier_map <- apply_sql(x = habitat_map_vector, sql = selection_criteria, adjacent = adjacency_criteria)
    
        end_time <- Sys.time()
        time_extraction <- difftime(end_time, start_time, units = "mins")
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(), "done"), tracking_file, append=TRUE)
        # Complete metadata file 
        info <- c(info, paste0("BARRIER SELECTION", "\n",
                               " * Selection criteria applied: ", "\n"))
        info <- c(info, paste0(selection_criteria, "\n"))
        info <- c(info, paste0(adjacency_criteria$SELECTED, " if adjacent to ", adjacency_criteria$SELECTOR, "\n"))
        info <- c(info, paste0(" * Computation time: ", time_extraction, " min", "\n\n"))
    # 2.2) Rasterize the reclassified habitat map ---------------------------------------
        write(paste(Sys.time(), "RASTERISING"), tracking_file, append=TRUE)
        start_time <- Sys.time()
    
        # Rasterize reclassified map with the study area raster as a canva (same resolution, extent and CRS)
        habitat_map_raster <- rasterize(selected_barrier_map, study_area_raster, background = 0) # by setting background to 0, we create a binary raster where barriers habitat = 1 and patch habitat = 0
        habitat_map_raster <- mask(habitat_map_raster, study_area_raster)
        
        # Extend the rasterized barriers map by Euclidean allocation to avoid edge effect
        extended_habitat_map <- extend_by_euclidean_allocation(habitat_map_raster, buffer_size)
        
        # Save temporary layers in scratchfold
        writeRaster(habitat_map_raster, file.path(scratch_folder, "mn_fragmentation.tif"), overwrite = TRUE)
        writeRaster(extended_habitat_map, file.path(scratch_folder, paste0("mn_fragmentation_buff", buffer_size, "m.tif")), overwrite = TRUE)

        # Free memory
        rm(habitat_map_raster)
        gc()

        end_time <- Sys.time()
        time_rasterizing <- difftime(end_time, start_time, units = "mins")
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(), "done"), tracking_file, append=TRUE)
        # Complete metadata file 
        info <- c(info, paste0("RASTERISATION AND EXTENSION BY EUCLIDEAN ALLOCATION", "\n",
                               " * Buffer size:", buffer_size, " m", "\n",
                               " * Computation time: ", time_rasterizing, " min", "\n\n"))
    # 2.3) Compute focal statistics ---------------------------------------------------------
        write(paste(Sys.time(), "COMPUTE ZONAL STATISTICS"), tracking_file, append=TRUE)
        start_time <- Sys.time()
    
        # Create a mooving circular window
        circular_window <- focalMat(extended_habitat_map, circular_window_size, "circle")
        # Attribute the same weight (1) for all cells in the window, those outside are not considered (NA)
        circular_window <- ifelse(circular_window == 0, NA, 1)
    
        # We want to obtain a map of natural patches, where each cell has a patch ID value. The cells corresponding to the barriers must have a value of 0, which is not assigned as a patch ID.
        # Creation of a 'mask' raster where the natural patches are NA and the barriers are 0. 
        mx <- matrix(c(0, NA,
                       1, 0), ncol=2, byrow=TRUE) 
        mask_raster <- classify(extended_habitat_map, mx)

        # Prepare patches map by setting non-zero values to NA
        patches_map <- extended_habitat_map
        patches_map[values(patches_map != 0)] <- NA 
        # Identify and assign patch ID values to natural habitat patches
        patches_map <- patches(patches_map, direction = 8) 
        # Combine the patches map with the mask raster
        patches_map <- sum(patches_map, mask_raster, na.rm=TRUE)
        # Free memory 
        rm(mask_raster)

        # Seperate habitat map raster in blocks to speed up MESH computation 
        blocks <- blockSize(raster(patches_map))
    
        # Calculate overlap size
        overlap <- round(circular_window_size / res(patches_map)[1])+1
    
        # Create list of raster blocks
        blocks_list <- lapply(1:blocks$n, function(i){
            create_block(i, raster(patches_map), mycrs, blocks, overlap)
        })

        # Detect number of cores and register parallel backend
        num_cores <- detectCores()-1
        registerDoParallel(num_cores)

        # Compute effective mesh size for each block in parallel
        results <- foreach(block = blocks_list, .combine = 'c', .packages = c('terra', 'raster')) %dopar% {
            focal(block, w = circular_window, fun = function(x){
                effective_mesh_size(x, na_threshold, res(patches_map)[1])}, na.rm=FALSE)
            }
        # Stop the parallel backend
        stopImplicitCluster()
    
        # Combine rasters from each block into a single raster by calling 'mosaic' function and apply the 'max' function in the overlaping area
        focal_raster <- do.call(mosaic, c(results, fun = max, na.rm = TRUE))
        focal_raster <- rast(focal_raster)

        # Save temporary layers in scratch folder
        writeRaster(focal_raster, file.path(scratch_folder, "mn_fragmentation_focal_raw.tif"), overwrite = TRUE)
    
        end_time <- Sys.time()
        time_focal_stat <- difftime(end_time, start_time, units = "mins")
    # -----------------------------------------------------------------------------------
        # Update progress tracking file
        write(paste(Sys.time(), "done"), tracking_file, append=TRUE)
        # Complete metadata file with computing and outputs' informations
        info <- c(info, paste0("ZONAL STATISTICS", "\n",
                               " * Applied function: ", "\n"))
        info <- c(info, deparse(effective_mesh_size))
        info <- c(info, paste0(" * Circular moving window size: ", circular_window_size, " m", "\n",
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
        writeRaster(final_raster, file.path(output_folder, paste0("mn_fragmentation_focal_", res(final_raster)[1], "m.tif")), overwrite = TRUE)
        if (!is.null(do_aggregate)) {
          writeRaster(aggregated_raster, file.path(output_folder, paste0("mn_fragmentation_focal_", unique(res(final_raster)*do_aggregate), "m.tif")), overwrite = TRUE)
        }

        end_time <- Sys.time()
        time_finalising <- difftime(end_time, start_time, units="mins")
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(), "done"), tracking_file, append=TRUE)
        # Complete metadata file with computing and outputs' informations
        info <- c(info, paste0("OUTPUTS ", "\n",
                               " * Final layers saved in the: ", output_folder, "\n",
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
