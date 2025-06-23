####################################################################################
################## Species distribution modelling with MaxEnt ######################
####################################################################################
#                                                                                  #
#                        PLEASE READ BEFORE RUNNING                                #
#                                                                                  #
#   This script allows you modelling the distribution of a list of given species   #
#   by calling up the Maxent software via R.                                       #
#   It takes as input a list of species and the coordinates of the observation     #
#   points, and relates these presences to environmental variables called          #
#   "predictors", in order to produce habitat suitability predictions. See here:   #
#    https://biodiversityinformatics.amnh.org/open_source/maxent/                  # 
#                                                                                  #
#   The first part of the script sets up the working environment and loads and     # 
#   prepares the data required for the process. In particular, it ensures that     #
#   all the data is correctly aligned with the reference raster.                   #
#   The second part of this script is used to generate a sampling bias raster,     #
#   which is used to produce pseudo-absences. These are used, along with the       # 
#   observations, to train MaXENT models. For each species, several models are     #
#   trained and then used to predict habitat suitability predictions for the       #
#   species by applying the model onto the predictors. All outputs are evaluated.  #
#   The  final predictions for each species are obtained by averaging the          #
#   predictions obtained for each model, and the evaluation metrics for all the    #
#   models are stored in a table. Final rasters is clipped to the study area       #
#   raster.                                                                        # 
#                                                                                  #
#   The results of this script are automatically saved in the "COMPOSITION"        #
#   subfolder of the project folder. Metadata file is automatically written in     #
#   output folder.                                                                 # 
#   For the script to run correctly, you only need to modify sections 'Paths' and  #
#   'Parameters' below.                                                            #
#   You don't need to make any other changes.                                      #
#                                                                                  # 
#   Make sure you select the kernel associated with the mainenv.yml environment    #
#   when you run the script.                                                       # 
#                                                                                  #
####################################################################################
####################################################################################

# Paths ----------------------------------------------------------------------------
    #  !! Make sure that the data you enter meets the following requirements : 
    #  ==>  Spatial data must share the same CRS and the same resolution for rasters 
    #  ==> 'NAME' = scientific name of species
    #  ==> 'PRIO' = protection status: either 'LR' (red list) or 'NLR' (non-red list)
    #  ==> 'X' and 'Y' = GPS coordinates in the CRS project
    #  ==> Respect the specified data extension

    # Study area raster (.asc or .tif)
    study_area_path <- "path/to/data/..." 
    # Directory with environmental predictors rasters (.asc or .tif)
    predictors_path <- "path/to/data/..."
    # Table with focus species and their conservation status (.csv, sep=","), with at least 'NAME' and 'PRIO' column
    species_list_path <- "path/to/data/..." 
    # Species observations (.csv, sep=",") with at least 'NAME', 'X' and 'Y' column
    species_obs_path <- "path/to/data/..." 
    
# Parameters ------------------------------------------------------------------------
    # Project
    # Name of the main shared project folder
    shared_directory <- "path/to/the/root/of/the/shared/folder"
    # Specify the name of an existing project or choose your new project's name
    # Please note that if you enter an existing project name, previously calculated results for this indicator may be overwritten.
    project_name <- "version name" 
    # Name of the pillar 
    pillar_name <- "COMPOSITION"
    # Name of the indicator 
    indicator_name <- "FAUNA"
    # Give a short descrition of the indicator
    description <- "Distribution of species (fauna) using biotic and abiotic criteria."

    # Datas and computing parameters 
    # CRS in which your data are
    mycrs <- "EPSG:2056" 
    
    # Modelling parameters  
    # Minimum observation threshold for a species. 
    minimum_obs <- 20 
    # Number of pseudo absence to generate
    num_abs <- 10000 
    # Number of repetitions for fitting model
    num_fit <- 10 
    # Proportion of species observations to be set aside for validation  
    # For example, set to 0.25 to keep 25% of the data as a test set, and use the other 75% to train the model
    test_proportion <- 0.25
    # Parameters for the MaXENT function
    # More details: https://www.rdocumentation.org/packages/dismo/versions/1.3-14/topics/maxent 
    maxent_parameters <- c(
        'defaultprevalence=0.5',
        'betamultiplier=2',
        'pictures=FALSE',
        'linear=TRUE',
        'quadratic=TRUE',
        'product=TRUE',
        'threshold=FALSE',
        'hinge=TRUE',
        'threads=2',
        'responsecurves=FALSE',
        'jackknife=FALSE')
    # Specify whether one of the environmental variables used as a predictor is categorical, by noting here the name of the raster as it appears in the predictors folder (e.g. : "name.asc")
    factors <- c()

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
        # list of packages required  
        required_packages <- c("raster", # Spatial data analysis
                               "sp", # Spatial data
                               "sf", # Support for simple features (spatial vector data)
                               "terra", # Spatial data analysis --> maxent() function
                               "dismo", # Species distribution modeling
                               "rJava", # Low-level interface to Java VM (require for maxent) 
                               "doParallel", # Provides a parallel backend for the %dopar% function
                               "foreach", #Support for the foreach looping construct 
                               "MASS", # Statistical analysis 
                               "stringr", # Common string opperations 
                               "data.table") # Enhanced data frame 
        # Executing install & load for each package
        sapply(required_packages, install.load.package)

        # Functions
        # remove_obs(df, var, n) = remove from a dataframe 'df' the rows corresponding to categories in the 'var' variable which have fewer than 'n' occurrences. The function returns the cleaned df. 
        remove_obs <- function(df, var, n){
            # Get the frequency of each category in 'df'
            freq <- table(df[[var]])
            # Select categories with a frequency greater than or equal to 'n'
            catselect <- names(which(freq >= n)) 
            # Create a new clean data frame with the rows of the selected categories
            dfclean <- df[df[[var]] %in% catselect,] 
            return(dfclean)
            }

        # align_raster(input_path, reference_raster, output_raster) = aligns a raster whose access path is given by 'input_path' with a 'reference_raster'. The raster thus aligned is then saved under 'output_path' with the same original name, and this updated path is returned. 
        align_raster <- function(input_path, reference_raster, output_path){
            raster <- rast(input_path)
            # Check if all values in the raster are integers to select the resampling method
            is_integer <- all(as.integer(values(raster)) == values(raster), na.rm = TRUE)
            # Choose the resampling method based on whether the raster has integer values
            method <- if (is_integer) "near" else "bilinear"
            # Resample the raster to match the reference raster using the chosen method
            raster_aligned <- resample(raster, reference_raster, method = method)
            # Construct the output path using the input filename
            output_path <- file.path(output_path, basename(input_path))
            # Save the resampled raster to the output path
            writeRaster(raster_aligned, output_path, overwrite = TRUE)
            return(output_path)
        }
    
        # calculate_bias(x,y,r) = estimate the sampling bias by using the 'x' and 'y' species observation coordinates, create a density raster using the kde2d (Two-Dimensional Kernel Density Estimation) function. The function returns a density raster based on the reference raster 'r'. 
        calculate_bias <- function(x, y, r) {
                # To avoid an error caused by interquartile ranges equal to 0, and therefore a non-positive bandwidth
                bandwidth_x <- ifelse(bandwidth.nrd(x) == 0, 0.1, bandwidth.nrd(x)) 
                bandwidth_y <- ifelse(bandwidth.nrd(y) == 0, 0.1, bandwidth.nrd(y))
                # Use kde2 from MASS to get a two-dimensional kernel density estimate
                dens <- kde2d(x, y, h = c(bandwidth_x, bandwidth_y),
                              n = c(nrow(r), ncol(r)), 
                              lims = c(ext(r)[1], ext(r)[2], ext(r)[3], ext(r)[4]))
                # Rasterize it
                d <- raster(dens)
                # Allign it with the reference raster
                d <- resample(d, r, method = "bilinear")
            return(d)
            }

        # split_data(df, p) = split observations in dataframe 'df' into train and test set, with respective proportions 'p' and '1-p'. Function returns a list containing 'train' and 'test' data frame.  
        split_data <- function(df, p){
            # On the sample of size 'nrow(df)' we randomly select a proportion of 'p' observations to make up the train set
            index_train <- sample(1:nrow(df), p*nrow(df))
            # Create the new object 'train' with selected obs 
            train <- df[index_train,]
            # Assign the remaining obs to the 'test' object
            test <- df[-index_train,]
            # Return a list object with 'train' and 'test' sets
            return(list(train=train, test=test))
            }

        # fit_pred_eval(sp, n presence, absence, predictor, args, factors, p, crs) = for species 'sp', trains 'n' models with maxent() by taking as input 'presence' and 'absence' data, a stack of 'predictor' rasters, and vector 'args' and 'fact' containing respectively model's parameters and any variables' name constituting predictors that are factors. Each fit is carried out with a new train and test group (separation of the data set according to the proportion 'p') by using the function split_data(). Then, evaluate and make habitat suitability predictions with each models. 'crs' is the CRS of spatial data. The function returns a list containing metrics (data.table) and predicted rasters for the 'n' models. 
        fit_pred_eval <- function (sp, n, presence, absence, predictor, args, factors, p, crs){
            results <- lapply(1:n, function(j) {
                split <- split_data(df = presence, p = p)
                train <- split$train
                test <- split$test
                fit <- maxent(x = predictor, p = train, a = absence, args = args, factors = factors)
                pred <- predict(fit, predictor)
                names(pred) <- sp
                crs(pred) <- crs
                ev <- evaluate(fit, p = test, a = absence, x = predictor)
                metrics <- data.table(NAME = sp, MODEL = j, AUC = ev@auc, COR = ev@cor, pCOR = ev@pcor)
                return(list(metrics = metrics, pred = pred))
            })
            return(results)
        }
    # -----------------------------------------------------------------------------------
    # 1.3) Working directories ----------------------------------------------------------
        # Path to the working directory
        work_directory <- file.path(shared_directory, "OUTPUTS","INDICATORS", project_name) 
        # Folder for the many intermediate results, which can be deleted at the end
        scratch_folder <- file.path(work_directory, "scratch", paste(pillar_name, indicator_name, date, sep="_"))
        # Directory for final outputs
        output_folder <- file.path(work_directory, pillar_name, indicator_name)
        # Subfolder for 'Red List' habitat (RL)
        LR <- file.path(output_folder, "RL") 
        # Subfolder for 'Non Red List' habitat (NRL)
        NLR <-file.path(output_folder, "NRL") 
    
        # Create the project folder and all its subfolders (if needed)
        dir.create(work_directory, showWarnings = FALSE, recursive = TRUE)
        dir.create(scratch_folder, showWarnings = FALSE, recursive = TRUE)
        dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
        dir.create(LR, showWarnings = FALSE, recursive = TRUE)
        dir.create(NLR, showWarnings = FALSE, recursive = TRUE)
        
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
    
        # Study area raster
        study_area_raster <- raster(study_area_path)
        crs(study_area_raster) <- mycrs
    
        # Environmental predictors
        # Get a list of all the predictors contained in 'predictors_path'
        predictors_list <- list.files(predictors_path, full.names=TRUE) 
        # Make sure to align predictors raster with study area
        predicors_aligned_list <- lapply(predictors_list, function(raster_path){
            align_raster(raster_path, rast(study_area_raster), scratch_folder)
        })
        # Create a stack of all the aligned predictors
        predictors <- stack(predicors_aligned_list)
        # Mask it with the study area raster
        predictors <- mask(predictors, study_area_raster)
        crs(predictors) <- mycrs

        # Species observations
        # Load species occurrences
        species_obs <- fread(species_obs_path, sep=",", header=TRUE)
        # Keep only occurences in 'species_obs' that are inside the study area
        species_obs <- terra::mask(vect(species_obs, geom=c("X","Y"), crs=mycrs), as.polygons(rast(study_area_raster), aggregate=TRUE, na.rm=TRUE)) 
        # Come back to a data frame
        species_obs <- as.data.table(species_obs, geom = "XY") 
        # Load list of species to focus on 
        species_list <- fread(species_list_path, sep=",", header=TRUE)
        # Keep only species which are both in 'species_obs' and 'species_list'
        species_obs <- species_obs[species_obs$NAME %in% unique(species_list[["NAME"]]),] 
        # Remove species with fewer than nlim occurences
        species_obs <- remove_obs(df = species_obs, var = "NAME", n = minimum_obs)
        # Get a final list of species names
        selected_species <- sort(unique(species_obs$NAME))
    
        end_time <- Sys.time()
        time_loading <- difftime(end_time, start_time, units="mins")
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(), "done"), tracking_file, append=TRUE)
        # Complete metadata file with inputs' informations 
        info <- c(info, paste0("INPUTS: ", "\n",
                               " * Study area: ", study_area_path, "\n",
                               " * Species observations file: ", species_obs_path, "\n",
                               " * Species list file: ", species_list_path, "\n",
                               " * Predictors list: ","\n"))
        info <- c(info, paste0("  ", predictors_list, "\n"))
        info <- c(info, paste0(" * Reading time: ", time_loading, " min", "\n\n"))
    # -----------------------------------------------------------------------------------
    # 2) Modelling ----------------------------------------------------------------------
    # 2.1) Bias file for pseudo-asbence -------------------------------------------------
        write(paste(Sys.time(), "BIAS"), tracking_file, append=TRUE)
        start_time <- Sys.time()

        # Initialise the parallelisation of the calculation process 
        # Detect the number of cores
        num_cores <- detectCores() - 1
        # Separate the list into num_cores groups, so that they can be processed in parallel. 
        data_chunks <- split(selected_species, rep(1:num_cores, length.out = length(selected_species)))
        # Register parallel backend 
        cl <- makeCluster(num_cores)
        registerDoParallel(cl)
        # Calculate the biasfile by applying the 'calculate_bias()' function to all the species and summing the rasters obtained for each of them
        bias_list <- foreach(chunk = data_chunks, .combine = 'c', .packages = c('raster', 'MASS', 'terra')) %dopar% {
            lapply(chunk, function(species_name){
                
            # Get coordinates
            obsx <- species_obs[species_obs$NAME==species_name, "x"]
            obsy <- species_obs[species_obs$NAME==species_name, "y"]

            # Compute bias 
            calculate_bias(x=obsx, y=obsy, r=study_area_raster)
            })
        }
        # Stop the parallel backend
        stopCluster(cl)

        # Combine all the bias rasters stored in 'bias_list' by first stacking and then adding them together 
        bias_staked <- do.call(stack, bias_list)
        bias_raster <- calc(bias_staked, sum, na.rm=TRUE)
        crs(bias_raster) <- mycrs

        # Mask it with study area raster
        bias_raster <- mask(bias_raster, study_area_raster)

        # Save the biasfile in the scratch folder
        writeRaster(bias_raster,filename=paste(scratch_folder,"biasfile.tif", sep="/"),format="GTiff", overwrite=TRUE) 

        # Free memory 
        length_bias <- length(bias_list)
        rm(bias_list)
    
        end_time <- Sys.time()
        time_bias <- difftime(end_time, start_time, units="mins")
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(), "done"), tracking_file, append=TRUE)
        # Complete metadata file with running time 
        info <- c(info, paste0("BIAS ", "\n",
                               " * Bias based on ", length_bias, " layers.", "\n",
                               " * Computation time: ", time_bias, " min", "\n\n"))
    # -----------------------------------------------------------------------------------
    # 2.2) Fit, predict and evaluate ----------------------------------------------------
        write(paste(Sys.time(), "FIT, PREDICT AND EVALUATE"), tracking_file, append=TRUE)
        start_time <- Sys.time()
        
        # Create pseudo-absences linked to biasfile (biasfile's values are interpreted as probability weights)
        pseudo_abs <- randomPoints(bias_raster, num_abs, prob=TRUE)
        colnames(pseudo_abs)=c("x", "y") 

        # Register parallel backend 
        cl <- makeCluster(num_cores)
        registerDoParallel(cl)
        # Run the fit_pred_eval function in parallel for each species
        all_metrics <- foreach(chunk = data_chunks, .combine='c', .packages = c('raster', 'dismo', 'terra', 'data.table', 'stringr')) %dopar% {
            lapply(chunk, function(species_name){   
                # Get informations for the species 
                # Get presence 
                presence <- species_obs[species_obs$NAME == species_name, c("x", "y")]
                # Get priority statut
                priority <- species_list[species_list$NAME == species_name, "PRIO"]
                
                # Call the fit_pred_eval function
                maxent_result <- fit_pred_eval(sp = species_name, 
                                               n = num_fit, 
                                               presence = presence, 
                                               absence = pseudo_abs, 
                                               predictor = predictors, 
                                               args = maxent_parameters,
                                               factors = factors, 
                                               p = 1-test_proportion, 
                                               crs = mycrs)

                # Combine the prediction rasters obtained for the num_fit models into a single raster, by averaging. 
                pred_stack <- stack(lapply(maxent_result, function(x) x$pred))
                final_pred <- calc(pred_stack, mean)

                # Store CRS and resolution
                crs_info <- crs(final_pred)
                res_info <- res(final_pred)

                # Combine in a data frame the metrics for each models 
                metrics <- rbindlist(lapply(maxent_result, function(x) x$metrics))
                
                # Create a folder to save the final prediction 
                # First, ensure that the file name is free of spaces
                nam <- str_replace_all(species_name, " ", "_")
                nam <- paste0(nam, ".tif")
                # Save the raster in the correct folder, according to its priority status 
                file_path <- file.path(output_folder, priority, nam)
                writeRaster(final_pred, file_path, overwrite=TRUE)

                # Return in a list the metrics' data frame, CRS and resolution 
                return(list(metrics = metrics, crs = crs_info, res = res_info))
                })
            }
        # Stop the parallel backend
        stopCluster(cl)

       # Extract CRS and resolution information
       crs_list <- lapply(all_metrics, function(x) x$crs)
       res_list <- lapply(all_metrics, function(x) x$res)
       # Check if all CRS are the same and get the unique value
       check_crs <- unique(crs_list)
       if (length(check_crs) == 1) {
           output_crs <- check_crs[[1]]
           } else {
           output_crs <- check_crs
           write("CRS values are not identical in final predictions.", tracking_file, append=TRUE)
           }
       # Check if all resolutions are the same and get the unique value
       check_resolution <- unique(res_list) 
       if (length(check_resolution) == 1) {
           output_resolution <- check_resolution[[1]][1]
           } else {
           output_resolution <- check_resolution
           write("Resolution values are not identical in final predictions.", tracking_file, append=TRUE)
       }                          

        # Combine all metrics tables stored in 'all_metrics' 
        total_metrics <- do.call(rbind, lapply(all_metrics, function(x) x$metrics))
                                               
        # Save the final data frame in the output folder
        write.csv(total_metrics, file.path(output_folder, "METRICS.csv"), row.names = FALSE) 

        end_time <- Sys.time()
        time_prediction <- difftime(end_time, start_time, units="mins")
    # -----------------------------------------------------------------------------------
        # Update the progress tracking file
        write(paste(Sys.time(), "done"), tracking_file, append=TRUE)
        # Complete metadata file with modelling parameters and running time 
        info <- c(info, paste0("MAXENT", "\n",
                               " * Applied function: "))
        info <- c(info, deparse(fit_pred_eval))
        info <- c(info, paste0(" * Maxent parameters used: "))
        info <- c(info, paste0("  ", maxent_parameters, "\n"))
        info <- c(info, paste0(" * Number of pseudo-absences generated: ", num_abs, "\n",
                               " * Minimum number of observations to make the prediction: ", minimum_obs, "\n",
                               " * Percentage of data retained for evaluation: ", test_proportion*100, " %", "\n",
                               " * Training and evaluation of the model with ", num_fit, " repetitions.", "\n", 
                               " * Computation time: ", time_prediction, " min", "\n\n"))
        info <- c(info, paste0("OUTPUTS", "\n",
                               " * Final predictions saved in: ", output_folder, "\n",
                               " * Resolution: ", output_resolution, " m", "\n",
                               " * CRS : ", output_crs, "\n",
                               " * Total number of species modelled: ", length(selected_species), "\n",
                               " * RL species:"))
        info <- c(info, paste0("  ", list.files(LR), "\n"))
        info <- c(info, paste0(" * NRL species:"))
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
write(paste0("##### Total duration : ", total_run_time, " min #####"), file.path(output_folder, "METADATA.txt"), append =TRUE) 
                                           
# Close error file
close(err)
                                 
# Manage scratch folder content
if(scratch_to_trash=="YES"){ # if the user has chosen "YES" to allow the contents of the sracth folder to be deleted 
    unlink(scratch_folder, recursive =TRUE) # deletion of files in the scratch folder
    }
# Manage files progress and error tracking 
if(track_to_trash=="YES" & file.size(error_file)==0){ # if the user has chosen "YES" to allow deletion AND that no errors have occurred 
    unlink(tracking_file) # deletion of the progress tracking file
    unlink(error_file) # deletion of the error tracking file 
    }
