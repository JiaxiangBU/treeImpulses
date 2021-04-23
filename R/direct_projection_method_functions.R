# File: direct_projection_method_functions.R
# Author: Tyler Pike
# Date: 2/5/2020
# Note(s): Houses functions to calculate different types of Jorda (2005) style direct projections

#-------------------------------------------------------------#
# Function to produce IRF                                  ####
#-------------------------------------------------------------#
#' Estimate local projections
#'
#' The primary function of the treeImpulses package. Calculates local projection style impulse responses using either linear regression
#' or tree-based estimators. The function supports state-dependent and instrumental variable routines. The user may choose between
#' trees based on Athey and Wager (2019) [\url{https://arxiv.org/abs/1510.04342}] or Oprescu et al (2019) [\url{http://proceedings.mlr.press/v97/oprescu19a.html}].
#'
#' @param data Dataframe of covariate
#' @param shock String denoting variable to shock
#' @param target String denoting variable betas to collect
#' @param horizons Horizons to forecast out to
#' @param engine String declaring how to estimate responses (must be 'OLS', 'AW', or 'OSW')
#' @param confidence Alpha used for two-sided confidence intervals
#' @param lags Lags to include in regressions. Must be a deterministic number of lag or an information criterion used to chose the number of lags (AIC or BIC). Chosen based on information criteria applied to OLS - future functionality may support tree-based ic.
#' @param lags_max Maximum number of lags
#' @param standardize Standardize inputs (mean = 0, variance = 1)
#' @param states Column name used to define states
#' @param NW Newey-West correction on variance-covariance matrix
#' @param NW_lags Number of lags to use in Newey-West correction
#' @param NW_prewhite Prewhite option for Newey-West correction (see sandwich::NeweyWest function)
#' @param IV Perform a two-stage instrumental variable regression (currently only supported for trees)
#' @param treeBag Number of tree-based models to bag over
#' @param Cholesky Append the initial structural shock from a cholesky decomposition [primarly for replication purposes]
#' @param cores Number of cores available to use (must be integer) set -1 to indicate use all available resources
#' @param seed Random seed to be used (must be integer)
#'
#' @return An object of class 'IRF', a list of \code{impulseResponse} (a matrix of impulse responses and standard errors),
#' \code{engine} (the method used to estimate impulses),
#' \code{target} (the variable that was responding to a shock),
#' \code{shock} (the variable that was shocked), and
#' \code{lags} (the number of lags used in estimation). 
#' Standard errors are not reported for OSW trees-based impulse responses, only upper and lower confidence intervals.
#'
#' @export

localProjectionIRF <-
  function(data,                   # dataframe of covariate
           shock,                  # string denoting variable to shock
           target,                 # string denoting variable betas to collect
           horizons,               # horizons to forecast out to
           engine = 'OLS',         # string declaring how to estimate responses (must be 'OLS', 'AW', or 'OSW')
           confidence = 0.1,       # alpha used for two-sided confidence intervals
           # regressor specifications
           lags = 'bic',           # lags to include in regressions
                                   #   (must be a deterministic number of lag, or an information criterion used to chose the number of lags)
                                   #   (chosen based on information criteria applied to OLS or AW tree - trees only support the mse criterion)
           lags_max = 12,          # maximum number of lags
           standardize = FALSE,    # standardize inputs (mean = 0, variance = 1)
           # Threshold/State dependent parameters
           states = NULL,          # column name used to define states
           # OLS-based IRF parameters
           NW = FALSE,             # Newey-West correction on variance-covariance matrix
           NW_lags = NULL,         # number of lags to use in Newey-West correction
           NW_prewhite = NULL,     # prewhite option for Newey-West correction (see sandwich::NeweyWest function)
           # Tree-based IRF parameteres
           IV = FALSE,             # perform a two-stage instrumental variable regression (currently only supported for trees)
           treeBag = 5,            # number of tree-based models to bag over
           Cholesky = FALSE,       # append the initial structural shock from a cholesky decomposition [primarly for replication purposes]
           # technical options
           cores = NULL,           # number of cores available to use (must be integer), set -1 to indicate use all available resources
           seed  = NULL            # random seed to be used (must be integer)
          ){

    #-------------------------------------------------#
    # warnings                                     ####
    #-------------------------------------------------#
    # check that a valid engine was chosen
    if(length(intersect(c('AW','OLS','OSW'), engine)) != 1){
      errorCondition('Please enter valid forecasting engine.')
    }
    # check that the states variable is found in the data matrix
    if(!is.null(states) & length(intersect(states, colnames(data))) != 1){
      errorCondition('States variable are not found in data provided.')
    }
    # check that the shock and target variables are found in the data matrix
    if(length(intersect(c(shock,target), colnames(data))) != 2){
      errorCondition('Shock and target variables are not found in data provided.')
    }
    # check that the horizons are a positive integer
    if(as.numeric(horizons) %% 1 != 0 | as.numeric(horizons) <= 0){
      errorCondition('Horizons must be a positive integer.')
    }
    # check that the lags are a positive integer or zero
    if((as.numeric(lags) %% 1 != 0 | as.numeric(lags) < 0) & (!lags %in% c('aic','bic','mse'))){
      errorCondition('Lags must be an non-negative integer or declare an information criterion for lag selection.')
    }
    # check the forbidden variable names are not used
    if('X' %in% names(data) | 'i' %in% names(data)){
      errorCondition('X and i are forbidden variable names')
    }

    #-------------------------------------------------#
    # python                                       ####
    #-------------------------------------------------#
    if(engine == 'OSW'){
       warning('Using reticulate to access python packages. \nRemember to set your python enviroment.')
       np = reticulate::import('numpy')
       sklearn = reticulate::import('sklearn')
       econML = reticulate::import('econml.ortho_forest')
    }

    #-------------------------------------------------#
    # establish parallel backend                   ####
    #-------------------------------------------------#
    # set random seed
    if(is.null(seed)){seed = 1116}
    # set number of cores available
    if(is.null(cores)){
      print('Cores note supplied so only 1 is used.')
      cores <- 1
    }

    #-------------------------------------------------#
    # prepare data for analysis                    ####
    #-------------------------------------------------#
    # first, create the proper table and variable names
    data <- data %>%
      data.frame() %>%
      dplyr::rename(target = target,
                    states = states) %>%
      na.omit()

    # edge case when the shock is the target variable
    if(shock == target){
      shock = 'target'
    }

    # exogenous variable names
    # data.variables <- colnames(select(data, -contains('states')))
    data.variables <- 
      data %>% names %>% 
      str_subset("states", negate = TRUE)
    # browser()
    # second, standardize inputs, if desired
    if(standardize == TRUE){
    data <- data %>%
      dplyr::mutate_at(vars(-contains('states')), function(X){return((X-mean(X, na.rm = T))/sd(X, na.rm = T))})
    }

    #-------------------------------------------------#
    # prepare and generate variable lags           ####
    #-------------------------------------------------#
    # pre-determined number of lags
    if(is.numeric(lags)){
      data <- list(data) %>%
        purrr::map(
          function(X, lags.n = lags){
              recipes::recipe(~ ., data = X) %>%
              recipes::step_lag(tidyselect::everything(), -tidyselect::contains('states'), lag = 1:lags) %>%
              recipes::step_naomit(tidyselect::contains('lag')) %>%
              recipes::prep(X) %>%
              recipes::bake(X)
          }
        ) %>%
        purrr::map(.f = na.omit) %>%
        as.data.frame()

    # choose the optimal number of lags
    }else if(lags %in% c('bic', 'aic', 'mse')){
      data <- list(data) %>%
        purrr::map(
          function(X, lags.m = lags_max){
            # information storage
            accuracy.temp <- vector(length = lags.m)

            # iterate through potential lags
            for(lags.temp in 1:lags.m){

              # prep lags
              data.temp <- X %>%
                recipes::recipe(~ ., data = X) %>%
                recipes::step_lag(tidyselect::everything(), -tidyselect::contains('states'), lag = 1:lags.temp) %>%
                recipes::step_naomit(tidyselect::contains('lag')) %>%
                recipes::prep(X) %>%
                recipes::bake(X)

              # model
              if(lags %in% c('aic', 'bic', 'mse') & engine == 'OLS'){
                model.temp <- lm(target ~., data = data.temp)
              }else if(lags == 'mse' & engine != 'OLS'){
                model.temp <-
                  grf::regression_forest(
                    X = data.temp %>% select(-target) %>% as.matrix(),
                    Y = data.temp %>% select(target) %>% as.matrix()
                )
              }

              # information criterion
              if(lags == 'aic'){
                accuracy.temp[lags.temp] <- AIC(model.temp)
              }else if(lags == 'bic'){
                accuracy.temp[lags.temp] <- BIC(model.temp)
              }else if(lags == 'mse'){
                accuracy.temp[lags.temp] <-
                  sum(
                    (dplyr::select(data.temp, target)[,1] - predict(model.temp)[,1])^2,
                    na.rm = TRUE) /
                  nrow(data.temp)
              }

            }

            # choose optimal lag
            lags_optimal <- which(accuracy.temp == min(accuracy.temp))
            if(is.vector(lags_optimal)){lags_optimal <- lags_optimal[1]}
            print(paste0('Optimal lags chosen: ', lags_optimal))

            # set final data
            data.temp <- X %>%
              recipes::recipe(~ ., data = X) %>%
              recipes::step_lag(tidyselect::everything(), -tidyselect::contains('states'), lag = 1:lags_optimal) %>%
              recipes::prep(X) %>%
              recipes::bake(X)

            return(data.temp)
          }
        ) %>%
        purrr::map(.f = na.omit) %>%
        as.data.frame()
    }

    #-------------------------------------------------#
    # prepare and generate variable leads          ####
    #-------------------------------------------------#
    # create leads
     for(horizon in 1:horizons){
       data[,paste0('lead_',horizon,'_target')] <- dplyr::lead(data$target, horizon)
    }
    # clear missing
    data = data %>% na.omit()

    #-------------------------------------------------#
    # Create state-dependence                      ####
    #-------------------------------------------------#
    # split by state
    if(!is.null(states)){
      # recast as dataframe [UPDATE]
      data = as.data.frame(data)
      # extract state names
      data <- data %>% dplyr::arrange(states)
      data.states <- paste0('state_',unique(data$states))
      # split into list
      data <- data %>% dplyr::group_split(states)
      # set list names
      names(data) <- data.states
      # remove state form design matrix
      data <- data %>%
        purrr::map(
          function(X){
            return(X %>% select(-states))
          }
        )
    }else{
      data = list(data)
    }

    #-------------------------------------------------#
    # perform local impulse response regressions   ####
    #-------------------------------------------------#
    # create storage matrix
    irfData <- matrix(ncol = 5, nrow = horizons)
    colnames(irfData) <- c('Horizon','Coef','Std.dev','lowerBound','upperBound')
    irfData[,1] <- c(1:horizons)

    # use OLS engine
    if(engine == 'OLS'){
    output <- data %>%
      map(
        function(X){
          # calculate regressions
          for(i in 1:horizons){

            # set experimental data
            Exogenous <- X %>% dplyr::select(-matches('lead_\\d+_target'))
            Y <- X %>% dplyr::select(target.lead = paste0('lead_',i,'_target'))
            designMatrix = data.frame(Exogenous, Y)

            # generate the projections
            directProjection <- lm(target.lead ~., data = designMatrix)

            # Newey-west correction, if desired
            if(NW == TRUE){
              out <- lmtest::coeftest(directProjection,
                                     vcov = sandwich::NeweyWest(directProjection, lags = NW_lags, prewhite = NW_prewhite))
            }else{
              out <- lmtest::coeftest(directProjection)
            }

            # store the data
            shockIndex <- match(shock,rownames(out))
            irfData[i,2] <- out[shockIndex,1]
            irfData[i,3] <- out[shockIndex,2]
            irfData[i,4] <- irfData[i,2] - irfData[i,3]*qnorm(confidence/2)
            irfData[i,5] <- irfData[i,2] + irfData[i,3]*qnorm(confidence/2)
          }
          return(data.frame(irfData))
        }
      )
    # use Tree-based engine
    }else if(engine %in% c('AW','OSW')){
      output <- data %>%
        map(
          function(X){
            for(i in 1:horizons){

                # set experiment data
                Treatment <- X %>% dplyr::select(shock) %>% as.vector()
                High_dim_controls <- X %>% dplyr::select(-shock, -matches('lead_\\d+_target')) %>% as.matrix()
                Exogenous <- X %>% dplyr::select(-shock, -matches('lead_\\d+_target')) %>% as.matrix()
                Y <- X %>% dplyr::select(paste0('lead_',i,'_target')) %>% as.vector()

                # set aside data to calc treatment effect when all other variables are at their historical average
                Simulation <- matrix(rep(0, ncol(Exogenous)), ncol = ncol(Exogenous))
                names(Simulation) <- colnames(Exogenous)

                # create storage for model outputs
                treatment_effects <- vector(length = treeBag)
                treatment_se <- vector(length = treeBag)
                treatment_interval <- matrix(ncol = 2, nrow = treeBag)

                for(j in 1:treeBag){

                  # Use Microsoft's econML engine for instrumental variable tree-based approach
                  if(engine == 'OSW'){

                    # define model
                    est <-
                      econML$ContinuousTreatmentOrthoForest(
                        bootstrap = TRUE,
                        model_T=sklearn$linear_model$LinearRegression(),
                        model_Y=sklearn$linear_model$LinearRegression(),
                        n_jobs = as.integer(cores),
                        random_state = as.integer(seed))

                    # fit model
                    est$fit(Y, Treatment, Exogenous, High_dim_controls, inference = 'blb')

                    # calculate treatment effect and confidence interval
                    treatment_effect = est$effect(np$asarray(Simulation))
                    intervals = est$effect_interval(X = Simulation, alpha = confidence)

                    # store results intervals
                    treatment_effects[j] <- mean(treatment_effect)
                    treatment_se <- NA
                    treatment_interval[j, 1] <- mean(intervals[[1]])
                    treatment_interval[j, 2] <- mean(intervals[[2]])

                  # Use Athey and Wager engine for unconfounded or instrumental variable tree-based approach
                  }else if(engine == 'AW'){

                    # define model and fit model
                    if(IV == FALSE){
                      treeModel <-
                        grf::regression_forest(
                          X = Exogenous,
                          Y = Y %>% as.matrix()
                          )
                    }else if(IV == TRUE){
                      treeModel <-
                        grf::causal_forest(
                          X = Exogenous,
                          Y = Y %>% as.matrix(),
                          W = Treatment %>% as.matrix()
                        )
                    }

                    # calculate treatment effect and confidence interval
                    predictTree <-
                      predict(treeModel,
                              estimate.variance = T,
                              newdata = as.matrix(Simulation))

                    # store results intervals
                    treatment_effects[j] <- predictTree$predictions
                    treatment_se[j] <- predictTree$variance.estimates
                    treatment_interval[j, 1] <- predictTree$predictions - predictTree$variance.estimates
                    treatment_interval[j, 2] <- predictTree$predictions + predictTree$variance.estimates

                  }

                  # store results
                  irfData[i,2] <- mean(treatment_effects)
                  irfData[i,3] <- mean(treatment_se)
                  irfData[i,4] <- mean(treatment_interval[,1])
                  irfData[i,5] <- mean(treatment_interval[,2])
                  irfData <- as.data.frame(irfData)
                }

            }

            # note the initial structural shock, via cholesky decomposition
            if(Cholesky == TRUE){
              A <- shocksMatrix(select(X, data.variables), lags_endog_lin = ncol(select(X, contains('lag'))))
              colnames(A) = rownames(A) = data.variables
              A <- A['target', shock]
              initial <- c('Horizon' = 1,'Coef' = A,'Std.dev' = NA, 'lowerBound' = A, 'upperBound'= A)
              irfData$Horizon <- irfData$Horizon + 1
              irfData <- dplyr::bind_rows(initial, irfData)
            }

            return(irfData)

          }
      )
    }


    #-------------------------------------------------#
    # finalize data and return                     ####
    #-------------------------------------------------#
    # cast single state output as data frame
    if(is.null(states)){output <- data.frame(output)}

    # return proper shock label
    if(shock == 'target'){
      shock <- target
    }

    # package info and cast as an irf object
    information <-
      list(
        impulseResponse = output,
        engine = engine,
        target = target,
        shock = shock,
        lags = lags
      )

    class(information) <- "irf"

    # return output
    return(information)

}

#-------------------------------------------------#
#  Wrapper to produce IRFs of a system ####
#-------------------------------------------------#
#' Estimate local projections
#'
#' [Description]
#'
#'
#' @param data Dataframe of covariate
#' @param shocks Vector of strings denoting variable to shock
#' @param targets Vector of string denoting variable betas to collect
#' @param horizons Horizons to forecast out to
#' @param engine String declaring how to estimate responses (must be 'OLS' 'AW' or 'OSW')
#' @param confidence Alpha used for two-sided confidence intervals
#' @param lags Lags to include in regressions. Must be a deterministic number of lag or an information criterion used to chose the number of lags. Chosen based on information criteria applied to OLS - future functionality may support tree-based ic.
#' @param lags_max Maximum number of lags
#' @param standardize Standardize inputs (mean = 0 variance = 1)
#' @param states Column name used to define states
#' @param NW Newey-West correction on variance-covariance matrix
#' @param NW_lags Number of lags to use in Newey-West correction
#' @param NW_prewhite Prewhite option for Newey-West correction (see sandwich::NeweyWest function)
#' @param IV Perform a two-stage instrumental variable regression (currently only supported for trees)
#' @param treeBag Number of tree-based models to bag over
#' @param Cholesky Append the initial structural shock from a cholesky decomposition [primarly for replication purposes]
#' @param cores Number of cores available to use (must be integer) set -1 to indicate use all available resources
#' @param seed Random seed to be used (must be integer)
#'
#' @return impulseResponse Matrix of impulse responses and standard errors
#' @return engine Method used to estimate impulses
#' @return target Variable that was responding to a shock
#' @return shock Variable that was shocked
#' @return lags Number of lags used in estimation
#'
#' @export

localProjectionIRF_VAR <-
  function(data,                  # dataframe of covariate
           shocks,                 # string denoting variable to shock
           targets,                # string denoting variable betas to collect
           horizons,              # horizons to forecast out to
           engine = 'OLS',        # string declaring how to estimate responses (must be 'OLS', 'AW', or 'OSW')
           confidence = 0.1,      # alpha used for two-sided confidence intervals
           # regressor specifications
           lags = 'bic',          # lags to include in regressions
                                  #   (must be a deterministic number of lag, or an information criterion used to chose the number of lags)
                                  #   (chosen based on information criteria applied to OLS - future functionality may support tree-based ic)
           lags_max = 12,         # maximum number of lags
           standardize = FALSE,   # standardize inputs (mean = 0, variance = 1)
           # Threshold/State dependent parameters
           states = NULL,         # column name used to define states
           # OLS-based IRF parameters
           NW = FALSE,            # Newey-West correction on variance-covariance matrix
           NW_lags = NULL,        # number of lags to use in Newey-West correction
           NW_prewhite = NULL,    # prewhite option for Newey-West correction (see sandwich::NeweyWest function)
           # Tree-based IRF parameteres
           IV = FALSE,            # perform a two-stage instrumental variable regression (currently only supported for trees)
           treeBag = 5,           # number of tree-based models to bag over
           Cholesky = FALSE,       # append the initial structural shock from a cholesky decomposition [primarly for replication purposes]
           # technical options
           cores = NULL,          # number of cores available to use (must be integer), set -1 to indicate use all available resources
           seed  = NULL           # random seed to be used (must be integer)
  ){

  impulses <- list()
  for(s in 1:length(shocks)){
    for(t in 1:length(targets)){

      # identify shock and target variables
      shock <- shocks[s]
      target <- targets[t]

      # keep track of progress
      cat(target,' response to ',shock,'\n')

      # calculate impulse reponses
      impulses[[paste0(target,'.',shock)]] <-
        localProjectionIRF(
          data = data,
          shock = shock,
          target = target,
          horizons = horizons,
          engine = engine,
          confidence = confidence,
          lags = lags,
          lags_max = lags_max,
          standardize = standardize,
          states = states,
          NW = NW,
          NW_lags = NW_lags,
          NW_prewhite = NW_prewhite,
          IV = IV,
          treeBag = treeBag,
          Cholesky = Cholesky,
          cores = cores,
          seed  = seed
        )
    }
  }


  # set object class
  class(impulses) <- "irf_var"

  # return information
  return(impulses)
}


#-------------------------------------------------#
#  Wrapper for lpirf::get_math_chol            ####
#-------------------------------------------------#
shocksMatrix = function(endog_data,
                        lags_endog_lin){
  specs <- list()
  specs$lags_endog_lin <- lags_endog_lin
  specs$lags_criterion <- 'AIC'
  specs$max_lags <- lags_endog_lin
  specs$trend <- 0
  specs$shock_type <- 0
  specs$exog_data <- NULL
  specs$lags_exog <- NULL
  specs$use_twosls <- FALSE
  specs$model_type <- 0
  specs$starts <- 1
  specs$ends <- dim(endog_data)[1]
  specs$column_names <- names(endog_data)
  specs$endog <- ncol(endog_data)
  data_lin <- lpirfs::create_lin_data(specs, endog_data)
  y_lin <- data_lin[[1]]
  x_lin <- data_lin[[2]]
  specs$y_lin <- y_lin
  specs$x_lin <- x_lin
  d <- lpirfs::get_mat_chol(y_lin, x_lin, endog_data, specs)
  return(d)
}

