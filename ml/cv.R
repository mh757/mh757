############################################################
## Code for HW4 of MATH 574M - Statistical ML, Fall 2024.
##
## Perform CV for model selection of k-NN.
## Perform CV for estimating test error for LDA and GLM.
##
## Author: Marcel Hudiani
############################################################


library(caret)
library(class)
library(MASS)
library(ggplot2)

####################
## Simulation Setup
####################

# level of verbosity = 0, 1, or 2.
VERBOSITY <- 1

# Enablers to run specific parts of the code.
RUN_LDA <- TRUE
RUN_GLM <- TRUE
RUN_KNN <- TRUE

################################################
## Compute prediction error for a classifier
##
## Input:
## 1. Predict output.
## 2. Expected output data (as factor).
################################################
get_classifier_error <- function(predictions, ytest_as_factor) {
  
  USE_CONFUSION_MATRIX <- FALSE
  prediction_error <- 0
  
  if (USE_CONFUSION_MATRIX) {
    # Create a confusion matrix
    confusion <- confusionMatrix(predictions$class,
                                 reference = ytest_as_factor)
    
    # View the confusion matrix
    if (VERBOSITY > 1) {
      print(confusion)
    }
    
    # Calculate accuracy
    accuracy <- confusion$overall['Accuracy']
    print(paste("Accuracy:", round(accuracy, 4)))
    
    # Compute prediction error (1 - accuracy)
    prediction_error <- 1 - accuracy
  } else {
    #print("Use conventional method to compute test error.")
    prediction_error <- mean(predictions$class != ytest_as_factor)
  }
  
  if (VERBOSITY > 0) {
    print(paste("Prediction Error:", round(prediction_error, 4)))
  }
  
  return(prediction_error)
}

#####################################################
# Wrapper for k-NN prediction.
# WARNING: binary classification only.
#
# Input:
# 1. Fitted k-NN model
# 2. Test data
#####################################################
test_knn <- function(fittedModel, xytest) {
  
  # Input/output setup
  ytest <- subset(xytest, select = y)
  
  # Initialize
  predictions <- data.frame(class = fittedModel, posterior = fittedModel)
  
  return(get_classifier_error(predictions, ytest$y))
}

#####################################################
# Wrapper for LDA prediction.
# WARNING: binary classification only.
#
# Input:
# 1. Fitted LDA model
# 2. Test data
#####################################################
test_lda <- function(fittedModel, xytest) {
  
  # Input/output setup
  xtest <- subset(xytest, select = -y)
  ytest <- subset(xytest, select = y)
  
  # Predict using the model
  predictions <- predict(fittedModel, newdata = xtest)
  
  return(get_classifier_error(predictions, ytest$y))
}

#####################################################
# Wrapper for GLM prediction.
# WARNING: binary classification only.
#
# Input:
# 1. Fitted GLM model
# 2. Test data
# 3. Labels for output data.
#####################################################
test_glm <- function(fittedModel, xytest, ylabels) {
  
  # Input/output setup
  xtest <- subset(xytest, select = -y)
  ytest <- subset(xytest, select = y)
  
  # Predicted values (GLM outputs probabilities)
  probs <- predict(fittedModel, newdata = xtest)
  dlen <- length(probs)
  
  # Initialize
  predictions <- data.frame(class = rep(0, dlen), posterior = probs)
  
  # Wrap
  for (i in 1:dlen) {
    if (probs[i] > 0) {
      predictions$class[i] <- 1
    }
  }
  
  # Formatting
  predictions$class <- factor(predictions$class, labels = ylabels)
  
  if (is.factor(ytest$y)) {
    return(get_classifier_error(predictions, ytest$y))
  } else {
    stop("ytest is not a factor type.")
  }
}

#####################################################
# Wrapper for fitting LDA to be used in a CV.
# CV passes both training and validation data.
# Therefore, one must introduce a wrapper since LDA
# does not require validation data for fitting.
#####################################################
lda_wrapper <- function(formula, data, test, prior) {
  
  # Remove 'test' or validation data
  args <- list(formula = formula,
               data = data,
               prior = prior)
  
  fit_lda <- do.call("lda", args)
  
  return(fit_lda)
}

#####################################################
# Wrapper for fitting GLM to be used in a CV.
# CV passes both training and validation data.
# Therefore, one must introduce a wrapper since LDA
# does not require validation data for fitting.
#####################################################
glm_wrapper <- function(formula, data, test, family) {
  
  # Remove 'test' or validation data
  args <- list(formula = formula,
               data = data,
               family = family)
  
  fit_glm <- do.call("glm", args)
  
  return(fit_glm)
}


#####################################################
# Wrapper for fitting KNN to be used in a CV.
#
#####################################################
knn_wrapper <- function(data, test, cl, k) {
  
  #'data' is the training data without y-value.
  xtrain <- subset(data, select = -y)
  
  #'test' is the testing data.
  xtest <- subset(test, select = -y)
  
  #'cl' is the classlist corresponding to the training data
  cl <- subset(data, select = y)
  
  if (VERBOSITY > 1) {
    print(paste("knn_wrapper: k-value: ", k))
  }
  
  fit_knn <- knn(xtrain, xtest, cl$y, k)
  
  return(fit_knn)
}


#####################################################
# Performs KNN model selection by using CV
#
#####################################################
knn_model_selection <- function(data_cv, k_cv, k_values,
                                fit_args, test_args) {
  
  predError <- rep(0, max(k_values))
  
  for (k in k_values) {
    # Update k-value
    fit_args$k <- k
    
    # Run the CV
    predError[k] <- runCV(data_cv,
                          "knn_wrapper", fit_args,
                          "test_knn", test_args,
                          kfold = k_cv,
                          seeds = c(1))
  }
  
  # FIX: Need nested CV for unbiased CV error.
  return(predError)
}

#####################################################
# Outer wrapper for CV
#
# Input:
# 1. CV data (D),
# 2. Fitting command/arguments,
# 3. Test command/arguments,
# 4. Number of folds (kfold),
# 5. List of seed numbers (seeds).
#####################################################
runCV <- function(D, fit_cmd, fit_args, test_cmd, test_args,
                  kfold = 5, seeds = c(1)) {
  
  # Initialize
  cvErr <- rep(0, length(seeds))
  
  # Outer wrapper to estimate CV error
  for (i in 1:length(seeds)) {
    set.seed(seeds[i])
    
    # Shuffle data
    D_randomized <- D[sample(nrow(D)),]
    
    # Build error matrix
    errVec <- rep(0, kfold)
    
    # Inner CV
    for (k in 1:kfold) {
      dlen <- nrow(D)
      
      # Upper and lower test (left-out) index
      ltestidx <- 1 + floor((k - 1) / kfold * dlen)
      utestidx <- floor(k/kfold * dlen)
      
      # Assign training and validation data
      tr_data <- D_randomized[-(ltestidx : utestidx),]
      val_data <- D_randomized[(ltestidx : utestidx),]
      
      # Fit model
      fit_args$data <- tr_data
      fit_args$test <- val_data
      fit_model <- do.call(fit_cmd, fit_args)

      # Get the error
      test_args$fittedModel <- fit_model
      test_args$xytest <- val_data
      errVec[k] <- do.call(test_cmd, test_args)
    }
    
    # Calculate and report CV error
    cvErr[i] <- mean(errVec)
    
    if (VERBOSITY > 0) {
      print(paste("Model: ", fit_cmd,
                  "Seed: ", seeds[i],
                  "CV error: ", cvErr[i]))
    }
  }
  
  # Return CV error data
  return(cvErr)
}

###############################
## MAIN
###############################
setwd("C:/Users/mhudi/doc/b/[STAT]/stat_574M_f24/hw4")

s1_traindata <- read.table("../hw3/binary_class_scenario_1_tr_data.csv")
s2_traindata <- read.table("../hw3/binary_class_scenario_2_tr_data.csv")

s1_testdata <- read.table("../hw3/binary_class_scenario_1_tst_data.csv")
s2_testdata <- read.table("../hw3/binary_class_scenario_2_tst_data.csv")


#######################################
### Pre-conditioning
#######################################
cv_kfold <- 5
cv_data <- subset(s1_traindata, select = c(y,x1,x2))

## Format CV data into a factor
##
## red = 0  : disable, no-go
## green = 1: enable, go
##
lvs <- c("red", "green")
cv_data$y <- factor(cv_data$y, labels = lvs)

#######################################
### LDA
#######################################
if (RUN_LDA) {
  
  ## problem 6.b.
  print("LDA: Running CV to predict the test error...")
  
  fit_args <- list(formula = y ~ .,
                 data = cv_data,
                 test = NULL,
                 prior = c(1,1)/2)
  
  test_args <- list(fittedModel = NULL, xytest = cv_data)
  
  cvErrData <- runCV(cv_data,
                   "lda_wrapper", fit_args,
                   "test_lda", test_args,
                   kfold = cv_kfold,
                   seeds = c(1,15,16,47,199,432,243,174,341,943,
                             2132,1558,6000,5234,8762,2356,8479,9847,
                             56075,12523))

  ## problem 6.c.
  hist(cvErrData,
      main = "Histogram of CV Error Data over 20 seeds",
      xlab = "CV Error",
      ylab = "Frequency",
      col = "lightblue",
      border = "black")

  ## problem 6.d.
  print("LDA: Computing the test error directly...")
  
  test_data <- subset(s1_testdata, select = c(y,x1,x2))
  test_data$y <- factor(test_data$y, labels = lvs)
  
  fit_lda <- lda(y ~ ., cv_data, fit_args$prior)
  test_error <- test_lda(fit_lda, xytest = test_data)
}

#######################################
### GLM
#######################################
if (RUN_GLM) {
  
  ## problem 7.b.
  print("GLM: Running CV to predict the test error...")
  
  fit_args <- list(formula = y ~ .,
                   data = cv_data,
                   test = NULL,
                   family=binomial(link="logit"))
  
  test_args <- list(fittedModel = NULL, xytest = cv_data, ylabels = lvs)
  
  runCV(cv_data,
        "glm_wrapper", fit_args,
        "test_glm", test_args,
        kfold = cv_kfold,
        seeds = c(1))
}

#######################################
### KNN
#######################################
if (RUN_KNN) {
  
  print("k-NN: Running CV for model selection...")
  
  # Specify fit arguments for KNN
  fit_args <- list(data = NULL,
                   test = NULL,
                   cl = NULL,
                   k=7)
  
  # Specify test arguments for KNN
  test_args <- list(fittedModel = NULL, xytest = cv_data)
  
  # Specify KNN k-values
  k_vals = 1:100
  knn_test_data <- subset(s2_testdata, select = c(y,x1,x2))
  knn_test_data$y <- factor(knn_test_data$y, labels = lvs)
  
  # Perform model selection
  knn_ms_err <- knn_model_selection(cv_data, cv_kfold,
                                  k_vals, fit_args, test_args)
  
  ms_k_opt <- which.min(knn_ms_err)
  
  # Report on best "k"
  print(paste("k-NN: Best k (from CV):", ms_k_opt))
  
  # Perform regular testing for each k
  print("k-NN: Computing test error...")
  
  testError <- rep(0, max(k_vals))
  
  for (k in k_vals) {
    fit_knn <- knn_wrapper(cv_data, knn_test_data, cv_data$y, k)
    testError[k] <- test_knn(fit_knn, knn_test_data)
  }
  
  # Plot CV error and test error
  results <- data.frame(k = k_vals,
                        testError = testError,
                        cvError = knn_ms_err)
  
  ggplot(results, aes(x = k)) +
    geom_line(aes(y = testError, color = "Test Error"), linewidth = 1) +
    geom_line(aes(y = cvError, color = "CV Error"), linewidth = 1) +
    geom_point(aes(y = testError, color = "Test Error")) +
    geom_point(aes(y = cvError, color = "CV Error")) +
    labs(title = "Test Error vs. k in k-NN",
         x = "Number of Neighbors (k)",
         y = "Error",
         color = "legend") +
    theme_minimal()
}