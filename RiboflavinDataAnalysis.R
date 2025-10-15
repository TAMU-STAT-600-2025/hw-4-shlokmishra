# Load the riboflavin data

# Uncomment below to install hdi package if you don't have it already; 
# install.packages("hdi") 
library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene expression
dim(riboflavin$x) # n = 71 samples by p = 4088 predictors
?riboflavin # this gives you more information on the dataset

# This is to make sure riboflavin$x can be converted and treated as matrix for faster computations
class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]


# Get matrix X and response vector Y
X = as.matrix(riboflavin$x)
Y = riboflavin$y

# Source your lasso functions
source("LassoFunctions.R")

# [ToDo] Use your fitLASSO function on the riboflavin data with 60 tuning parameters
cat("Fitting LASSO on riboflavin data with 60 tuning parameters...\n")
riboflavin_fit <- fitLASSO(X, Y, n_lambda = 60)
cat("LASSO fitting completed!\n")
cat("Data dimensions: n =", nrow(X), "samples, p =", ncol(X), "predictors\n")
cat("Lambda range:", min(riboflavin_fit$lambda_seq), "to", max(riboflavin_fit$lambda_seq), "\n")

# [ToDo] Based on the above output, plot the number of non-zero elements in each beta versus the value of tuning parameter
# Calculate sparsity (number of non-zero coefficients) for each lambda
sparsity <- colSums(riboflavin_fit$beta_mat != 0)

# Create sparsity plot
par(mfrow = c(1, 2))
plot(riboflavin_fit$lambda_seq, sparsity, 
     type = "l", log = "x",
     xlab = "Lambda", ylab = "Number of non-zero coefficients",
     main = "Sparsity vs Lambda",
     col = "blue", lwd = 2)
grid()

# Also plot the coefficient paths
matplot(riboflavin_fit$lambda_seq, t(riboflavin_fit$beta_mat), 
        type = "l", log = "x",
        xlab = "Lambda", ylab = "Coefficient values",
        main = "LASSO Coefficient Paths",
        col = "gray", lty = 1)
grid()

# [ToDo] Use microbenchmark 10 times to check the timing of your fitLASSO function above with 60 tuning parameters
library(microbenchmark)

cat("\nTiming fitLASSO with microbenchmark (10 runs)...\n")
timing_results <- microbenchmark(
  fitLASSO(X, Y, n_lambda = 60),
  times = 10
)

print(timing_results)
median_time <- median(timing_results$time) / 1e9  # Convert to seconds
cat("Median timing:", round(median_time, 2), "seconds\n")


# [ToDo] Report your median timing in the comments here: (~5.8 sec for Irina on her laptop)
# Median timing: 2.24 seconds (on my machine) 

# [ToDo] Use cvLASSO function on the riboflavin data with 30 tuning parameters (just 30 to make it faster)
cat("\nPerforming 5-fold cross-validation with 30 tuning parameters...\n")
riboflavin_cv <- cvLASSO(X, Y, n_lambda = 30, k = 5)
cat("Cross-validation completed!\n")
cat("Selected lambda_min:", riboflavin_cv$lambda_min, "\n")
cat("Selected lambda_1se:", riboflavin_cv$lambda_1se, "\n")

# [ToDo] Based on the above output, plot the value of CV(lambda) versus tuning parameter. Note that this will change with each run since the folds are random, this is ok.
# Create CV plot
par(mfrow = c(1, 1))
plot(riboflavin_cv$lambda_seq, riboflavin_cv$cvm, 
     type = "l", log = "x",
     xlab = "Lambda", ylab = "CV Error",
     main = "Cross-Validation Error vs Lambda",
     col = "red", lwd = 2)

# Add error bars
arrows(riboflavin_cv$lambda_seq, 
       riboflavin_cv$cvm - riboflavin_cv$cvse,
       riboflavin_cv$lambda_seq, 
       riboflavin_cv$cvm + riboflavin_cv$cvse,
       length = 0.05, angle = 90, code = 3, col = "red")

# Mark selected lambdas
abline(v = riboflavin_cv$lambda_min, col = "blue", lty = 2, lwd = 2)
abline(v = riboflavin_cv$lambda_1se, col = "green", lty = 2, lwd = 2)

# Add legend
legend("topright", 
       legend = c("CV Error", "lambda_min", "lambda_1se"),
       col = c("red", "blue", "green"),
       lty = c(1, 2, 2), lwd = c(2, 2, 2))

grid()

# Show sparsity at selected lambdas
lambda_min_idx <- which.min(abs(riboflavin_cv$lambda_seq - riboflavin_cv$lambda_min))
lambda_1se_idx <- which.min(abs(riboflavin_cv$lambda_seq - riboflavin_cv$lambda_1se))

cat("\nSparsity analysis:\n")
cat("At lambda_min:", sum(riboflavin_cv$beta_mat[, lambda_min_idx] != 0), "non-zero coefficients\n")
cat("At lambda_1se:", sum(riboflavin_cv$beta_mat[, lambda_1se_idx] != 0), "non-zero coefficients\n")

# Show some statistics
cat("\nSummary statistics:\n")
cat("CV error range:", min(riboflavin_cv$cvm), "to", max(riboflavin_cv$cvm), "\n")
cat("CV error at lambda_min:", riboflavin_cv$cvm[lambda_min_idx], "\n")
cat("CV error at lambda_1se:", riboflavin_cv$cvm[lambda_1se_idx], "\n")

