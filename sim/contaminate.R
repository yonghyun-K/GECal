library(mice)

generate_synthetic_data <- function(confidential_data, missing_fraction = 0.2) {
  
  # Ensure it's a data.frame
  confidential_data <- as.data.frame(confidential_data)
  
  # Introduce NAs randomly but do not replace zeros
  mask <- matrix(runif(nrow(confidential_data) * ncol(confidential_data)) < missing_fraction, 
                 nrow = nrow(confidential_data), ncol = ncol(confidential_data))
  
  # Replace with NA only where values are not zero
  confidential_data_with_na <- confidential_data
  # confidential_data_with_na[mask & confidential_data != 0] <- NA
  confidential_data_with_na[mask] <- NA
  
  # Apply the mice algorithm for imputation
  imputation_result <- mice(confidential_data_with_na, m = 1, method = 'pmm', maxit = 5, printFlag = FALSE)
  
  # Extract the imputed (synthetic) data
  synthetic_data <- complete(imputation_result)
  
  return(synthetic_data)
}

tmpdat <- generate_synthetic_data(Xs1, missing_fraction = 0.25)
tmpdat <- generate_synthetic_data(tmpdat, missing_fraction = 0.25)
tmpdat <- generate_synthetic_data(tmpdat, missing_fraction = 0.25)
tmpdat <- generate_synthetic_data(tmpdat, missing_fraction = 0.25)

sum(tmpdat); sum(Xs1)
colSums(tmpdat); 
colSums(Xs1)

summary(tmpdat)
summary(Xs1)

head(tmpdat * 2)
head(Xs1)

# Save both objects in the same .rda file
save(my_vector, my_dataframe, file = "data/my_objects.rda")

View(Xs)
# Function to contaminate a vector
contaminate_vector <- function(vec, missing_fraction = 0.2, noise_level = 1) {
  # Ensure vec is numeric
  if (!is.numeric(vec)) {
    stop("The input vector must be numeric.")
  }
  
  # Step 2: Add noise to the non-missing values
  noise <- rnorm(length(vec), mean = 0, sd = noise_level * sd(vec, na.rm = TRUE))
  contaminated_vec <- ifelse(is.na(vec), vec, vec + noise)
  
  return(contaminated_vec)
}

tmpvec <- contaminate_vector(y_S)
tmpvec[tmpvec < 0] <- 0

summary(y_S)
summary(tmpvec)



tmpdat <- round(tmpdat)
acresCRD <- tmpdat
total <- total1
pesticideusage <- tmpvec
pimat <- pimat_S
d_S <- d_S0
save(acresCRD, total, pesticideusage, pimat, d_S, file = "data/pesticides.rda")


