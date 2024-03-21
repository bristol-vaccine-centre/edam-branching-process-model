

# Required packages
library(distcrete)
library(epitrix)
library(here)
library(tidyverse) # total overkill, we merely use tibble and dplyr
library(randtoolbox)

# Source files
source("make_disc_gamma.R")
source("branching_model_basic.R")

# Enter the number of observed cases here:
n_cases <- 

# We will accept simulation parameters where simulated case numbers are within
# 20% of the observed case number

percent_tolerance <- 

# Hence, the maximum and minimum simulated number of cases that
# the ABC algorithm will accept are:

n_cases_max <- round(n_cases * (1 + percent_tolerance/100))   
n_cases_min <- round(n_cases * (1 - percent_tolerance/100))  

# Enter serial interval mean and s.d. here:

SI_mean <- 
SI_sd <- 


# Enter other parameters here
max_duration <- 
mean_intro_rate <- 
break_point <- 


# Sample randomly from a grid

min_R_samp <- 
max_R_samp <- 
npoints <- 10000

# Draw R_basic from a uniform distribution from 0 to 3
# (ideally, should use Sobol sequences)

R_sample_location <- runif(n = npoints,
                         min = min_R_samp,
                         max = max_R_samp)

# Also draw the intervention efficacy from a uniform distribution from 0 to 1
points_sampled <- data.frame(R_basic = R_sample_location,
                             intervention_efficacy = runif(n = npoints))

# Number of points we'd like in our posterior distribution
n_posterior <- 

# Set up a dataframe for the posterior distribution
posterior_distribution <- data.frame(R_basic = double(),
                                     intervention_efficacy = double()               
                                    )

# Initiate counters for use in the while loop
posterior_distribution_size <- 0
count <- 1

# We need a while loop for the ABC process: only stop when we reach the
# required number of posterior estimates

while(posterior_distribution_size <= n_posterior) {
  
                                      out <- branching_model_basic(points_sampled[count,1], # This is R_basic
                                                                   points_sampled[count,2], # This is the intervention efficacy
                                                                   serial_interval_mean = SI_mean, 
                                                                   serial_interval_sd = SI_sd,
                                                                   r_daily_intro = mean_intro_rate,
                                                                   max_duration = max_duration,
                                                                   #n_sim = n_sim,
                                                                   break_point = break_point)
      
       # Record the number of cases from this simulation
       n_cases <- nrow(out)

        
       # If the number of cases falls within the acceptable limits,
       # add the pair of values to our posterior distribution
       
       if((n_cases < n_cases_max) & (n_cases > n_cases_min) ){
         
        posterior_sample <- data.frame(R_basic = points_sampled[count,1],
                                       intervention_efficacy = points_sampled[count,2])
         
        posterior_distribution <- rbind(posterior_distribution,
                                        posterior_sample)
          
          # increment the counter for posterior distribution size
        posterior_distribution_size = posterior_distribution_size + 1
        
                                                              }
                      
       # Increment the original counter     
                     count <- count + 1
                                                   }

# Plot the posterior distribution

ggplot(posterior_distribution, 
       aes(x = R_basic,
           y = intervention_efficacy)) +
  
  geom_density_2d_filled() +
  xlab("R for Undetected Cases") +
  ylab("Intervention Efficacy")


