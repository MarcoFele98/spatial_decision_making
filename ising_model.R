library(ggplot2)
library(reshape)

rm(list = ls())

# Functions ____________________________________________________________________________________________________________________________________________________________________________________________________----
heading <- function(coordinates, coordinates_goal) { # calculate heading from target (vectorized)
 # browser()
  mat_direction <- coordinates_goal - coordinates
  mat_direction_squared <- mat_direction * mat_direction
  if(is.matrix(coordinates)) { # for vectorized
    vec_direction <- mat_direction_squared[, 1] + mat_direction_squared[, 2]
    return(mat_direction / sqrt(cbind(vec_direction, vec_direction)))
  }
  else { # for non vectorized
    return(mat_direction/sqrt(sum(mat_direction_squared)))
  }
}

theta <- function(coordinates_1, coordinates_2) { # angle difference between two vectors (vectroized)
  # Vectorized all of this https://stackoverflow.com/questions/1897704/angle-between-two-vectors-in-r
  coordinates_combination <- coordinates_1 * coordinates_2
  coordinates_combination_1 <- coordinates_1 * coordinates_1
  coordinates_combination_2 <- coordinates_2 * coordinates_2
  
  length_vectors <- nrow(coordinates_1)
  if(is.null(length_vectors)) {
    coordinates_combination <- matrix(coordinates_1 * coordinates_2, nrow = 1) 
    coordinates_combination_1 <- matrix(coordinates_1 * coordinates_1, nrow = 1)
    coordinates_combination_2 <- matrix(coordinates_2 * coordinates_2, nrow = 1)
  }
 
  angle_difference_vector <- (coordinates_combination[, 1] + coordinates_combination[, 2]) / (sqrt(coordinates_combination_1[, 1] + coordinates_combination_1[, 2]) * 
                                                                                                sqrt(coordinates_combination_2[, 1] + coordinates_combination_2[, 2]))
  
  angle_difference_vector <- pmin(angle_difference_vector, 1) # to eliminate computing errors
  angle_difference_vector <- pmax(angle_difference_vector, -1) # to eliminate computing errors
  
  return(acos(angle_difference_vector))
}

H <- function(neurons_p, neurons_activation) { # calculate Hamiltonian (energy function). I use a complete network
  #browser()
  n_neurons <- length(neurons_activation)
  
  neurons <- rep(neurons_activation, times = n_neurons) # repeat this to do fast vectorized operations
  neurons_2 <- rep(neurons_activation, each = n_neurons)
  
  neurons_state <- matrix(rep( t(neurons_p), n_neurons), ncol =  2, byrow = TRUE) # repeated the rows of a matrix https://stackoverflow.com/questions/19590541/r-duplicate-a-matrix-several-times-and-then-bind-by-rows-together
  neurons_state_2 <- neurons_p[rep(1:n_neurons, each = n_neurons), ]
  
  j <- cos(pi * (abs(theta(neurons_state, neurons_state_2))/pi)^v)
  return(-k/n * sum(j * neurons_2 * neurons) / 2) # divide by 2 because of the single edge between nodes (instead of double)
}

metropolis_hastings <- function(encoded_headings, neurons_activation, MH_samples) {
  for(i in 1:MH_samples) {
    # find current energy
    old_neurons_activation <- neurons_activation
    current_H <- H(encoded_headings, old_neurons_activation)
    
    # find alternative state
    neuron_to_change <- sample(1:n, 1)
    neurons_activation[neuron_to_change] <- !neurons_activation[neuron_to_change]
    new_H <- H(encoded_headings, neurons_activation)
    
    delta_H <- new_H - current_H
    if(delta_H >= 0) {
      if(rbinom(1, 1, prob = 1 - exp(-(delta_H) / t) )) { # do not change state
        neurons_activation <- old_neurons_activation
      }
    }
  }
  return(neurons_activation)
}

simulate <- function(max_time, MH_samples) {
  #browser()
  # Initial conditions
  position <- c(0, 0)
  neurons_activation <- rbinom(n, 1, 0.5) # start from random neurons activation pattern
  neurons_p <- heading(matrix(position, 
                              byrow = T,
                              ncol = 2, # number of spatial dimensions
                              nrow = n), # number of neurons
                       k_positions[rep(c(1, 2), 
                                       each = n/2), ]) # target 
  # Add noise
  neurons_p[, 1] <- neurons_p[, 1] + rnorm(n, sd = sigma_e)
  neurons_p[, 2] <- neurons_p[, 2] + rnorm(n, sd = sigma_e)
  neurons_p[, 1] <- neurons_p[, 1] / (neurons_p[, 1] * neurons_p[, 1] + neurons_p[, 2] * neurons_p[, 2]) # standardize again
  neurons_p[, 2] <- neurons_p[, 2] / (neurons_p[, 1] * neurons_p[, 1] + neurons_p[, 2] * neurons_p[, 2]) # standardize again
  
  data_positions <- matrix(0, ncol = 2, nrow = max_time) # to keep data
  for(time in 1:max_time) {
    #if(time == 800) browser()
    print(time)
    
    # Find neuronal activation through metropolis hastings optimization
    neurons_activation <- metropolis_hastings(encoded_headings = neurons_p,
                                              neurons_activation = neurons_activation, 
                                              MH_samples = MH_samples) # dynamic scoping of r makes this unnecessary but still more clear in my opinion to pass it as function parameter
   
    # Move
    position <- position + c(mean(neurons_p[, 1] * neurons_activation), 
                             mean(neurons_p[, 2] * neurons_activation)) * v0
    
    # update preferred directions
    neurons_p <- heading(matrix(position, 
                                byrow = T,
                                ncol = 2, # number of spatial dimensions
                                nrow = n), # number of neurons
                         k_positions[rep(c(1, 2), 
                                         each = n/2), ]) # target
    # Add noise
    neurons_p[, 1] <- neurons_p[, 1] + rnorm(n, sd = sigma_e)
    neurons_p[, 2] <- neurons_p[, 2] + rnorm(n, sd = sigma_e)
    neurons_p[, 1] <- neurons_p[, 1] / (neurons_p[, 1] * neurons_p[, 1] + neurons_p[, 2] * neurons_p[, 2]) # standardize again
    neurons_p[, 2] <- neurons_p[, 2] / (neurons_p[, 1] * neurons_p[, 1] + neurons_p[, 2] * neurons_p[, 2]) # standardize again
    
    # save results
    data_positions[time, ] <- position
  }
  
  return(data_positions)
}

# Parameters __________________________________________________________________________________________________________________________________________________________________________________________________----

# Network parameters
n <- 60 # number of neurons
k <- 2 # number of choices
v <- 0.5 # keep it euclidean bro (?)
t <- 0.2 # noise
v0 <- 0.01 # 0.001 # unit of movement
sigma_e <- 0.0002 # 0.02

# Metropolis-Hastings parameters
MH_samples <- 300 

# Set up
set.seed(666)
max_time <- 750
n_replicates <- 10
k_positions <- matrix(c(4.33, 4.33, # x coordinates
                        2.50, -2.50), # y coordinates
                      nrow = 2) # target positions

# Visualize one run of the Metropolis-Hastings algorithm ___________________________________________________________________________________________________________________________________________________________________________----

# example network coordinates for sub-critical: c(1, 0), example for super-critical: c(3, 0)
encoded_headings <- heading(matrix(c(3, # network x coordinate 
                                     0),  # network y coordinate
                            byrow = T,
                            ncol = 2, # number of spatial dimensions
                            nrow = n), # number of neurons
                     k_positions[rep(c(1, 2), 
                                     each = n/2), ]) # target

neurons_activation <- rbinom(60, 1, 0.5)
res <- matrix(neurons_activation, nrow = 1)
for(i in 1:300) {
  print(i)
  # find current energy
  old_neurons_activation <- neurons_activation
  current_H <- H(encoded_headings, old_neurons_activation)
  
  # find alternative state
  neuron_to_change <- sample(1:n, 1)
  neurons_activation[neuron_to_change] <- !neurons_activation[neuron_to_change]
  new_H <- H(encoded_headings, neurons_activation)
  
  delta_H <- new_H - current_H
  if(delta_H >= 0) {
    if(rbinom(1, 1, prob = 1 - exp(-(delta_H) / t) )) { # do not change state
      neurons_activation <- old_neurons_activation
    }
  }
  res <- rbind(res, neurons_activation)
}

rownames(res) <- NULL
test <- melt(res)

ggplot(test) +
  geom_tile(aes(X1, # time
                X2, # neuron
                fill = value)) # on or off


# Simulation __________________________________________________________________________________________________________________________________________________________________________________________________________________-----

data <- cbind(simulate(max_time, MH_samples), 
              1, 
              1:max_time)

for(replicate in 2:n_replicates) {
  print(replicate)
  # Simulate
  data <- rbind(data,
                cbind(simulate(neurons_p, neurons_activation, max_time, MH_samples),
                replicate,
                1:max_time))
}

colnames(data) <- c("x", "y", "replicate", "time")
  
ggplot(as.data.frame(data)) +
  geom_path(aes(x, y, color = replicate, group = replicate)) +
  #geom_point(aes(x, y, color = time)) +
  geom_point(aes(k_positions[1,1], k_positions[1,2]), color = "red", size = 5) +
  geom_point(aes(k_positions[2,1], k_positions[2,2]), color = "red", size = 5)

