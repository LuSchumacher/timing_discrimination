wienerProcess2 <- function(v=1, a=1, z=0.5, ndt=0.3){
  
  ##-- Standard DDM Model Parameters --##
  # v          = drift rate for controlled process
  # a          = boundary separation
  # z          = starting point; = 0.5 corresponds to no bias
  # ndt        = non-decision time in s
  
  # dt         = timesteps
  dt <- 0.001
  sqrt_dt <- sqrt(dt)
  # maxTime    = maximum process duration in s
  maxTime <- 15
  
  # initialize diffusion path for current trial
  path <-  a * z
  noise <- rnorm(maxTime/dt, 0, 1)
  i <- 1
  
  while (path>0 & path<a) {
    path <- path + v*dt + sqrt_dt * noise[i]
    i <- i+1
  }
  
  return(c(as.numeric(path>a),ndt+i*dt))
  
}

