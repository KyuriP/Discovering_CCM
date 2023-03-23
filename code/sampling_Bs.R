
# function creating random B matrix
randomB <- function(B){
  # do I need to explicitly exclude 0?
  samplevals <- runif(sum(B!=0), -0.8, 0.8)
  B[B!=0] <- samplevals  
  return(B)
}


# generate data 
# specify the sample sizes
N <- c(50, 150, 500, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 10000)
# specify replication number
n <- 500
# specify alpha level
alpha <- 0.05


B5dense
ranB <- randomB(B5dense)
# equilibrium check
equilibrium_check(ranB)

# generate data (sample size as specified above)
simdat <- N %>% future_map(function(z) {
  replicate(n = n,
            expr = gen_dat(ranB, N = z),  
            simplify = FALSE)
}, .options = furrr_options(seed=123))