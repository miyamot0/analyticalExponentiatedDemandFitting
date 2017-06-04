###
#
#   Shawn Gilroy
#   06/04/2017
#
#   GPLv2+
#
#   Proofs for solving exponentiated demand model via analytic gradient
#   
#   Requires: minpack.lm
#
#   Purposes:
#   1) Simulates demand model, with noise for testing with various conditions (enable seed to disable)
#   2) Manually specify out analytic derivatives and hand build jacobian matrix
#   3) Scale individual, low-level errors for e-3,4 parameters
#
###

set.seed(65536)      #  <---- Seed here

require(minpack.lm)

###
#  
#   Demand Function (Exponentiated out form)
#
###
exponentiatedDemandFunc <- function(x, q0, alpha, k) 
{
  func <- expression(q0 * 10^(k * (exp(-alpha * q0 * x) - 1)))
  eval(func)
}

###
#  
#   Jacobian of the demand function (Exponentiated out form)
#
###
exponentiatedGradient <- function(x, q0, alpha, k) 
{
  func <- expression(q0 * 10^(k * (exp(-alpha * q0 * x) - 1)))
  c(eval(D(func, "q0")), 
    eval(D(func, "alpha" )),
    eval(D(func, "k"  )))
}

###
#  
#   Residual helper function
#
###
exponentiatedResidualFunc <- function(params, x, consumption, demandFunction, jacobianFunction)
{
  consumption - do.call("demandFunction", c(list(x = x), as.list(params)))
}

###
#  
#   Jacobian helper function
#
###
exponentiatedJacobianMatrix <- function(params, x, consumption, demandFunction, jacobianFunction) 
{
  -do.call("jacobianFunction", c(list(x = x), as.list(params)))
}  

# Simulated price points
pricePoints <- seq(0.01, 
                   20, 
                   length=20)

# Simulated base parameters
seedParams <- c(q0 = 8, 
                alpha = 0.006, 
                k = 1)

# Fitted base model
demandWithNoise <- do.call("exponentiatedDemandFunc", 
                           c(list(x = pricePoints), 
                           as.list(seedParams)))

# Fitted base model, with noise for simulation
simulatedConsumption <- demandWithNoise + rnorm(length(demandWithNoise), mean = demandWithNoise, sd = 0.2 * max(demandWithNoise))

# nls starting parameters (vague)
startingParameters <- c(q0 = 10, alpha = 0.001, k = 1)

# non-linear LM fitting
#
#   Comments: refer to diag = c(...)
#   This vector applies a scaling constant to the residual function
#   This should even the load for parameters 2-3+ orders away from the others
#   Example: alpha parameter generally 3-4 orders away from highest other parameter
#   
#   Note: diag is not a permanent fix, this is an oddly behaving function and may need tweaking 
#   in certain, more extreme cases
#
fit <- nls.lm(par = startingParameters, 
              fn = exponentiatedResidualFunc, 
              jac = exponentiatedJacobianMatrix,
              demandFunction = exponentiatedDemandFunc, 
              jacobianFunction = exponentiatedGradient,
              x = pricePoints, 
              consumption = simulatedConsumption, 
              control = nls.lm.control(maxiter = 1000,
                                       diag = c(q0 = 1, alpha = 1000, k = 1)))

# Map fitted params to base model
fittedDemand <- do.call("exponentiatedDemandFunc", 
                        c(list(x = pricePoints), 
                        fit$par))  

# Plot out simulated data with fit
plot(pricePoints, 
     simulatedConsumption, 
     cex = 0.5,
     log = "x",
     xlab = "Unit Price",
     ylab = "Consumption",
     main = paste("Analytically solved and scaled LM for Exponentiated Demand\n", 
                  "Q0:", round(fit$par[1], 5),
                  "  alpha:", round(fit$par[2], 5),
                  "  K:", round(fit$par[3], 5),
                  sep = " "))

# draw out fitted solution
lines(pricePoints, fittedDemand, col="red", lwd=1)   
