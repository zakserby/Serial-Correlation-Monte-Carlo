#===============================================================================
# AEM 6850
# Naive vs Newey-West standard errors under serial correlation
#===============================================================================

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 1). Preliminary -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Clean up workspace and load or install necessary packages if necessary
rm(list=ls())
want <- c("sandwich")
need <- want[!(want %in% installed.packages()[,"Package"])]
if (length(need)) install.packages(need)
lapply(want, function(i) require(i, character.only=TRUE))
rm(want, need)

# Working directories
dir <- list()
dir$root <- dirname(getwd())
dir$output_figure <- paste(dir$root,"/output_figure",sep="")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 2). Main code -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Generate data
n_vec <- c(25,50,100,200,400,800,1600)
beta <- 1
rho_vec <- c(0,0.25,0.5,0.75)  # autocorrelation coefficient

# Settings
r <- 1000
alpha <- 0.05

# Create a grid to input the rejection rate results
results <- expand.grid(n=n_vec, rho=rho_vec)
results$reject_naive <- NA
results$reject_nw <- NA

# Loop over sample sizes and autocorrelation coefficients
set.seed(123)
for (i in 1:nrow(results)) {
  
  n <- results$n[i] # Using a single value for every iteration
  rho <- results$rho[i]
  
  # Vectors to store rejections for each simulation
  reject_naive_vec <- numeric(r)
  reject_nw_vec <- numeric(r)
  
  for (sim in 1:r) {
    
    # Generate serial correlated X
    x <- numeric(n)
    x[1] <- rnorm(1)
    for (t in 2:n) x[t] <- rho * x[t-1] + rnorm(1)
    
    # Error term
    e <- numeric(n)
    e[1] <- rnorm(1)
    for (t in 2:n) e[t] <- rho * e[t-1] + rnorm(1)
    
    # DGP
    y <- beta * x + e
    
    # Estimate OLS
    ols <- lm(y ~ x)
    
    # SEs
    naive_se <- sqrt(diag(vcov(ols)))
    nw_cov <- NeweyWest(ols, lag=4)
    nw_se <- sqrt(diag(nw_cov))
    
    # t-stats
    naive_t <- (coef(ols)[2] - beta)/naive_se[2]
    nw_t <- (coef(ols)[2] - beta)/nw_se[2]
    
    # 2-sided rejection
    crit <- qt(1-alpha/2, df=n-2)
    reject_naive_vec[sim] <- abs(naive_t) > crit
    reject_nw_vec[sim] <- abs(nw_t) > crit
  }
  
  # Store average rejection rates
  results$reject_naive[i] <- mean(reject_naive_vec)
  results$reject_nw[i] <- mean(reject_nw_vec)
}


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# 3). Plot -----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Save as png
png(paste0(dir$output_figure, "/Rejection_rates_Serial_Dependence.png"), width=1400, height=850, res = 150)

par(mfrow=c(1,2), mar=c(5,5,4,2))

rho_colors_naive <- c("#f6a5a5", "#f58c8c", "#f26b6b", "#ef3b3b")
rho_colors_nw    <- c("#aecbeb", "#89b4e5", "#5f9fe0", "#2f7fcf")
rho_pch <- c(0, 1, 2, 3)

# Plotting Naive SEs
# Generate plot
plot(NULL, xlim = range(n_vec), ylim = c(0, 0.40), xlab = "Sample size", ylab = "Rejection rate",
     log  = "x", cex.lab = 1.25, cex.main = 1.75, axes = F)

# Axis labels
axis(2, at = c(0, 0.05, 0.10, 0.20, 0.30, 0.40),
     labels = c("0", "0.05", "0.10", "0.20", "0.30", "0.40"), las = 1) # y-axis
axis(1, at = n_vec, labels = n_vec) # x-axis

box()
abline(h = 0.05, lty = 2, lwd = 1.35)

# Plot Naive SEs
for (i in 1:length(rho_vec)) {
  indiv_lines <- subset(results, rho == rho_vec[i])
  lines(indiv_lines$n, indiv_lines$reject_naive, col = rho_colors_naive[i], lwd = 1.9)
  points(indiv_lines$n, indiv_lines$reject_naive, pch = rho_pch[i], col = rho_colors_naive[i], cex = 1.35)
}

legend("top", title = expression(rho), legend = rho_vec, col = rho_colors_naive, pch = rho_pch, 
       horiz = TRUE, cex = 0.85, pt.cex = 0.75, lty = 1, lwd = 1.15, inset = c(0, 0.02))

mtext(text = "Naive SEs", side = 3, line = 1.2, font = 2, col = "#b32020", cex = 1.45)


# Plotting Newey-West SEs
# Generate plot
plot(NULL, xlim = range(n_vec), ylim = c(0, 0.40), xlab = "Sample size", ylab = "Rejection rate",
     log  = "x", cex.lab = 1.25, cex.main = 1.75, axes = F)

# Axis labels
axis(2, at = c(0, 0.05, 0.10, 0.20, 0.30, 0.40),
     labels = c("0", "0.05", "0.10", "0.20", "0.30", "0.40"), las = 1) # y-axis
axis(1, at = n_vec, labels = n_vec) # x-axis

box()
abline(h = 0.05, lty = 2, lwd = 1.35)

# Plot Newey-West SEs
for (i in 1:length(rho_vec)) {
  indiv_lines <- subset(results, rho == rho_vec[i])
  lines(indiv_lines$n, indiv_lines$reject_nw, col = rho_colors_nw[i], lwd = 1.9)
  points(indiv_lines$n, indiv_lines$reject_nw, pch = rho_pch[i], col = rho_colors_nw[i], cex = 1.35)
}

legend("top", title = expression(rho), legend = rho_vec, col = rho_colors_nw, pch = rho_pch, 
       horiz = TRUE, cex = 0.85, pt.cex = 0.75, lty = 1, lwd = 1.15, inset = c(0, 0.02))

mtext(text = "Newey-West SEs", side = 3, line = 1.2, font = 2, col = "#0f4c81", cex = 1.45)

dev.off()

