### TSOL MOZAMBIQUE
### -- helper functions

##
## FIT GAMMA
##

## fit a gamma distribution to a 95% confidence interval
fit_gamma <-
function(lwr, upr) {
  f <-
    function(par, target) {
      ## calculate sum of squared deviations
      sum((qgamma(c(.025, .975), par[1], par[2]) - target) ^ 2)
    }
    
  optim(c(1, 1), f, target = c(lwr, upr))
}


##
## RESIDUAL LIFE EXPECTANCY
##

## define ages for which RLE is available
age <- c(0, 1, 5 * 1:19)

## GBD2010 standard life expectancy table
le_gbd <-
  c(86.02, 85.21, 81.25, 76.27, 71.29, 66.35, 61.40,
    56.46, 51.53, 46.64, 41.80, 37.05, 32.38, 27.81,
    23.29, 18.93, 14.80, 10.99,  7.64,  5.05,  3.31)

## WHO/GHE standard life expectancy table
le_who <-
  c(91.94, 91.00, 87.02, 82.03, 77.04, 72.06, 67.08,
    62.11, 57.15, 52.20, 47.27, 42.36, 37.49, 32.65,
    27.86, 23.15, 18.62, 14.41, 10.70,  7.60,  5.13)

## function to interpolate life expectancy table369
rle <-
function(x, le = c("gbd", "who")) {
  le <- match.arg(le)
  LE <- switch(le,
               gbd = le_gbd,
               who = le_who)
  approx(age, LE, x)$y
}


##
## BURDEN
##

## function to calculate integral
f <-
function(x, K, C = .1658, beta = .04, r, a) {
  K * C * x * exp(-beta * x) * exp(-r * (x - a)) +
    (1 - K) * exp(-r * (x - a))
}

## burden calculation function
burden <-
function(N, DW, A, L, K, r, a) {
  N * DW * integrate(f, lower = A, upper = A + L, K = K, r = r, a = a)$value
}


###
### SENSITIVITY ANALYSIS FUNCTIONS
###

## partial correlation coefficients
sa_pcc <-
function(y, x) {
  out <- matrix(ncol = 2, nrow = ncol(x))
  colnames(out) <- c("rho", "p")
  rownames(out) <- colnames(x)

  for (i in seq(ncol(x))){
    lm_y <- lm(y ~ x[, -i])      # regress y to other x's
    lm_x <- lm(x[, i] ~ x[, -i]) # regress x to other x's
    out[i, ] <-
      unlist(cor.test(lm_y$residuals, lm_x$residuals)[4:3],
             use.names = FALSE)
  }

  return(out[order(abs(out[, "rho"]), decreasing = TRUE), ])
}

## ggplot2 tornado graph
tornado <-
function(coef, names = NULL) {
  ## copy names if needed
  if (is.null(names)) names <- rownames(coef)

  ## create data frame
  df <- data.frame(est = coef[, "rho"],
                   order = order(abs(coef[, "rho"])),
                   name = names)

  ## sort data frame
  df <- df[df$order, ]

  ## create ggplot
  ggplot(df, aes(x = order, y = est)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_x_continuous(element_blank(),
                       breaks = seq(nrow(df)),
                       labels = df$name) +
    scale_y_continuous("partial correlation coefficient",
                       limits = c(min(0, min(df$est) - 0.1),
                                  max(0, max(df$est) + 0.1))) +
    geom_text(aes(x = order, y = est, label = formatC(est, 3, form = "f")),
              size = 3,
              hjust = ifelse(df$est > 0, -0.1, 1.1),
              vjust = 0.4) +
    theme_bw()
}
