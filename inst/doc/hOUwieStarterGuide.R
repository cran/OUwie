## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(cache=FALSE)

## ----echo=FALSE, message=FALSE, warning=FALSE---------------------------------
library(OUwie)
library(corHMM)
library(ggplot2)
library(reshape2)

## -----------------------------------------------------------------------------
data(tworegime)
print(head(trait))

## -----------------------------------------------------------------------------
data(tworegime)
dat <- data.frame(sp = tree$tip.label,
                  X = sample(c(0, 1), length(tree$tip.label), replace = TRUE),
                  Y = sample(c(0, 1), length(tree$tip.label), replace = TRUE),
                  FS = rnorm(length(tree$tip.label), 10, 3))
print(head(dat))

## -----------------------------------------------------------------------------
p <- c(0.01670113, 0.39489947, 0.18619839, 1.67259459, 0.16817414) # my fixed set of parameters
pp_oum <- hOUwie(tree, trait, rate.cat = 1, discrete_model = "ER", 
                continuous_model = "OUM", nSim = 25, p = p) # you likely won't use this p argument
print(pp_oum)

## ----message=FALSE, warning=FALSE---------------------------------------------
p_bm1 <- c(0.01663838, 0.13324222, 1.32800027) 
pp_bm1 <- hOUwie(tree, trait, rate.cat = 1, discrete_model = "ER", 
                continuous_model = "BM1", nSim = 25, p = p_bm1)
p_ou1 <- c(0.01686173, 0.16642882, 0.17165229, 1.33223112)
pp_ou1 <- hOUwie(tree, trait, rate.cat = 1, discrete_model = "ER", 
                continuous_model = "OU1", nSim = 25, p = p_ou1) 
p_bmv <- c(0.01662943, 0.06926977, 0.19887583, 1.46411994)
pp_bmv <- hOUwie(tree, trait, rate.cat = 1, discrete_model = "ER", 
                continuous_model = "BMV", nSim = 25, p = p_bmv) 

## -----------------------------------------------------------------------------
model_set <- list(bm1_fit = pp_bm1, ou1_fit = pp_ou1, bmv_fit = pp_bmv, oum_fit = pp_oum)
print(getModelTable(model_set))

## -----------------------------------------------------------------------------
p_bmv_cid <- c(2.915276e+02, 5.319353e-03, 1.728115e-04, 5.529171e+00, 
               2.053390e-01, 1.666330e-01, 1.318453e+00)
pp_bmv_cid <- hOUwie(tree, trait, rate.cat = 2, discrete_model = "ER", null.model = TRUE,
                continuous_model = "BMV", nSim = 25, p = p_bmv_cid)
p_oum_cid <- c(3.206646e+02, 1.639380e-02, 1.172179e-03, 1.639380e+00, 
               1.977996e-01, 2.080072e-01, 3.560215e-01, 1.518203e+00)
pp_oum_cid <- hOUwie(tree, trait, rate.cat = 2, discrete_model = "ER", null.model = TRUE,
                continuous_model = "OUM", nSim = 25, p = p_oum_cid)

## -----------------------------------------------------------------------------
model_set <- list(bm1_fit = pp_bm1, ou1_fit = pp_ou1, bmv_fit = pp_bmv,
                  oum_fit = pp_oum, bmv_cid_fit = pp_bmv_cid, oum_cid_fit = pp_oum_cid)
print(getModelTable(model_set))

## -----------------------------------------------------------------------------
model_avg_pars <- getModelAvgParams(model_set)
print(head(model_avg_pars))

## -----------------------------------------------------------------------------
plot_data <- melt(model_avg_pars)
ggplot(plot_data, aes(x = tip_state, y = value, color = tip_state)) +
  geom_point(size = 5, shape = 21) +
  stat_summary(fun=mean,geom="point",aes(group=1, size = 2)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group=1), width = 0.15, color = "black") +
  theme_classic() +
  facet_wrap(~variable, scales = "free")

## ----echo = FALSE-------------------------------------------------------------
disc_strcuture <- getFullMat(list(getRateCatMat(2),getRateCatMat(2)), getRateCatMat(2))
disc_strcuture[disc_strcuture > 0] <- 1
Q <- disc_strcuture
diag(Q) <- -rowSums(Q)
rownames(Q) <- colnames(Q) <- c("00", "01", "10", "11")
set.seed(1985)
dat <- hOUwie.sim(tree, Q, c(1,0,0,0), c(3,3,3,3), c(1.5,1.5,1.5,1.5), 5, c(5, 5, 10, 10))
flower_data <- data.frame(sp = dat$data$sp, 
                          scent = ifelse(dat$data$reg == 1 | dat$data$reg == 2, 0, 1),
                          color = ifelse(dat$data$reg == 1 | dat$data$reg == 3, "red", "blue"), 
                          size = dat$data$x)

## -----------------------------------------------------------------------------
# simulated data
print(head(flower_data))

## -----------------------------------------------------------------------------
plot(tree, show.tip.label = FALSE, x.lim = c(0, 4), no.margin = TRUE)
tiplabels(pch = 21, col = flower_data$scent, cex = 0.75, offset = 0.02)
tiplabels(pch = 16, col = flower_data$color, cex = 0.75, offset = 0.05)
proportions <- dat$data$x/max(dat$data$x)
plotting_matrix <- matrix(c(rep(3.6, length(tree$tip.label)), 
                    1:length(tree$tip.label), 3.6 + (0.4 * proportions), 
                    1:length(tree$tip.label)), ncol = 4)
for(i in 1:64){
  lines(t(matrix(plotting_matrix[i,], 2, 2)), lwd = 2, col = "darkgrey")
}

## -----------------------------------------------------------------------------
ou1_structure <- getOUParamStructure(model = "OU1", nObsState = 4, rate.cat = 1, 
                                     null.model = FALSE)
print(ou1_structure)

## -----------------------------------------------------------------------------
p_ou1 <- c(1, 3, 1.5, 7.5) # note: not the mle
ou1_fit <- hOUwie(tree, flower_data, rate.cat = 1, discrete_model = "ER", 
                continuous_model = ou1_structure, nSim = 50, p = p_ou1)
print(ou1_fit)

## -----------------------------------------------------------------------------
print(getStateMat4Dat(flower_data[,-4])$legend)

## -----------------------------------------------------------------------------
m1_structure <- ou1_structure
print("Original OU1 structure")
print(m1_structure)
m1_structure[3,] <- c(3, 3, 4, 4)
print("Modified to allow differing optima based on flower scent")
print(m1_structure)

## -----------------------------------------------------------------------------
m3_structure <- m2_structure <- ou1_structure
print("Original OU1 structure")
print(m1_structure)
m2_structure[3,] <- c(3, 4, 3, 4)
print("M2: Depends on flower color")
print(m2_structure)
m3_structure[3,] <- c(3, 4, 5, 6)
print("M3: Depends on flower scent and flower color")
print(m3_structure)


## -----------------------------------------------------------------------------
p_m1 <- c(1, 1, 1.5, 5, 10) # note: not the mle
m1_fit <- hOUwie(tree, flower_data, rate.cat = 1, discrete_model = "ER", 
                continuous_model = m1_structure, nSim = 50, p = p_m1)
p_m2 <- c(1, 3, 1.5, 5, 10) # note: not the mle
m2_fit <- hOUwie(tree, flower_data, rate.cat = 1, discrete_model = "ER", 
                continuous_model = m2_structure, nSim = 50, p = p_m2)
p_m3 <- c(1, 3, 2, 5, 11, 6, 12) # note: not the mle
m3_fit <- hOUwie(tree, flower_data, rate.cat = 1, discrete_model = "ER", 
                continuous_model = m3_structure, nSim = 50, p = p_m3)

## -----------------------------------------------------------------------------
print(getModelTable(list(ou1 = ou1_fit, m1 = m1_fit, m2 = m2_fit, m3 = m3_fit)))

