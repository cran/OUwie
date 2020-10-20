## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(cache=FALSE)
knitr::opts_chunk$set(tidy.opts=list(width.cutoff=75),tidy=TRUE)

## ---- echo=FALSE, message=FALSE, warning=FALSE---------------------------
require(OUwie)

## ---- eval=TRUE----------------------------------------------------------
data(tworegime)
pp <- OUwie(tree, trait, model="OUM", root.station=TRUE, scaleHeight=TRUE, shift.point=0, algorithm="invert", quiet=TRUE)
pp

## ---- eval=TRUE----------------------------------------------------------
data(tworegime)
pp <- OUwie(tree, trait, model="OUM", root.station=FALSE, scaleHeight=TRUE, shift.point=0, algorithm="invert", quiet=TRUE)
pp

## ---- eval=TRUE, echo=FALSE, fig.height=10, fig.width=6.5, fig.cap = "The edges are painted by regime, assuming an optimum $\\theta_i$ for each color. As shown in Ho and Ane (2014a) the left shows an case of unidentifiability case because every regime forms a connected component. The tree on the right shows a case of identifiability because the black regime covers two completely disconnected parts in the tree."----

par(mfcol=c(1,2), mar=c(4,4,0.5,0.5), oma=c(1,2,1,1))

load("simUnidentifiable.Rsave")
regimes <- OUwie:::GetParameterPainting(phy=phy, data=data, rates=c(3,1,2), k.pars=3)
nb.tip <- Ntip(phy)
nb.node <- Nnode(phy)
comp <- numeric(Nedge(phy))
comp[match(1:(Ntip(phy)), phy$edge[,2])] <- as.factor(regimes[1:nb.tip])
comp[match((2+Ntip(phy)):(Nedge(phy)+1), phy$edge[,2])] <- as.factor(regimes[(nb.tip+1):(nb.tip+nb.node)][-1])
co <- c("black", "blue", "orange")
plot.phylo(phy, edge.color=co[comp], show.tip.label=FALSE)

load("simIdentifiable.Rsave")
regimes <- OUwie:::GetParameterPainting(phy=phy, data=data, rates=c(3,1,2), k.pars=3)
nb.tip <- Ntip(phy)
nb.node <- Nnode(phy)
comp <- numeric(Nedge(phy))
comp[match(1:(Ntip(phy)), phy$edge[,2])] <- as.factor(regimes[1:nb.tip])
comp[match((2+Ntip(phy)):(Nedge(phy)+1), phy$edge[,2])] <- as.factor(regimes[(nb.tip+1):(nb.tip+nb.node)][-1])
co <- c("black", "blue", "orange")
plot.phylo(phy, edge.color=co[comp], show.tip.label=FALSE)


## ---- eval=TRUE----------------------------------------------------------
load("simUnidentifiable.Rsave")
check.identify(phy, data)

## ---- eval=TRUE----------------------------------------------------------
load("simIdentifiable.Rsave")
check.identify(phy, data)

## ---- eval=FALSE---------------------------------------------------------
#  load("simsOUidentify_8")
#  surfaceDatThetaR_2 <- OUwie.contour(oum.root, focal.params=c("theta_Root", "theta_2"), focal.params.upper=c(10,10), nreps=1000, n.cores=4)
#  surfaceDatTheta1_2 <- OUwie.contour(oum.root, focal.params=c("theta_1", "theta_2"), focal.params.upper=c(10,10), nreps=1000, n.cores=4)

## ---- eval=FALSE---------------------------------------------------------
#  plot(surfaceDatThetaR_2, mle.point="red", levels=c(0:20*0.1), xlab=expression(theta[root]), ylab=expression(theta[2]) , xlim=c(0,5,1), ylim=c(0,5,1), col=grey.colors(21, start=0, end=1))

## ---- eval=TRUE, echo=FALSE, fig.height=3.5, fig.width=6.5, fig.cap = "A contour plot of a OUM model, with the $\\theta_{root}$ included in the model. (A) The contour for when $\\theta_{root}$ and $\\theta_{2}$ are specified as the focal parameters, and (B) shows the likelihood surface for when $\\theta_{1}$ and $\\theta_{2}$ are specified. In both cases, the likelihood surface appears as a ridge, indicating that the regimes are not separately identifiable."----

par(mfcol=c(1,2), mar=c(4,4,0.5,0.5), oma=c(1,2,1,1))

load("surfaceDatRIndentR_2.Rsave")
plot(surfaceDatRIndent, mle.point="red", levels=c(0:20*0.1), xlab=expression(theta[root]), ylab=expression(theta[2]) , xlim=c(0,5,1), ylim=c(0,5,1), col=grey.colors(21, start=0, end=1))

load("surfaceDatRIndent1_2.Rsave")
plot(surfaceDatRIndent, mle.point="red", levels=c(0:20*0.1), xlab=expression(theta[1]), ylab=expression(theta[2]) , xlim=c(0,5,1), ylim=c(0,5,1), col=grey.colors(21, start=0, end=1))

## ---- eval=TRUE, echo=FALSE, fig.height=3.5, fig.width=6.5, fig.cap = "A contour plot of a OUM model, with the $\\theta_{root}$ removed from the model. (A) The contour for when $\\theta_{1}$ and $\\theta_{2}$ are specified as the focal parameters, and (B) for when $\\theta_{1}$ and $\\theta_{3}$ are specified, the likelihood surface is sufficiently peaked. In other words, the likelihood surface no longer appears as a ridge and regimes are separately identifiable."----

par(mfcol=c(1,2), mar=c(4,4,0.5,0.5), oma=c(1,2,1,1))

load("surfaceDatNRIndent1_2.Rsave")
plot(surfaceDatNRIndent, mle.point="red", levels=c(0:20*0.1), xlab=expression(theta[1]), ylab=expression(theta[2]) , xlim=c(0,5,1), ylim=c(0,5,1), col=grey.colors(21, start=0, end=1))

load("surfaceDatNRIndent1_3.Rsave")
plot(surfaceDatNRIndent, mle.point="red", levels=c(0:20*0.1), xlab=expression(theta[1]), ylab=expression(theta[3]) , xlim=c(0,5,1), ylim=c(0,5,1), col=grey.colors(21, start=0, end=1))

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
data(tworegime)
set.seed(42)
ouwiefit <- OUwie(tree, trait, model="OUM", scaleHeight=TRUE, root.station=FALSE, shift.point=0.5, algorithm="invert", quiet=TRUE, check.identify=FALSE)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  recon <- OUwie.anc(ouwiefit)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
recon <- OUwie.anc(ouwiefit, knowledge=TRUE)

## ---- eval=TRUE, echo=FALSE, fig.height=10, fig.width=6.5, fig.cap = "A plot of the ancestral state reconstruction under an OUwie model."----
plot(recon, fsize=0.5)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
data(tworegime)
three.point <- OUwie(tree, trait, model="OUMV", root.station=FALSE, scaleHeight=TRUE, shift.point=0.5, algorithm="three.point", quiet=TRUE)
three.point

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
data(tworegime)
invert <- OUwie(tree, trait, model="OUMV", root.station=FALSE, scaleHeight=TRUE, shift.point=0.5, algorithm="invert", quiet=TRUE)
invert

