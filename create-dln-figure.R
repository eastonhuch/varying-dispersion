#########################
# Figure of DLN Behavior
# Sept 15, 2025
#########################

require(MASS)
require(ggplot2)
library(RColorBrewer)
source("discrete-log-normal.R")

### 
# Heat maps of mean and standard deviation as mu[z] and sigma[z] change

musdz <- expand.grid(muz=seq(-2, 8, length=90), sdz=seq(0.6, 3.23, length=25))

musdz$muy <- integrate_dln(musdz$muz, musdz$sdz, 1)
ey2 <- integrate_dln(musdz$muz, musdz$sdz, 2)
vary <- ey2 - (musdz$muy^2)
musdz$sdy <- sqrt(vary)

#Approximations from the paper
my <- exp(musdz$muz + musdz$sdz/2)
s2y <- exp(musdz$sdz^2 -1)*exp(2*musdz$muz + musdz$sdz^2)
plot(musdz$muy, my)
plot(musdz$sdy^2, s2y)


ggplot(musdz, aes(x=muz, y=sdz, fill=log(muy))) + 
	geom_tile() +
  scale_fill_distiller(palette = "OrRd", direction = 1) + 
  theme_minimal() + 
  labs(
  	x=expression(mu[z]), 
  	y=expression(sigma[z]), 
  	fill=expression("log("~mu[y]~")"), 
  	title="Log Mean of Y for DLN")
  	
ggplot(musdz, aes(x=muz, y=sdz, fill=log(sdy))) + 
	geom_tile() +
  scale_fill_distiller(palette = "OrRd", direction = 1) + 
  theme_minimal() + 
  labs(
  	x=expression(mu[z]), 
  	y=expression(sigma[z]), 
  	fill=expression("log("~sigma[y]~")"), 
  	title="Log Standard Deviation of Y for DLN")

######################
# pmf of y for a few values of mu[z] and sigma[z]

muzlist <- c(0, 1, 2)
y <- 0:10
mypal <- brewer.pal(6, "Dark2")

for(i in 1:3){
	thismuz <- muzlist[i]
	thissigmaz <- 1
	ddln <- pnorm( (log(y+1) - thismuz)/thissigmaz ) - pnorm((log(y) - thismuz)/thissigmaz )
	if(i==1){
		plot(y, ddln, pch=19, col=mypal[1], type='b', lty=2, lwd=1.5, ylab="Density", main=expression("Discrete Log Normal pmf:"~mu[z]))
	}else{
		lines(y, ddln, pch=19, col=mypal[i], type='b', lty=2, lwd=1.5)
	}
}
legend("topright", legend=as.expression(lapply(muzlist, function(m) bquote(mu[z] == .(m)  * "," ~ sigma[z] == 1))), lty=2, lwd=1.5, pch=19, col=mypal[1:3])



sigmazlist <- c(0.5, 1, 2)
y <- 0:10

for(i in 1:3){
	thismuz <- 1
	thissigmaz <- sigmazlist[i]
	ddln <- pnorm( (log(y+1) - thismuz)/thissigmaz ) - pnorm((log(y) - thismuz)/thissigmaz )
	if(i==1){
		plot(y, ddln, pch=19, type='b', lty=2, lwd=1.5, ylab="Density", col=mypal[4], main=expression("Discrete Log Normal pmf:"~sigma[z]))
	}else{
		lines(y, ddln, pch=19, type='b', lty=2, lwd=1.5, col=mypal[i+3])
	}
}
legend("topright", legend=as.expression(lapply(sigmazlist, function(m) bquote(mu[z] == 1 * "," ~ sigma[z] == .(m)))), lty=2, lwd=1.5, pch=19, col=mypal[4:6])



# # all combinations of mu and sigma (too noisy)

# muzlist    <- c(0, 1, 2)
# sigmazlist <- c(1, 2)
# y <- 0:8
# mypal <- paste0(brewer.pal(6, "Paired")[c(2,4,6,1,3,5)])

# params <- expand.grid(mu = muzlist, sigma = sigmazlist)

# plot(NULL, xlim = c(0, 8), ylim = c(0, 0.5),
     # xlab = "y", ylab = "Density")

# for(i in seq_len(nrow(params))){
  # thismu    <- params$mu[i]
  # thissigma <- params$sigma[i]
  
  # ddln <- pnorm( (log(y+1) - thismu)/thissigma ) - 
          # pnorm( (log(y)   - thismu)/thissigma )
  
  # lines(y, ddln, pch=19 - (thissigma-1)*2, col=mypal[i], type='b', lty=2 + (thissigma-1)*2, lwd=1.5)
# }

# # make labels like: μ = 0, σ = 1
# lablist <- as.expression(
  # apply(params, 1, function(row) {
    # bquote(mu == .(row["mu"]) * "," ~ sigma == .(row["sigma"]))
  # })
# )

# legend("topright", legend = lablist, lty=c(2,2,2,4,4,4), lwd=1.5, pch=c(19,19,19,17,17,17), col=mypal[1:nrow(params)])






library(ggplot2)
library(RColorBrewer)
library(patchwork)

# Heatmap 1
p1 <- ggplot(musdz, aes(x=muz, y=sdz, fill=log(muy))) + 
  geom_tile() +
  scale_fill_distiller(palette = "OrRd", direction = 1) + 
  theme_minimal() + 
  labs(
    x = expression(mu[z]), 
    y = expression(sigma[z]), 
    fill = expression("log("~mu[y]~")"), 
    title = "Log Mean of Y for DLN"
  )

# Heatmap 2
p2 <- ggplot(musdz, aes(x=muz, y=sdz, fill=log(sdy))) + 
  geom_tile() +
  scale_fill_distiller(palette = "BuGn", direction = 1) + 
  theme_minimal() + 
  labs(
    x = expression(mu[z]), 
    y = expression(sigma[z]), 
    fill = expression("log("~sigma[y]~")"), 
    title = "Log SD of Y for DLN"
  )

# PMF varying mu
muzlist <- c(2, 2.5, 3)
y <- 0:30
mypal <- brewer.pal(6, "Dark2")

pmf_mu <- data.frame()
sigfix <- 0.2
for(i in seq_along(muzlist)){
  ddln <- pnorm((log(y+1) - muzlist[i])/sigfix) - pnorm((log(y) - muzlist[i])/sigfix)
  pmf_mu <- rbind(pmf_mu, data.frame(y, ddln, mu = muzlist[i], sigma = 1))
}

p3 <- ggplot(pmf_mu, aes(y, ddln, color=factor(mu))) +
  geom_line(linetype=2, linewidth=0.5) +
  geom_point(size=2.5) +
  scale_color_brewer(palette="Dark2") +
  labs(x="y", y="Density",
       color=expression(mu[z]),
       title=expression("DLN pmf (varying"~mu[z]*")")) +
  theme_minimal() + 
  ylim(c(0, 0.32))
p3

# PMF varying sigma
sigmazlist <- c(0.1, 0.2, 0.4)
mufix <- 2.5
pmf_sigma <- data.frame()
for(i in seq_along(sigmazlist)){
  ddln <- pnorm((log(y+1)-mufix)/sigmazlist[i]) - pnorm((log(y)-mufix)/sigmazlist[i])
  pmf_sigma <- rbind(pmf_sigma, data.frame(y, ddln, mu=1, sigma = sigmazlist[i]))
}

p4 <- ggplot(pmf_sigma, aes(y, ddln, color=factor(sigma))) +
  geom_line(linetype=2, linewidth=0.5) +
  geom_point(size=2.5) +
  scale_color_brewer(palette="Set2") +
  labs(x="y", y="Density",
       color=expression(sigma[z]),
       title=expression("DLN pmf (varying"~sigma[z]*")")) +
  theme_minimal() + 
  ylim(c(0,0.32))

# Combine all 4
(p3 | p4) / (p1 | p2)

ggsave("./figures/DLNFigure.png")