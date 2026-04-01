library("viridis")
cols = viridis(4)
wf = read.table("wf_test.txt")
hist(wf[,1], prob=TRUE, density=20, col=cols[1],
     main="Wright-Fisher transition density started from 0",
     xlab="x", ylab="Density")
hist(wf[,2], prob=TRUE, density=20, col=cols[2], add=TRUE)
hist(wf[,3], prob=TRUE, density=20, col=cols[3], add=TRUE)
hist(wf[,4], prob=TRUE, density=20, col=cols[4], add=TRUE)
lines(c(0,1), c(1,1))
legend("topright",
       legend=c("t = 1.125", "t = 2.25", "t = 4.5", "t = 9", "t = infinite"),
       lwd=3, col=c(cols, "black"))

n = 100000
ode = read.table("ode_test.txt")$V1
for (i in 1:length(ode)) {
  ode[i] = max(ode[i], 0)
}
ode = ode[1:(length(ode) - 1)] - ode[2:length(ode)]
sde = read.table("sde_test.txt")
sde[,4] = sde[,4] / n
sde_proj = rep(0, max(sde[,1]) * 10)
for (i in 1:length(sde[,1])) {
  sde_proj[sde[i,1] * 10] = sde_proj[sde[i,1] * 10] + sde[i,4]
}
# Factors of 1/2 and 2 because ODE discretisation has twice the resolution.
plot((1:length(ode))/2, 2 * ode, type="l", ylim=c(0, max(c(ode, sde_proj))),
     xlim=c(0, max(length(ode) / 2, length(sde))), xlab="Depth (cm)",
     ylab="Stopping power / proton")
par(new=TRUE)
plot(1:length(sde_proj), sde_proj, type="l", ylim=c(0, max(c(ode, sde_proj))),
     xlim=c(0, max(length(ode) / 2, length(sde))), xlab = "", ylab = "", lty=2,
     xaxt="n", yaxt="n")
legend("topleft", legend = c("ODE", "SDE"), lty=1:2, lwd=1)
