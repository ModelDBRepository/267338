library(circular)

phminsYSTD = read.table("phmins_Esm75_YSTD.dat")
phmaxsYSTD = read.table("phmaxs_Esm75_YSTD.dat")

print(length(phminsYSTD$V1))
print(mean(phminsYSTD$V1))
print(sd(phminsYSTD$V1))
print(mean(phmaxsYSTD$V1))
print(sd(phmaxsYSTD$V1))

phminsNSTD = read.table("phmins_Esm75_NSTD.dat")
phmaxsNSTD = read.table("phmaxs_Esm75_NSTD.dat")

print(mean(phminsNSTD$V1))
print(sd(phminsNSTD$V1))
print(mean(phmaxsNSTD$V1))
print(sd(phmaxsNSTD$V1))

setEPS()
postscript("Fig8C2_STD.eps")
rose.diag(phminsYSTD,bins=24,col="blue",ticks=TRUE,border=NA)
rose.diag(phmaxsYSTD,bins=24,col="red",ticks=TRUE,border=NA,add=TRUE)
dev.off()

setEPS()
postscript("Fig8C2_NoSTD.eps")

rose.diag(phminsNSTD,bins=24,col="blue",ticks=TRUE,border=NA)
rose.diag(phmaxsNSTD,bins=24,col="red",ticks=TRUE,border=NA,add=TRUE)
dev.off()

DataMin <- watson.two.test(phminsYSTD, phminsNSTD, alpha=0.001)
DataMax <- watson.two.test(phmaxsYSTD, phmaxsNSTD, alpha=0.0001)
print(DataMin,DataMax)
