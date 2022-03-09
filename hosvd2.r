#HOSVD of the matrix A

#Construction of the data tensor

tenseur <- array(dim=c(511, 7, 9000)) #g,e,t

for(i in 0:8999){
	wavefunction <- read.table(paste0("wavefunction/vector", i, ".txt"))
	wavefunction <- complex(real=wavefunction$V1, imag=wavefunction$V2)
	
	mat.a <- matrix(wavefunction, nrow=511, ncol=7)
	
	tenseur[ , , (i + 1)] <- mat.a
	
#	print(i)
}

#Unfolding of the tensor

unfolding1 <- matrix(tenseur, nrow=511) #g, t * e_max + e

unfolding1 <- matrix(aperm(tenseur, perm=c(1,2,3)), nrow=511) #g, t * e_max + e
#aperm is useless here, but it allows to understand the pattern for the other unfoldings

unfolding2 <- matrix(aperm(tenseur, perm=c(2, 1, 3)), nrow=7) #e, g * t_max + t

unfolding3 <- matrix(aperm(tenseur, perm=c(3, 2, 1)), nrow=9000) #t, e * g_max + g

#SVD decomposition
decomp1 <- svd(unfolding1)
decomp2 <- svd(unfolding2)
decomp3 <- svd(unfolding3)

#Matrix product function
prod.mat <- function(matrice, vecteur){
	return(matrice %*% vecteur)
}

#Computation of the core tensor
kroneck <- decomp3$u %x% decomp2$u
core.tensor <- Conj(t(decomp1$u)) %*% unfolding1 %*% kroneck

sortie <- sort(Mod(core.tensor), decreasing=TRUE)

write.table(sortie, "code/resultatHOSVD3.txt", row.names=FALSE, col.names=FALSE)
