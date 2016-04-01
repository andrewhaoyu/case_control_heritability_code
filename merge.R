fileDir=getwd()
files=dir(fileDir,patter="coeff",full.names = T)
a1 <- NULL
s1 <- NULL
a2 <- NULL
s2 <- NULL
a3 <- NULL
s3 <- NULL
for(file in files){
	temp <- read.table(file)
	a1 <- c(a1,temp[1,1])
	s1 <- c(s1,temp[2,1])
	a2 <- c(a2,temp[1,2])
	s2 <- c(s2,temp[2,2])
	a3 <- c(a3,temp[1,3])
	s3 <- c(s3,temp[2,3])
}
library(data.table)
library(hexbin)
try <- fread("relation_matrix1")
plot_size <- 0.01*nrow(try)
idx <- sample(1:nrow(try),plot_size)
png("u_with_g.png",width=8,height=6,units="in",res=600)
plot(hexbin(try$V4,try$V2),xlab="g_ij",ylab="(u_i-u_j)^2")
dev.off()
png("u_with_d.png",width=8,height=6,units="in",res=600)
plot(hexbin(try$V3,try$V2),xlab="d_ij",ylab="(u_i-u_j)^2")
dev.off()
png("g_with_d.png",width=8,height=6,units="in",res=600)
plot(hexbin(try$V3,try$V4),xlab="d_ij",ylab="(g_i-g_j)^2")
dev.off()



fileDir <- getwd()
files <- dir(fileDir,patter="result",full.names=T)
a1 <- NULL
s1 <- NULL
s2 <- NULL
for(file in files){
	load(file)
	a1 <- c(a1,result[1])
	s1 <- c(s1,result[2])
	s2 <- c(s2,result[3])

}


fileDir <- getwd()
files <- dir(fileDir,patter="test_coeff",full.names=T)
s1 <- NULL
s2 <- NULL
s3 <- NULL
s4 <- NULL
for(file in files){
	load(file)
	s1 <- c(s1,test_coeff[1])
	s2 <- c(s2,test_coeff[2])
	s3 <- c(s3,test_coeff[3])
	s4 <- c(s4,test_coeff[4])
}
