rm(list=ls())
commanarg <- commandArgs(trailingOnly = T)
i1 <- as.numeric(commanarg[1])



load(paste0("/dcl01/leek/nchatter/hzhang1/simulation1/result/relation_matrix_ideal",i1))
model <- lm(g~d,data=relation_matrix)
beta <- c(model$coefficients[1],model$coefficients[2]*5000)
save(beta,file=paste0("/dcl01/leek/nchatter/hzhang1/simulation1/result/sigma",i1))
