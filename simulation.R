rm(list=ls())
commanarg <- commandArgs(trailingOnly = T)
i1 <- as.numeric(commanarg[1])

library(Rcpp)

relation_matrix_f <- function(disease,u,genetics){
  sample_size <- length(disease)/2
  relation_u <- NULL 
  relation_d <- NULL
  relation_g <- NULL
  for(i in 1:(2*sample_size-1)){
    rest_disease <- disease[(i+1):(2*sample_size)]
    temp <- (disease[i]-rest_disease)^2
    relation_d <- c(relation_d,temp)
    rest_u <- disease[(i+1):(2*sample_size)]
    temp <- (u[i]-rest_u)^2
    relation_u <- c(relation_u,temp)
    
  }
  for(i in 1:(2*sample_size-2)){
    rest_g <- genetics[(i+1):(2*sample_size),]
    temp_g <- sweep(rest_g,2,genetics[i,])
    temp <- apply(temp_g^2,1,mean)
    relation_g <- c(relation_g,temp)
  }
  i <- 2*sample_size-1
  rest_g <- genetics[(i+1):(2*sample_size),]
  temp_g <- rest_g-genetics[i,]
  temp <- mean(temp_g^2)
  relation_g <- c(relation_g,temp)
  result <- cbind(relation_u,relation_d,relation_g)
  return(result)
}

cppFunction('NumericVector relation_vector(NumericVector x){
            int n = x.size();
            NumericVector out(0.5*n*(n-1));
            int temp=0;
            for(int i=0;i< n-1;++i){
            for(int j=i+1;j<n;++j){
            out[temp] = pow(x[i]-x[j],2);
            temp +=1;
            }
            }
            return out;
            }')

cppFunction('double vector_distance(NumericVector x, NumericVector y){
            int n=x.size();
            double out=0;
            for(int i=0;i < n; ++i){
            out += pow(x[i]-y[i],2);
            }
            out = out/n;
            return out;
            }'
)

cppFunction('NumericVector relation_matrix(NumericMatrix x){
            int nrow=x.nrow(), ncol=x.ncol();
            NumericVector out(0.5*nrow*(nrow-1));
            int temp=0;
            double result=0;
            for(int i=0;i<nrow-1;++i){
            for(int j=i+1;j<nrow;++j){
            result=0;
            for(int k=0;k < ncol; ++k){
            result += pow(x(i,k)-x(j,k),2);
            }
            result = result/ncol;
            out[temp]=result;
            temp+=1;
            }
            }
            return out;
            }')


n <- 5000
sample_size <- 5000
f <- 0.25
alpha <- -4.6
sigma2 <- 2*log(2)

beta <- rnorm(n,0,sqrt(sigma2/(n)))
control_simulate_f<- function(n,sample_size,f){
  x <- n*sample_size
  result <- rbinom(x,1,f)+rbinom(x,1,f)
  result <- (result-2*f)/sqrt(2*f*(1-f))
  controls <- matrix(result,sample_size,n)
  return(controls)
}
controls <- control_simulate_f(n,sample_size,f)


case_simulate <- function(sample_size,beta,f){
  result <- NULL
  c1 <- exp(beta*2)*f^2
  c2 <- exp(beta)*2*f*(1-f)
  c3 <- (1-f)^2
  p1 <- c1/(c1+c2+c3)
  p2 <- c2/(c1+c2+c3)
  p3 <- 1-p1-p2
  p <- cbind(p1,p2,p3)
  result <- apply(p,1,onesnp_all)
  result <- (result-2*f)/sqrt(2*f*(1-f))
  return(result) 
 
}
onesnp_all <- function(y){
  G <- rmultinom(n=sample_size,size=1,prob=y)
  idx <- apply(G,2,function(x){which(x==1)})
  G <- 3-idx
  return(G)
}
cases <- case_simulate(sample_size,beta,f)

disease<- c(rep(0,sample_size),rep(1,sample_size))
genetics <- rbind(controls,cases)
u <- genetics%*%beta
population <- cbind(disease,u,genetics)
relation_u <- as.numeric(dist(u)^2)
relation_d <- as.numeric(dist(disease)^2)
relation_g <- as.numeric(dist(genetics)^2/n)

population <- cbind(u,disease,genetics)


relation_matrix <- cbind(relation_u,relation_d,relation_g)
relation_matrix <- as.data.frame(relation_matrix)
colnames(relation_matrix) <- c("u","d","g")

model1 <- lm(u~d,data=relation_matrix)
model2 <- lm(u~g,data=relation_matrix)
model3 <- lm(g~d,data=relation_matrix)




colnames(relation_matrix) <- c("u","d","g")
relation_matrix <- as.data.frame(relation_matrix)
coeff=cbind(model1$coefficients,model2$coefficients,model3$coefficients)
write.table(relation_matrix,file=paste0("relation_matrix",i1))
write.table(population,file=paste0("population",i1))
write.table(coeff,file=paste0("coeff",i1))

library(microbenchmark)
n <- 10
a <- matrix(rnorm(n*n),n,n)
dist1 <- function(a,n){
  as.numeric(dist(a)^2/n)
}
microbenchmark(
  dist1(a,n),
  relation_matrix(a)
  
)

time1 <- proc.time()
dist1(a,n)
time1 <- proc.time()-time1
