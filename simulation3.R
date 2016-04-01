rm(list=ls())
commanarg <- commandArgs(trailingOnly = T)
i1 <- as.numeric(commanarg[1])
library(Rcpp)
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
sigma2 <- 2*log(2)

beta <- rnorm(n,0,sqrt(sigma2/(n)))

geno_simulate_f<- function(n,sample_size,f){
  x <- n*sample_size
  result <- rbinom(x,1,f)+rbinom(x,1,f)
  result <- (result-2*f)/sqrt(2*f*(1-f))
  result_continous <- matrix(result,sample_size,n)
  return(result_continous)
}

geno <- geno_simulate_f(n,sample_size,f)

y <- geno%*%beta
relation_y <- relation_vector(y)
relation_g <- relation_matrix(geno)
model1 <- lm(relation_y~relation_g)
model2 <- lm(relation_y~relation_g-1)

result <- c(coef(model1)[1],coef(model1)[2],coef(model2))
save(result,file=paste0("result",i1,".Rdata"))
