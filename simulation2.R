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

case_control_simulate <- function(f,n,sample_size,beta,alpha)
{
  allsnp <- rbinom(150*n*sample_size,1,f)+rbinom(150*n*sample_size,1,f)
  x <- matrix(allsnp,150*sample_size,n)
  G <- (x-2*f)/sqrt(2*f*(1-f))
  u <- G%*%beta
  p <- exp(alpha+u)/(1+exp(alpha+u))
  all_d <- apply(p,1,function(x){rbinom(1,1,x)})
  idx <- which(all_d==1)
  cases_id <- sample(idx,sample_size)
  cases <- cbind(u[cases_id],all_d[cases_id],G[cases_id,])
  jdx <- which(all_d==0)
  controls_id <- sample(jdx,sample_size)
  controls <- cbind(u[controls_id],all_d[controls_id],G[controls_id,])
  result <- rbind(controls,cases)
  return(result)
  
}
n <- 10
sample_size <- 10
f <- 0.25
alpha <- -4.6
sigma2 <- 2*log(2)

beta <- rnorm(n,0,sqrt(sigma2/(n)))





population <- case_control_simulate(f,n,sample_size,beta,alpha)


disease<- population[,2]
genetics <- population[,3:(n+2)]
u <- population[,1]

relation_id<- t(combn(disease,2))
relation_u <- relation_vector(u)
relation_d <- relation_vector(disease)
relation_g <- relation_matrix(genetics)

relation_matrix <- cbind(relation_id,relation_u,relation_d,relation_g)
relation_matrix <- as.data.frame(relation_matrix)
colnames(relation_matrix) <- c("id1","id2","u","d","g")

model1 <- lm(u~d,data=relation_matrix)
model2 <- lm(u~g,data=relation_matrix)
model3 <- lm(g~d,data=relation_matrix)


idx1 <- which(relation_matrix$id1==0&relation_matrix$id2==0)
idx2 <- which(relation_matrix$id1==1&relation_matrix$id2==1)
idx3 <- which(relation_d==1)
model4 <- lm(u~g-1,data=relation_matrix)
model5 <- lm(u[idx1]~g[idx1]-1,data=relation_matrix)
model6 <- lm(u[idx2]~g[idx2]-1,data=relation_matrix)
model7 <- lm(u[idx3]~g[idx3]-1,data=relation_matrix)


coeff=cbind(model1$coefficients,model2$coefficients,model3$coefficients)
test_coeff <- c(model4$coefficients,model5$coefficients,model6$coefficients,model7$coefficients)
save(relation_matrix,file=paste0("/dcl01/leek/nchatter/hzhang1/simulation1/relation_matrix_ideal",i1))
save(population,file=paste0("/dcl01/leek/nchatter/hzhang1/simulation1/population_ideal",i1))
save(coeff,file=paste0("/dcl01/leek/nchatter/hzhang1/simulation1/coeff_ideal",i1))
save(test_coeff,file=paste0("/dcl01/leek/nchatter/hzhang1/simulation1/test_coeff_ideal",i1))





