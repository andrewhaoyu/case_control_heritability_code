sigma2 <- 2*log(2)
n <- 5000
population_size <- 60000
beta <- rnorm(n,0,sd = sqrt(sigma2/n))
f <- 0.25
x <- matrix(rbinom(n*population_size,2,f),population_size,n)
G <- (x-2*f)/sqrt(2*f*(1-f))
alpha <- -5.3
u <- G%*%beta
p <- exp(alpha+u)/(1+exp(alpha+u))
disease <- sapply(p,function(x){rbinom(1,1,x)})
population <- cbind(disease,u,p,G)
sample_size <- 500
case_control_sample_f <- function(population,sample_size){
  population_case <- which(population[,1]==1)
  cases <- population[sample(population_case,sample_size),]
  population_control <- which(population[,1]==0)
  controls <- population[sample(population_control,sample_size),]
  result <- rbind(cases,controls)
  return(result)
}
case_control_sample <- case_control_sample_f(population,sample_size)
population <- as.data.frame(population)
random_effect <- matrix(0,0.5*2*sample_size*(2*sample_size-1),4)
temp <- 1
for(i in 1:(2*sample_size-1)){
  for(j in (i+1):(2*sample_size)){
    random_effect[temp,1] <- (case_control_sample[i,2]-case_control_sample[j,2])^2
    random_effect[temp,2] <- (case_control_sample[i,1]-case_control_sample[j,1])^2
    random_effect[temp,3] <- mean((case_control_sample[i,4:(n+3)]-case_control_sample[j,4:(n+3)])^2)
    random_effect[temp,4] <- mean(beta^2%*%(case_control_sample[i,4:(n+3)]-case_control_sample[j,4:(n+3)])^2)
    temp <- temp+1
  }
}
colnames(random_effect) <- c("u","d","gd","g")
random_effect <- as.data.frame(random_effect)
model1 <- lm(u~d,data=random_effect)
model2 <- lm(u~gd-1,data=random_effect)
model3 <- lm(gd~d,data = random_effect)
model4 <- lm(gd~d,data = random_effect)




