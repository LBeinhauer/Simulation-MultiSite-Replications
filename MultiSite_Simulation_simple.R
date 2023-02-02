# library loading and installing as necessary

packages <- c("metafor", "here")

# check, whether library already installed or not - install and load as needed:
apply(as.matrix(packages), MARGIN = 1, FUN = function(x) {
  
  pkg_avail <- nzchar(system.file(package = x))   # check if library is installed on system
  
  if(pkg_avail){
    require(x, character.only = TRUE)             # load the library, if already installed
    
  }else{
    install.packages(x)                           # install the library, if missing
    require(x, character.only = TRUE)             # load after installation
  }
})




sim_heterogeneity <- function(n = 100, k = 100, mean_effect = .5, tau = .1,
                              pooled_sd = 10, mu1 = 10){
  
  sample_effect_size <- rnorm(k, mean = mean_effect, sd = tau)
  
  
  sim_d.L <- lapply(as.matrix(1:k), FUN = function(sample){
    
    SMD <- sample_effect_size[sample]
    
    mu2 <- mu1 - (SMD*pooled_sd)
    
    sample1 <- rnorm(n, mean = mu1, sd = pooled_sd)
    
    sample2 <- rnorm(n, mean = mu2, sd = pooled_sd)
    
    return(data.frame(score = c(sample1, sample2),
                      group = c(rep(1, length(sample1)),
                                rep(2, length(sample2))),
                      sample = sample))
  })
  
  
  
  return(data.frame(score = unlist(lapply(sim_d.L, FUN = function(x){x$score})),
                    group = unlist(lapply(sim_d.L, FUN = function(x){x$group})),
                    sample = unlist(lapply(sim_d.L, FUN = function(x){x$sample}))))
  
}







SMD_simdata = function(sim_heterogeneity_output){
  
  ES <- lapply(unique(sim_heterogeneity_output$sample), FUN = function(sample){
    
    dat.sample <- sim_heterogeneity_output[sim_heterogeneity_output$sample == sample,]
    
    sample1 <- dat.sample$score[dat.sample$group == 1]
    sample2 <- dat.sample$score[dat.sample$group == 2]
    
    mean1 = mean(sample1)
    mean2 = mean(sample2)
    
    s2_1 = mean((sample1 - mean(sample1))^2)
    s2_2 = mean((sample2 - mean(sample2))^2)
    
    n1 = length(sample1)
    n2 = length(sample2)
    
    s_pooled = sqrt(((n1-1)*s2_1 + (n2-1)*s2_2) / (n1 + n2 - 2))
    
    SMD = (mean1 - mean2)/s_pooled
    
    V_SMD = ((n1+n2) / (n1*n2)) + (SMD^2 / (2 * (n1 + n2)))
    
    SE_SMD = sqrt(V_SMD)
    
    return(list(SMD = SMD,
                SE_SMD = SE_SMD,
                sample = sample))
  })
  
  return(data.frame(SMD = unlist(lapply(ES, FUN = function(x){x$SMD})),
                    SE_SMD = unlist(lapply(ES, FUN = function(x){x$SE_SMD})),
                    sample = unlist(lapply(ES, FUN = function(x){x$sample}))))
  
}




test <- sim_heterogeneity()



ES.df <- SMD_simdata(test)
metafor::rma(measure = "GEN", yi = SMD, sei = SE_SMD, data = ES.df,
             method = "REML")



