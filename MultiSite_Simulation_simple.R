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




p_hacking <- function(sim_heterogeneity_output){
  
  hacking_logical <- sapply(unique(sim_heterogeneity_output$sample), function(sample){
    
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
    
    
    t <- (mean1 - mean2) / (s_pooled / sqrt(n1+n2-2))
    
    p <- pt(t, df = (n1 + n2 - 2), lower.tail = FALSE)
    
    return(p < .05)
    
  })
  
  
  hacked_samples <- lapply(unique(sim_heterogeneity_output$sample), function(sample){
    
    if(hacking_logical[sample]){
    
      
      dat.sample <- sim_heterogeneity_output[sim_heterogeneity_output$sample == sample,]
      
      
      sample1 <- dat.sample$score[dat.sample$group == 1]
      sample2 <- dat.sample$score[dat.sample$group == 2]
      
      
      s1_sorted <- sort(sample1, decreasing = FALSE)
      s2_sorted <- sort(sample2, decreasing = TRUE)
      
      hacked_s1 <- sample1[sample1 %in% s1_sorted[-(1:10)]]
      hacked_s2 <- sample2[sample2 %in% s2_sorted[-(1:10)]]
      
      
      hacked_sample <- data.frame(scores = c(hacked_s1, hacked_s2),
                                  group = c(rep(1, length(hacked_s1)),
                                            rep(2, length(hacked_s2))),
                                  sample = sample)
      
      
    }else{
      hacked_sample <- sim_heterogeneity_output[sim_heterogeneity_output$sample == sample,]
    }
    
    
    
    return(hacked_sample)
    
  })
  
  return(data.frame(score = unlist(lapply(hacked_samples, FUN = function(x){x$score})),
                    group = unlist(lapply(hacked_samples, FUN = function(x){x$group})),
                    sample = unlist(lapply(hacked_samples, FUN = function(x){x$sample}))))
  
  
}



test <- sim_heterogeneity(mean_effect = 0)

ES.df <- SMD_simdata(test)

metafor::rma(measure = "GEN", yi = SMD, sei = SE_SMD, data = ES.df,
             method = "REML")



test_hacked <- p_hacking(test)

ES_hacked.df <- SMD_simdata(test_hacked)

metafor::rma(measure = "GEN", yi = SMD, sei = SE_SMD, data = ES_hacked.df,
             method = "REML")
