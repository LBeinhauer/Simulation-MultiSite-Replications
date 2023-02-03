# library loading and installing as necessary

packages <- c("metafor")

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



# function to simulate data
# General structure: 2 groups in each sample, mean_effect = mean1 - mean2 / sd_pooled
 # n participants in each group (therefore 2n in each sample), k samples drawn in total
 # heterogeneity in mean effect defined using tau
  # score scaling defined using pooled_sd and mu1

# CURRENTLY NO HETEROGENEITY IN POOLED_SD OR MEAN DIFFERENCE ALONE
# HETEROGENEITY INTRODUCED BY LETTING SMD VARY ACROSS SAMPLES

sim_heterogeneity <- function(n = 100, k = 100, mean_effect = .5, tau = .1,
                              pooled_sd = 10, mu1 = 10){
  
  # draw k sample effect sizes
  sample_effect_size <- rnorm(k, mean = mean_effect, sd = tau)
  
  # draw sample data, according to sampled effect sizes (loop usnig apply structure)
  sim_d.L <- lapply(as.matrix(1:k), FUN = function(sample){
    
    # extract effect size for this specific sample
    SMD <- sample_effect_size[sample]
    
    # find second mean using SMD, mu1 and pooled_SD
    mu2 <- mu1 - (SMD*pooled_sd)
    
    # randomly draw observations in group 1
    sample1 <- rnorm(n, mean = mu1, sd = pooled_sd)
    
    # randomly draw observations in group 2
    sample2 <- rnorm(n, mean = mu2, sd = pooled_sd)
    
    return(data.frame(score = c(sample1, sample2),
                      group = c(rep(1, length(sample1)),
                                rep(2, length(sample2))),
                      sample = sample))
  })
  
  
  # return a single data frame, containing in col 1 the scores, in col 2 the group indicator
   # and in col 3 the sample indicator
  return(data.frame(score = unlist(lapply(sim_d.L, FUN = function(x){x$score})),
                    group = unlist(lapply(sim_d.L, FUN = function(x){x$group})),
                    sample = unlist(lapply(sim_d.L, FUN = function(x){x$sample}))))
  
}






# function to extract standardized effect size (standardized mean difference) and respective
 # standard error from previously simulated data
SMD_simdata = function(sim_heterogeneity_output){
  
  # generate effect sizes for each sample respectively (using apply to loop across samples)
  ES <- lapply(unique(sim_heterogeneity_output$sample), FUN = function(sample){
    
    # specific sample data
    dat.sample <- sim_heterogeneity_output[sim_heterogeneity_output$sample == sample,]
    
    # separate groups
    sample1 <- dat.sample$score[dat.sample$group == 1]
    sample2 <- dat.sample$score[dat.sample$group == 2]
    
    # extract group means
    mean1 = mean(sample1)
    mean2 = mean(sample2)
    
    # extract group sample variances
    s2_1 = mean((sample1 - mean(sample1))^2)
    s2_2 = mean((sample2 - mean(sample2))^2)
    
    # extract group sizes
    n1 = length(sample1)
    n2 = length(sample2)
    
    # calculate pooled standard deviation
    s_pooled = sqrt(((n1-1)*s2_1 + (n2-1)*s2_2) / (n1 + n2 - 2))
    
    # caclulate SMD
    SMD = (mean1 - mean2)/s_pooled
    
    # generate squared standard error (variance) of SMD-estimate
    V_SMD = ((n1+n2) / (n1*n2)) + (SMD^2 / (2 * (n1 + n2)))
    
    # generate standard error of SMD-estimate
    SE_SMD = sqrt(V_SMD)
    
    return(list(SMD = SMD,
                SE_SMD = SE_SMD,
                sample = sample))
  })
  
  # again, return all estimates in a single data-frame.
  # col 1 contains SMD, col 2 contains standard error of SMD, col 3 contains sample indicator
  return(data.frame(SMD = unlist(lapply(ES, FUN = function(x){x$SMD})),
                    SE_SMD = unlist(lapply(ES, FUN = function(x){x$SE_SMD})),
                    sample = unlist(lapply(ES, FUN = function(x){x$sample}))))
  
}



# function to initiate a very simple standardized p-hacking procedure
# The function takes the simulated data and checks for each sample respectively, if a 
 # t-test would lead to significant values (at alpha = .05)
# If the t-test is NOT significant, the function will automatically remove the lowest 10
 # observation from group 1 and the highest 10 observations from group 2 (imagine a very
 # biased outlier-removal procedure)
# Subsequently, the function provides the same data format as it received
p_hacking <- function(sim_heterogeneity_output){
  
  # check wheter t-test would lead to a significant result (again looping across samples)
  hacking_logical <- sapply(unique(sim_heterogeneity_output$sample), function(sample){
    
    # extract respective sample
    dat.sample <- sim_heterogeneity_output[sim_heterogeneity_output$sample == sample,]
    
    # extract groups
    sample1 <- dat.sample$score[dat.sample$group == 1]
    sample2 <- dat.sample$score[dat.sample$group == 2]
    
    # extract means
    mean1 = mean(sample1)
    mean2 = mean(sample2)
    
    # extract sample variances
    s2_1 = mean((sample1 - mean(sample1))^2)
    s2_2 = mean((sample2 - mean(sample2))^2)
    
    # extract group sizes
    n1 = length(sample1)
    n2 = length(sample2)
    
    # calculate pooled standard deviation
    s_pooled = sqrt(((n1-1)*s2_1 + (n2-1)*s2_2) / (n1 + n2 - 2))
    
    # calculate t-value
    t <- (mean1 - mean2) / (s_pooled / sqrt(n1+n2-2))
    
    # check for p of t-value
    p <- pt(t, df = (n1 + n2 - 2), lower.tail = FALSE)
    
    # return whether the p-value is below .05
    return(p < .05)
    
  })
  
  
  # manipulate the samples which need hacking
  hacked_samples <- lapply(unique(sim_heterogeneity_output$sample), function(sample){
    
    # check if sample needs hacking
    if(hacking_logical[sample]){
    
      # extract sample
      dat.sample <- sim_heterogeneity_output[sim_heterogeneity_output$sample == sample,]
      
      # extract groups
      sample1 <- dat.sample$score[dat.sample$group == 1]
      sample2 <- dat.sample$score[dat.sample$group == 2]
      
      # sort values by their size
      s1_sorted <- sort(sample1, decreasing = FALSE)
      s2_sorted <- sort(sample2, decreasing = TRUE)
      
      # remove bottom 10 of group1, remove top 10 of group 2
      hacked_s1 <- sample1[sample1 %in% s1_sorted[-(1:10)]]
      hacked_s2 <- sample2[sample2 %in% s2_sorted[-(1:10)]]
      
      # combine values in data frame, similar to initial simulated data
      hacked_sample <- data.frame(scores = c(hacked_s1, hacked_s2),
                                  group = c(rep(1, length(hacked_s1)),
                                            rep(2, length(hacked_s2))),
                                  sample = sample)
      
      # if no hacking necesary, return initial simulated data
    }else{
      hacked_sample <- sim_heterogeneity_output[sim_heterogeneity_output$sample == sample,]
    }
    
    
    
    return(hacked_sample)
    
  })
  
  # again, return all data in a single dataframe
  # col 1 is the scores, col 2 the group indicator, col 3 the sample indicator
  return(data.frame(score = unlist(lapply(hacked_samples, FUN = function(x){x$score})),
                    group = unlist(lapply(hacked_samples, FUN = function(x){x$group})),
                    sample = unlist(lapply(hacked_samples, FUN = function(x){x$sample}))))
  
  
}

