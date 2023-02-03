
# library loading and installing as necessary

packages <- c("metafor", "ggplot2")

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



# load the functions defined in a separate R file
source(file.path("FUNCTIONS_MultiSite_Simulation_simple.R"))


## Apply the functions to simulate data once

test <- sim_heterogeneity(mean_effect = 0)

ES.df <- SMD_simdata(test)

metafor::rma(
             measure = "GEN"
             , yi = SMD
             , sei = SE_SMD
             , data = ES.df,
             method = "REML"
             )



test_hacked <- p_hacking(test)

ES_hacked.df <- SMD_simdata(test_hacked)

metafor::rma(
             measure = "GEN"
             , yi = SMD
             , sei = SE_SMD
             , data = ES_hacked.df
             , method = "REML"
             )





## preparing a larger sampling scheme, leading to multiple "replication projects"



prep.Effect_Sizes <- seq(
                         from = 0
                         , to = 2
                         , by = .1
                         )



prep.Tau_Variations <- seq(
                           from = 0
                           , to = 1
                           , by = .1
                           )

prep.Conditions <- expand.grid(prep.Effect_Sizes, prep.Tau_Variations)

names(prep.Conditions) <- c("ES", "tau")


# Replicate multiple replication projects (iteratively looping using lapply)
Multi_Replications <- lapply(
                             1:nrow(prep.Conditions)
                             , FUN = function(it){
                             
                              # simulating data from a single replication project (k = 100)
                              single.rep_proj <- sim_heterogeneity(n = 100
                                                                   , k = 100
                                                                   , mean_effect = prep.Conditions$ES[it]
                                                                   , tau = prep.Conditions$tau[it]
                                                                   , pooled_sd = 10
                                                                   , mu1 = 10)
                              }) # end of lapply


# Calculate standard errors and standardized mean differences (SMD) for each replication in 
 # all simulated replication projects
# Using lapply, we are looping across the list of replication projects
Multi_ES <- lapply(
                   Multi_Replications
                   , FUN = function(x){
                     
                     # apply function to generate SMD and standard error to simulated labs
                      # in single replication project
                     single.rep_SMDs <- SMD_simdata(x)
                   })



# use random-effects meta-analytic models to estimate a meta-analytic effect size of
 # SMD, and estimate its heterogeneity
# Again,we use lapply to iterator over the list of replication projects
Multi_MA_estimates <- lapply(
                             Multi_ES
                             , FUN = function(x){
                               
                               # use metafor to model a random-effects meta-analysis
                               single.rma <- metafor::rma(
                                                          measure = "GEN"
                                                          , yi = SMD
                                                          , sei = SE_SMD
                                                          , data = x
                                                          )
                               
                               # format the results
                               return(data.frame(
                                                 rma_est = single.rma$b[1]
                                                 , rma_tau = sqrt(single.rma$tau2)
                                                 )
                                      )
                               
                             })


# format the results from the random-effects meta-analysis in a single data frame
# add the "true" mean effect sizes and heterogeneity tau

df.Multi_MA_estimates <- data.frame(
                                    rma_est = unlist(lapply(Multi_MA_estimates
                                                            , FUN = function(x){x$rma_est}))
                                    ,rma_tau = unlist(lapply(Multi_MA_estimates
                                                             , FUN = function(x){x$rma_tau}))
                                    , sim_SMD = prep.Conditions$ES
                                    , sim_tau = prep.Conditions$tau
                                    )



# Visualise some of the results

# estimated vs simulated SMD
plot(df.Multi_MA_estimates$sim_SMD, df.Multi_MA_estimates$rma_est)
# estimated vs simulated tau
plot(df.Multi_MA_estimates$sim_tau, df.Multi_MA_estimates$rma_tau)

# estimated vs simulated SMD, differentiated by levels of heterogeneity
ggplot(
       df.Multi_MA_estimates
       , aes(x = sim_SMD
             , y = rma_est
             )
       ) +
  geom_point() +
  geom_smooth(method = "lm"
              , se = FALSE) +
  facet_wrap(df.Multi_MA_estimates$sim_tau)

# estimated vs simulated tau, differentiated by levels of effect size
ggplot(  
      df.Multi_MA_estimates
      , aes(x = sim_tau
            , y = rma_tau
            )
      ) +
  geom_point() +
  geom_smooth(
              method = "lm"
              , se = FALSE
              ) +
  facet_wrap(df.Multi_MA_estimates$sim_SMD)







## Lets try p-hacking the simulations done above

# looping along the simulated replication projects, we apply the p_hacking function
hack.Multi_Replications <- lapply(Multi_Replications
                                  , FUN = function(x){
                                    
                                    # apply p_hacking function to single replication proj.
                                    hacked.rep <- p_hacking(x)
                                  })




# Calculate standard errors and standardized mean differences (SMD) for each hacked replication 
# Using lapply, we are looping across the list of hacked creplication projects
hack.Multi_ES <- lapply(hack.Multi_Replications
                        , FUN = function(x){
                          
                          # apply function to generate SMD and standard error to simulated labs
                          # in single replication project
                          single.rep_SMDs <- SMD_simdata(x)
                        })


# use random-effects meta-analytic models to estimate a meta-analytic effect size of
# hacked SMD, and estimate its heterogeneity
# Again,we use lapply to iterator over the list of hacked replication projects
hack.Multi_MA_estimates <- lapply(hack.Multi_ES
                                  , FUN = function(x){
                                    
                                    # use metafor to model a random-effects meta-analysis
                                    single.rma <- metafor::rma(measure = "GEN"
                                                               , yi = SMD
                                                               , sei = SE_SMD
                                                               , data = x
                                    )
                                    
                                    # format the results
                                    return(data.frame(rma_est = single.rma$b[1]
                                                      , rma_tau = sqrt(single.rma$tau2)
                                                      )
                                    )
                                  })


# format the results from the random-effects meta-analysis in a single data frame
# add the "true" mean effect sizes and heterogeneity tau

df.hack.Multi_MA_estimates <- data.frame(rma_est = unlist(lapply(hack.Multi_MA_estimates
                                                                 , FUN = function(x){x$rma_est}))
                                         ,rma_tau = unlist(lapply(hack.Multi_MA_estimates
                                                                  , FUN = function(x){x$rma_tau}))
                                         , sim_SMD = prep.Conditions$ES
                                         , sim_tau = prep.Conditions$tau
)




# Visualise some of the results

# estimated vs simulated SMD
plot(df.hack.Multi_MA_estimates$sim_SMD
     , df.hack.Multi_MA_estimates$rma_est)
# estimated vs simulated tau
plot(df.hack.Multi_MA_estimates$sim_tau
     , df.hack.Multi_MA_estimates$rma_tau)

# estimated vs simulated SMD, differentiated by levels of heterogeneity
ggplot(df.hack.Multi_MA_estimates
       , aes(x = sim_SMD
             , y = rma_est
  )
) +
  geom_point() +
  geom_smooth(method = "lm"
              , se = FALSE) +
  facet_wrap(df.hack.Multi_MA_estimates$sim_tau)

# estimated vs simulated tau, differentiated by levels of effect size
ggplot(df.hack.Multi_MA_estimates
       , aes(x = sim_tau
        , y = rma_tau
        )
       ) +
  geom_point() +
  geom_smooth(method = "lm"
              , se = FALSE
              ) +
  facet_wrap(df.hack.Multi_MA_estimates$sim_SMD)





# compare hacked and non-hacked SMD estimates

# prepare data.frame, containing both hacked and non.hacked estimates
comb.df_MA_estimes <- data.frame(pure.rma_est = df.Multi_MA_estimates$rma_est
                                 , hack.rma_est = df.hack.Multi_MA_estimates$rma_est
                                 , pure.rma_tau = df.Multi_MA_estimates$rma_tau
                                 , hack.rma_tau = df.hack.Multi_MA_estimates$rma_tau
                                 , sim_SMD = prep.Conditions$ES
                                 , sim_tau = prep.Conditions$tau
           )



# plot non-hacked (x-axis) and hacked (y-axis) SMD estimates
plot(comb.df_MA_estimes$pure.rma_est
     , comb.df_MA_estimes$hack.rma_est)

# plot non-hacked (x-axis) and hacked (y-axis) SMD estimates - conditional on 
 # simulated SMD = 0 (p_hacking should have the strongest effect here)
plot(comb.df_MA_estimes$pure.rma_est[comb.df_MA_estimes$sim_SMD == 0]
     , comb.df_MA_estimes$hack.rma_est[comb.df_MA_estimes$sim_SMD == 0])


# plot non-hacked (x-axis) and hacked (y-axis) tau estimates
plot(comb.df_MA_estimes$pure.rma_tau
     , comb.df_MA_estimes$hack.rma_tau)

# plot non-hacked (x-axis) and hacked (y-axis) tau estimates - conditional on 
 # simulated SMD = 0 (p_hacking should have the strongest effect here)
plot(comb.df_MA_estimes$pure.rma_tau[comb.df_MA_estimes$sim_SMD == 0]
     , comb.df_MA_estimes$hack.rma_tau[comb.df_MA_estimes$sim_SMD == 0])



# plot non-hacked (x-axis) and hacked (y-axis) SMD estimates, conditional on sim. tau
ggplot(comb.df_MA_estimes
       , aes(x = pure.rma_est
             , y = hack.rma_est)
       ) +
  geom_point() +
  # geom_smooth(method = "lm"
  #             , se = FALSE
  #             ) +
  facet_wrap(vars(sim_tau))

# plot non-hacked (x-axis) and hacked (y-axis) tau estimates, conditional on sim. SMD
ggplot(comb.df_MA_estimes
       , aes(x = pure.rma_tau
             , y = hack.rma_tau)
) +
  geom_point() +
  # geom_smooth(method = "lm"
  #             , se = FALSE
  #             ) +
  facet_wrap(vars(sim_SMD))

