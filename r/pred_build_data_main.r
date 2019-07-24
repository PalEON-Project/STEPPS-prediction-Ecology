
res  = 3

version="v2.4"
AR     = TRUE
draw   = FALSE
ndraws = 10
nknots = 120
one_time=FALSE
lambda_fixed = TRUE
bchron = FALSE

# stat model flags
decomp     = TRUE
bt         = TRUE
mpp        = TRUE
save_plots = TRUE

bacon = TRUE
# draw  = TRUE
add_varves = TRUE
constrain  = FALSE
# how far to extrapolate past last geochron anchor
nbeyond = 5000

age_model  = 'bacon'

suff_dat = '12taxa_mid_comp_ALL_v0.3'


run_pl_Ka_Kgamma = list(suff_fit  = 'cal_pl_Ka_Kgamma_EPs_ALL_v0.4c1', 
                        suff_dat  = suff_dat,
                        kernel    = 'pl', 
                        one_a     = FALSE, 
                        one_b     = TRUE,
                        one_gamma = FALSE, 
                        EPs       = TRUE)

runs = list(run_pl_Ka_Kgamma)

# for (run in runs){
  # for (res in grids){
#     if (draw){
#       for (dr in 1:ndraws){
#         source('r/pred_build_data.r')
#       }
#     }
  # }
# }

if (draw){
  for (dr in 1:ndraws){
    source('r/pred_build_data.r')
  }
} else {
  dr = 1
  source('r/pred_build_data.r')
}
