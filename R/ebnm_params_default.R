
#'@title default ebnm parameters
#'@export
ebnm_params_default = function(){
  return(list(prior_family='normal_scale_mixture',
              mode='estimate',
              scale = "estimate",
              g_init = NULL,
              fix_g = FALSE,
              output = output_default(),
              optmethod = NULL))
}
