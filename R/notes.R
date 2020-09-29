stlm_construct <- function(estmethod) {
  new_estmethod <- structure(list(), class = estmethod)
}

stlm <- function(t) {
  UseMethod("stlm", object = t)
}

stlm.loglik <- function(t) {
  x <- TRUE
  return(x)
}

stlm("loglik")

stlm.gam <- function()
stlm.svls <- function()


# i think i should just use switch
# stlm_sv, stlm_loglik, stlm_gam
