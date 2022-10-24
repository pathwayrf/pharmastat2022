
#' Simulate trial data for multiple scenarios
#'
#' @param covObj an object of class \code{.covClass} generated in \code{\link{set_cov}}
#' @param seed the seed of R‘s random number generator. Default is 47
#' @param path file name for saving the output including folder path
#' @return a list of \code{matrix} containing simulated covariates information
#'
#' @export
#' @keywords method
setGeneric(name="simu_time", def=function(dt, eventObj, clinInt, clinExt, seed, path){standardGeneric("simu_time")})
setMethod(f="simu_time", signature(eventObj = ".eventClass"),
          definition=function(dt, eventObj, clinInt, clinExt, seed, path){

            if (missing(dt)) stop("Please provide a list of simulated data (dt).")
            if (missing(eventObj)) stop("Please provide eventObj.")
            if (missing(clinInt)) stop("Please provide clinInt.")
            if (missing(clinExt)) stop("Please provide clinExt.")
            if (missing(seed)){
              message("Set.seed(47)")
              seed = 47
            }

            seed_list <- seq(seed, seed + length(dt), by = 1)

            res_list <- lapply(seq(1, length(dt), by = 1), function(i){
              seed_i = seed_list[i]
              print(paste("seed =", seed_i))
              add_time(dt = dt[[i]], eventObj = eventObj, clinInt = clinInt, clinExt = clinExt, seed = seed_i)
            }) #loop foreach

            if (missing(path)) print("Simulated time are not saved.") else {
              save(res_list, file = path)
              print(paste("Simulated time are saved as", path))
            }
            res_list
          })
