#' View mixKernel User's Guide
#' 
#' Find the location of the mixKernel User's Guide and optionnaly opens it
#'
#' @param html logical. Should the document returned by the function be the
#' compiled PDF or the Rmd source. Default to \code{TRUE}
#' @param view logical. Should the document be opened using the default HTML
#' viewer? Default to \code{html}. It has no effect if \code{html = FALSE}
#' 
#' @return Character string giving the file location. If \code{html = TRUE} and
#' \code{view = TRUE}, the HTML document reader is started and the User's Guide
#' is opened in it.
#' 
#' @details
#' If the operating system is not Windows, then the HTML viewer used is that
#' given by \code{Sys.getenv("R_BROWSER")}. The HTML viewer can be changed using
#' \code{Sys.setenv(R_BROWSER = )}.
#' 
#' 
#' @author Jerome Mariette <jerome.mariette@@inrae.fr>
#' Nathalie Vialaneix <nathalie.vialaneix@@inrae.fr>
#' @export
#' @examples
#' mixKernel.users.guide(view = FALSE)
#' mixKernel.users.guide(html = FALSE)
#' \dontrun{mixKernel.users.guide()}
#' 
mixKernel.users.guide <- function(html = TRUE, view = html) {
  if (html) {
    f <- system.file("doc", "mixKernelUsersGuide.html", package = "mixKernel")
    if (view) {
      if (.Platform$OS.type == "windows")
        shell.exec(f)
      else browseURL(paste0("file://", f))
    }
  } else {
    f <- system.file("doc", "mixKernelUsersGuide.Rmd", package = "mixKernel")
    if (view) {
      warning("'mixKernelUsersGuide.Rmd' can not be viewed.
              However, the location of the file is returned by the function.",
              call. = FALSE)
    }
  }
  return(f)
}