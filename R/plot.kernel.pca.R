#############################################################################################################
# Authors:
#   Kim-Anh Le Cao, French National Institute for Agricultural Research and ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Leigh Coonan, Queensland Faculty for Advanced Bioinformatics, Australia
#
# created: 2010
# last modified: 19-04-2016
#
# Copyright (C) 2010
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#############################################################################################################
#' @export

plot.kernel.pca <- function(x, ncomp = min(10, x$ncomp), ...) {
    
  #-- checking general input parameters --#

  if (!inherits(x, "kernel.pca")) {
  	stop("'x' should be an instance of 'kernel.pca' object returned by the",
  	     " kernel.pca function.", call. = FALSE)
  }
    
  #-- ncomp
  if (!is.numeric(ncomp) || ncomp < 1 || !is.finite(ncomp)) {
    stop("Invalid value for 'ncomp'.", call. = FALSE)
  }
  if (ncomp > length(x$sdev)) {
  	stop("'ncomp' must be lower or equal than ", length(x$sdev), ".", call. = FALSE)
  }
  
  barplot(x$prop_expl_var$X[1:ncomp], names.arg = seq(1, ncomp), 
          xlab = "Principal Components", ylab = "Explained Variance", ...)
    
}
