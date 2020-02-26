#############################################################################################################
# Author :
#   Jerome Mariette, MIAT, Universite de Toulouse, INRA 31326 Castanet-Tolosan France
#   Nathalie Villa-Vialaneix, MIAT, Universite de Toulouse, INRA 31326 Castanet-Tolosan France
#
# Copyright (C) 2017
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


.onAttach <- function(libname, pkgname) {
  
  packageStartupMessage("\nmixKernel will soon be included in the mixOmics package",
    "\n\nVisit http://www.mixOmics.org for a tutorial to use our method.",
    "\nAny bug reports or comments? Notify us at jerome.mariette@inrae.fr",
    " or https://bitbucket.org/klecao/package-mixomics/issues"
  )
  
}
