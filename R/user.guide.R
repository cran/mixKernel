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