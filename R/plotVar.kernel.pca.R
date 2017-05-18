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

plotVar.kernel.pca <- function(object, blocks = unique(object$cc.blocks), 
                               ndisplay = 5, ncol = 2, ...) {
  
  df.values <- df.variables <- df.blocks <- vector(length = ndisplay * length(blocks))
  i <- 1
  for (block in blocks) {
    values.block <- object$cc.distances[which(object$cc.blocks == block)]
    variables.block <- object$cc.variables[which(object$cc.blocks == block)]
    ordered.ids <- order(values.block, decreasing = TRUE)
    next.i <- i + ndisplay - 1
    df.values[i:next.i] <- values.block[ordered.ids][1:ndisplay]
    df.variables[i:next.i] <- variables.block[ordered.ids][1:ndisplay]
    df.blocks[i:next.i] <- rep(block, ndisplay)
    i <- i + ndisplay
  }

  df <- data.frame("variables" = df.variables, "values" = df.values,
                   "blocks" = df.blocks)
  df$variables <- reorder(df$variables, -df$values)

  ggplot(df, aes_string(x="variables", y = "values", fill="blocks")) +
    geom_bar(stat = "identity") + theme_bw() + 
    theme(axis.text.x = element_text(angle=45, hjust=1)) + ylab("") + xlab("") +
    facet_wrap(~ blocks, ncol = ncol, scales = "free_x") + 
    theme(legend.position = "none")
  
}
