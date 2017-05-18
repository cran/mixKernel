\name{TARAoceans}
\alias{TARAoceans}

\docType{data}

\title{TARA ocean microbiome data}

\description{
The TARA Oceans expedition facilitated the study of plankton communities by 
providing oceans metagenomic data combined with environmental measures to the 
scientific community. This dataset focuses on 139 prokaryotic-enriched samples 
collected from 68 stations and spread across three depth layers: the surface 
(SRF), the deep chlorophyll maximum (DCM) layer and the mesopelagic (MES) zones. 
Samples were located in height different oceans or seas: Indian Ocean (IO), 
Mediterranean Sea (MS), North Atlantic Ocean (NAO), North Pacific Ocean (NPO), 
Red Sea (RS), South Atlantic Ocean (SAO), South Pacific Ocean (SPO) and South 
Ocean (SO). Here, only a subset of the original data is provided (1\% of the 
35,650 prokaryotic operational taxonomic units (OTUs) and of the 39,246 
bacterial genes (NOGs) (selected at random).
}

\usage{data(TARAoceans)}

\format{A list containing the following components:
  \describe{
    \item{\code{phychem}}{data matrix with 139 rows and 22 columns. Each row 
    represents a sample and each column an environmental variable.}
    \item{\code{pro.phylo}}{data matrix with 139 rows (samples) and 356 columns
    (prokaryotic OTUs).}
    \item{\code{taxonomy}}{data matrix with 356 rows (prokaryotic OTUs) and 6 
    columns indicating the taxonomy of each OTU.}
    \item{\code{phylogenetic.tree}}{a phylo object (see package 'ape') 
    representing the prokaryotic OTUs phylogenetic tree.}
    \item{\code{pro.NOGs}}{data matrix with 139 rows (samples) and 638 columns 
    (NOGs).}
    \item{\code{sample}}{a list containing three following entries (all three 
    are character vectors): \code{name} (sample name), \code{ocean} (oceanic 
    region of the sample) and \code{depth} (sample depth).}
}}

\source{
The raw data were downloaded from \url{http://ocean-microbiome.embl.de/companion.html}.
}

\section{References}{
Sunagawa S., Coelho L.P., Chaffron S., Kultima J.R., Labadie K., Salazar F., 
Djahanschiri B., Zeller G., Mende D.R., Alberti A., Cornejo-Castillo F., 
Costea P.I., Cruaud C., d'Oviedo F., Engelen S., Ferrera I., Gasol J., Guidi L.,
Hildebrand F., Kokoszka F., Lepoivre C., Lima-Mendez G., Poulain J., Poulos B., 
Royo-Llonch M., Sarmento H., Vieira-Silva S., Dimier C., Picheral M., Searson 
S., Kandels-Lewis S., \emph{Tara} Oceans coordinators, Bowler C., de Vargas C.,
Gorsky G., Grimsley N., Hingamp P., Iudicone D., Jaillon O., Not F., Ogata H., 
Pesant S., Speich S., Stemmann L., Sullivan M., Weissenbach J., Wincker P., 
Karsenti E., Raes J., Acinas S. and Bork P. (2015). Structure and function of 
the global ocean microbiome. \emph{Science}, \bold{348}, 6237.
}

\keyword{datasets}
