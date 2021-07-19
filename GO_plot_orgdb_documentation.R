## Xenopus
library(dplyr)
library(BiocManager)
library(clusterProfiler)
library(AnnotationHub)

hub<-AnnotationHub()
query(hub, "Xenopus tropicalis")
org.xenopus.eg.db<-hub[["AH94382"]]

## S. pombe
library(AnnotationForge)
library(biomaRt)
library(GO.db)
library(devtools)

if (! "org.Spombe.eg.db" %in% installed.packages()) {
  orgdb <- AnnotationForge::makeOrgPackageFromNCBI(
    version = "0.1", author = "atb <abelew@gmail.com>",
    maintainer = "atb <abelew@gmail.com>", tax_id = "4896",
    genus = "Schizosaccharomyces", species = "pombe")
  devtools::install_local("org.Spombe.eg.db")
}