install.packages("igraph")
install.packages("hash")
install.packages("magrittr")
install.packages("foreach")
install.packages("doParallel")
install.packages("tidyr")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("GGally") # ggnet2 is available through the GGally package.
install.packages("network") # ggnet2 depends on this package.
install.packages("sna")  # ggnet2 depends on this package.
install.packages("clue") # for hungarian algorithm.

source("https://bioconductor.org/biocLite.R")
biocLite("GOSemSim")
biocLite("org.Hs.eg.db")
biocLite("DOSE")
biocLite("clusterProfiler")