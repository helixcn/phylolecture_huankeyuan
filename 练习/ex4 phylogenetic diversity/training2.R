## 用picante计算系统发育多样性
## 
setwd("/Users/jinlong/Desktop/beijing2019/peking university/ex4 phylogenetic diversity")
library(picante)

comm <- read.csv("grassland.community.csv", header = TRUE, row.names = 1)
phy <- read.tree("grassland.phylogeny.newick")
traits <- read.csv("species.traits.csv", header = TRUE, row.names = 1)
metadata <- read.csv("plot.metadata.csv", header = TRUE, row.names = 1)

comm.pd <- pd(comm, phy)
head(comm.pd)
phy.dist <- cophenetic(phy)

# ses.mpd
comm.sesmpd <- ses.mpd(comm, phy.dist, null.model = "richness", 
                       abundance.weighted = FALSE, 
                       runs = 999)
head(comm.sesmpd)

# ses.mntd
comm.sesmntd <- ses.mntd(comm, phy.dist, null.model = "richness",
                         abundance.weighted = FALSE,
                         runs = 999)
head(comm.sesmntd)
# 
comm.mntd.dist <- comdistnt(comm, phy.dist, 
                            abundance.weighted = TRUE)

# 更多内容，敬请关注 http://blog.sciencenet.cn/u/zjlcas
