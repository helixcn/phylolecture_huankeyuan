### 第一部分，用RAxML建树
# 数据来源
# Kress, W. J., Erickson, D. L., Jones, F. A., Swenson, N. G., Perez, R., Sanjur, O., & Bermingham, E. (2009). Plant DNA barcodes and a community phylogeny of a tropical forest dynamics plot in Panama. Proceedings of the National Academy of Sciences, 106(44), 18621–18626.
# https://www.pnas.org/content/106/44/18621/tab-figures-data

### 从GenBank下载序列，用muscle比对，创建supermatrix，及RAxML的入门操作

# MacOS
setwd("/Users/jinlong/Desktop/beijing2019/huankeyuan20190916/ex1 data preparation and RAxML")
# Windows
setwd("C:\\Users\\sj\\Desktop\\beijing2019\\excercise")

# 读取论文附件中的Genbank Accession Number
library(openxlsx)
library(ape)
dat <- read.xlsx("BCI_barcodes.xlsx")
head(dat)

# 获取序列
# matk <- read.GenBank(access.nb = dat$matK)
# rbcLa <- read.GenBank(access.nb = dat$rbcLa)
# trnH_psbA <- read.GenBank(access.nb = dat$trnH_psbA)

# 保存序列
# write.FASTA(matk, "matk.fas")
# write.FASTA(rbcLa, "rbcLa.fas")
# write.FASTA(trnH_psbA, "trnH_psbA.fas")

# 重命名序列
# devtools::install_github("helixcn/phylotools")
library(phylotools)
matk      <- read.fasta("matk.fas") 
rbcla     <- read.fasta("rbcLa.fas")
trnH_psbA <- read.fasta("trnH_psbA.fas")

library(plantlist)
# 提取学名，不带命名人
res_taxa <- parse_taxa(dat$Taxon) 
# 增加“干净”的学名
dat$species <- gsub("__", "", paste(res_taxa$GENUS_PARSED, res_taxa$SPECIES_PARSED, res_taxa$INFRASPECIFIC_RANK_PARSED, res_taxa$INFRASPECIFIC_EPITHET_PARSED, sep = "_"))


# 保存一个副本，看匹配是否正确
# Save a copy, and examine
# write.xlsx(dat, "dat_xxxxx.xlsx")
# matK	rbcLa	trnH_psbA

# 生成序列名和物种对照表，用于序列重命名
matk_tab <- subset(dat, select = c("matK", "species"))
rbcla_tab <- subset(dat, select = c("rbcLa", "species"))     
trnH_psbA_tab  <- subset(dat, select = c("trnH_psbA", "species"))

# 序列重命名
rename.fasta("matk.fas", matk_tab, "matk_renamed.fasta") 
rename.fasta("rbcLa.fas", rbcla_tab, "rbcla_renamed.fasta")
rename.fasta("trnH_psbA.fas", trnH_psbA_tab, "trnH_psbA_renamed.fasta")

# system("mafft trnH_psbA_renamed.fasta > aligned_trnH_psbA_renamed.fasta")



# 每个片段，分别用Mafft比对序列
# （请先安装好mafft）
system("mafft matk_renamed.fasta > aligned_matk_renamed.fasta")
system("mafft rbcla_renamed.fasta > aligned_rbcla_renamed.fasta")

# 每个片段，分别用MUSCLE比对序列
# （请先将MUSCLE.exe）拷贝到工作文件夹
system("muscle.exe -in matk_renamed.fasta -out aligned_matk_renamed.fasta")
system("muscle.exe -in rbcla_renamed.fasta -out aligned_rbcla_renamed.fasta")

# 用bioedit或者AliView进行校对
# (本步操作从略）

# 用phylotools建立supermatrix 
supermat(infiles = c("aligned_matk_renamed.fasta", "aligned_rbcla_renamed.fasta"), outfile = "BCI2.phy")
         
# 用RAxML建树，获取BestTree （因RAxML的参数组合十分复杂，请注意命令选取）
## MacOS 注意， 这步在普通台式机上也最少要运行2h以上
system("./raxml -f a -x 12345 -p 12345 -# 100 -m GTRGAMMA -s BCI2.phy -n BCI_rbcla_matk_besttree -T 2")

## Windows 注意， 这步在普通台式机上也最少要运行2h以上
system("raxmlHPC.exe -f a -x 12345 -p 12345 -# 100 -m GTRGAMMA -s BCI.phy -n BCI_rbcla_matk_besttree")

## 若在本地机上建树太慢，可通过 The CIPRES Science Gateway V. 3.3 服务器上建树
### 该网站提供各软件的程序界面，可以选择相应参数，由于是在服务器上并行计算，速度很快。
## http://www.phylo.org/

# 用ape或Figtree将进化树保存为pdf文件
# 查看Best Tree的拓扑结构，核对进化树的准确性

# 用iqtree建树 （包括了模型筛选）
setwd("/Users/jinlong/Desktop/beijing2019/excercise")
system("iqtree -s BCI.phy")

# 用Paup*、Phylip、BEAST、MrBayes建树 从略

##################################################
########### 上述操作有哪些问题？ #####################
1. 建树时没有指定外类群
2. 没有root
3. 未经过模型筛选
4. 未使用分隔模型
5. 快速Bootstrap数量够吗？


#####################################################
################### APE 的基本操作#####################
#####################################################
#####################################################

# 用ape或者FigTree设定外类群（reroot），tree ladderize
setwd("/Users/jinlong/Desktop/beijing2019/peking university/ex1 data preparation and RAxML")
library(ape)
example_tree <- read.tree("example.tre") # 由Phylomatic生成的进化树，当然，这是一棵Ultrametric树
plot(example_tree)

plot(ladderize(example_tree, right = TRUE)) # 类群排序
plot(ladderize(example_tree, right = FALSE))

# 保存排序的树
laddered_tree <- ladderize(example_tree, right = FALSE)

# 重新设定外类群，看是否合理
rerooted_tree <- root(ladderize(example_tree, right = FALSE), "Pinus_massoniana")

# 查看拓扑结构
plot(rerooted_tree) # 这棵树有什么问题？ 外类群选取不当

# 用ape保存为nexus格式的进化树
write.nexus(laddered_tree, file = "laddered_tree.nex")

######################################################
############ ex2 r8s #################################
######################################################

# r8s的安装和使用方法
http://blog.sciencenet.cn/blog-255662-1144730.html
http://blog.sciencenet.cn/blog-255662-305898.html


# 将RAxML的best tree读取到R中
# read.tree
# 保存为nexus格式
# write.nexus

# 用Notepad++，在nexus进化树后，加上r8s 模块


# r8s模块的编辑，选取进化树内的最近共同祖先
# 在http://www.timetree.org/查询分化时间
# 最好是根据化石查询

# r8s软件的运行
# 查看 r8s.bat

# 从r8s的结果拷贝出newick文件，另存为BCI_dated_phylo.tre
