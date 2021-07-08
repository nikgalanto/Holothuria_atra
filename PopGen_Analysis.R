# to install strataG we had to comment out the FLIBS variable in Makeconf file located at /Library/Frameworks/R.framework/Resources/etc/Makeconf
# Following instructions in the second answer found here https://stackoverflow.com/questions/23916219/os-x-package-installation-depends-on-gfortran-4-8
devtools::install_github("stranda/rmetasim")
devtools::install_github('ericarcher/strataG', build_vignettes = TRUE)
library(strataG)
library(rentrez)
library(pegas)
library(ape)
library(stringr)
library(geomedb)
coelestis<-querySanger(marker="CR”,query="+genus:Pomacentrus+specificEpithet:coelestis)
coelestis<-querySanger(marker="CR",query=""+genus:Pomacentrus+specificEpithet:coelestis)
coelestis<-querySanger(marker="CR”,query="+genus:Pomacentrus+specificEpithet:coelestis)
coelestis<-querySanger(marker="CR”,query="+genus:Pomacentrus+specificEpithet:coelestis)
pcoel<-querySanger(locus="CR",query="genus = Pomacentrus AND specificEpithet = coelestis")
image.DNAbin(pcoel)
pcoel_aln <- muscle(pcoel)
image.DNAbin(pcoel_aln)
pcoel_aln_trim <- pcoel_aln[,25:360]
pcoel <- pcoel_aln_trim
image.DNAbin(pcoel)

#haplotype network lesson

pop<-gsub(names(pcoel), pattern = "_\\d+", replacement="")
d <-dist.dna(pcoel)
h <-pegas::haplotype(pcoel)
h <- sort(h, what = "label")
net <-pegas::haploNet(h)
i <-stack(setNames(attr(h, "index"), rownames(h)))
i <-i[order(i$values),]
ind.hap <-table(hap=i$ind, pop=pop)
length(pop)
plot(net,size=attr(net, "freq"), scale.ratio=10, pie=ind.hap, legend=F, labels=F, threshold=0, show.mutation=2)
names(pcoel)
rownames(pcoel)
library(stringr)
ind.names <-str_replace(string=rownames(pcoel), pattern=".+\\[tissueID = (.+?)\\].+", replacement = "\\1")
