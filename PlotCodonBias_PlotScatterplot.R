setwd("~/Documents/DeepSequencing")

#To plot codon bias from Test_uORF_Properties.html codon counts csv.
install.packages("ggplot2")
library(ggplot2)
codonbias = read.table("codonbias.csv",header=TRUE, stringsAsFactors = F)

codonbiasplot = matrix(ncol=3, nrow=length(codonbias[,1])*2, colnames = c("type", "freq", "codon"))
colnames(codonbiasplot) = c("type", "freq", "codon")
for (i in c(1:length(codonbias[,1]))){ 
  codonbiasplot[(i*2)-1,1] <- "uorf"
  codonbiasplot[(i*2),1] <- "CDS"
  codonbiasplot[(i*2)-1,2] <- codonbias[i,"uorf"]
  codonbiasplot[(i*2),2] <- codonbias[i,"CDS"]
  codonbiasplot[(i*2)-1,3] <- codonbias[i,"codon"]
  codonbiasplot[(i*2),3] <- codonbias[i,"codon"]
}

codonbiasplot[,"freq"] <- unlist(codonbiasplot[,"freq"])
codonbiasplot <- as.data.frame(codonbiasplot)
codonbiasplot$freq=as.numeric(levels(codonbiasplot$freq))[codonbiasplot$freq]

jpeg("CodonBias.jpg", width = 1500, height = 300)
ggplot(data=codonbiasplot, aes(x=codon, y=freq, fill=type)) +
  geom_bar(stat="identity", position=position_dodge()) + scale_y_continuous(name="Frequency", limits=c(0, .045))
dev.off()

#To plot one-to-one scatterplots between features and CDS_TE from previously made dataframes in linear regression RMD files. 
load("for_erica.dms")
uORFEff<-all_info

#Make an additional dataframe with only the transcripts with a single uORF. 
remove = c()
for (i in c(1:length(uORFEff[,1]))){ 
  match = which(uORFEff[,1] == uORFEff[i, 1])
  if (length(match) > 1){remove <- c(remove, i)
  }
  }
uORFEff_single <- uORFEff[which(!1:length(uORFEff[,1])%in%remove),]

#Find dominant uorfs
remove= c()
for (i in c(1:length(uORFEff[,1]))){ 
  match = which(uORFEff[,1] == uORFEff[i, 1])
  dom = c(which(uORFEff[match,"# Harr Reads"] == max(uORFEff[match,"# Harr Reads"])))
  match <- match[which(!1:length(match)%in%dom)]
  if (length(match) > 0){remove <- c(remove, match)
  i=i+length(match)-1
  }
}
uORFEff_dom <- uORFEff[which(!1:length(uORFEff[,1])%in%remove),]

#Make ECDFs for relevant properties. 

#Change main and name for each property of interest based on column names.
main = "uORF Frame vs. CDS Frame"
name = "Frame vs CDS"
jpeg(paste(main, ".jpg", sep=""))
plot(ecdf(uORFEff[,name]),do.points=T, main=main)
plot(ecdf(uORFEff_single[,name]), do.points=T, col="red", pch=16, add=T)
plot(ecdf(uORFEff_dom[,name]), do.points=T, col="blue", pch=16, add=T)
dev.off()

#Make scatterplot with line of best fit and correlation for dominant uORFs. One example:
jpeg("CDSTEvsuORFLen_dom.jpg", width = 500, height = 500)
plot(uORFEff_dom[,"Length"], uORFEff_dom[,"CDS_log2 TE"], main = "Log2(CDS TE) vs. uORF Length, Dominant uORF", xlab = "uORF Length", ylab = "Log2(CDS TE)")
abline(lm(uORFEff_dom[,"CDS_log2 TE"] ~ uORFEff_dom[,"Length"],  na.action=na.omit))
legend("topleft", bty="n", legend=cor.test(uORFEff_dom[,"CDS_log2 TE"], uORFEff_dom[,"Length"])$estimate)
dev.off()

#Make scatterplot with line of best fit and correlation for single uORFs. One example:
jpeg("CDSTEvsuORFLen_single.jpg", width = 500, height = 500)
plot(uORFEff_single[,"Length"], uORFEff_single[,"CDS_log2 TE"], main = "Log2(CDS TE) vs. uORF Length, Single uORF", xlab = "uORF Length", ylab = "Log2(CDS TE)")
abline(lm(uORFEff_single[,"CDS_log2 TE"] ~ uORFEff_single[,"Length"],  na.action=na.omit))
legend("topleft", bty="n", legend=cor.test(uORFEff_single[,"CDS_log2 TE"], uORFEff_single[,"Length"])$estimate)
dev.off()