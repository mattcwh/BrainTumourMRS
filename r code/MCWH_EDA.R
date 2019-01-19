library(XML)
library(car)
library(MASS)
library(ggplot2)
### explore xml file ###
xmlfile = xmlTreeParse("1TE-ShortMRS(corregido).xml")
class(xmlfile)
topxml = xmlRoot(xmlfile)
class(topxml)
xmlName(topxml)
xmlSize(topxml)
xmlAttrs(topxml)
topxml[[1]]
xmlName(topxml[[1]])
xmlAttrs(topxml[[1]])
xmlSApply(topxml[[1]], xmlName) 
xmlSApply(topxml[[1]], xmlAttrs)
xmlSApply(topxml[[1]], xmlSize) 
xmlName(topxml[[1]][[1]])
xmlAttrs((topxml[[1]][[1]]))
xmlName(topxml[[1]][[2]])
xmlSApply(topxml[[1]][[2]],xmlName)
topxml[[1]][[2]][['Parameters']]
topxml[[1]][[2]][['Points']]

###################################################

# load xml file
xmlfile1= xmlParse('1TE-ShortMRS(corregido).xml')   
class(xmlfile1)
# extract id
id =  do.call(rbind,xpathApply(xmlfile1,"//DATASET/Case", xmlAttrs)) 
class(id)
head(id)
# extract spectra points (to a list of strings)
points1=xpathApply(xmlfile1,"//DATASET/Case/Spectrum/Points", xmlValue)
class(points1)
# split and convert spectra points into matrix
pointmatrix=matrix(0,304,512)
for (i in 1:304) 
{pointmatrix[i,] = as.numeric(unlist(strsplit(points1[[i]]," ")))}
pointmatrix[2,]
points1[[2]] # check if converted correctly
# combine id and spectra 
spectra = data.frame(id,pointmatrix)
# extract tumour type
type = do.call(rbind,xpathApply(xmlfile1,"//DATASET/Case/Tissue", xmlAttrs))
# combine type and spectra
spectra = data.frame(type,pointmatrix,stringsAsFactors=FALSE)

#####################################################

# subset GL spectra
GLspectra = subset(spectra,Type=='gl')
# calculate df of spectra median and corresponding ppm
GLmedian = data.frame(ppm = seq(-2.7, 7.1, length.out=512), median = sapply(GLspectra[, 2:513], median))
head(GLmedian)
# add median absolute deviation to the df
GLmedian$mad = sapply(GLspectra[, 2:513], mad)
# add upper and lower MAD
GLmedian$upperMAD = GLmedian$median + GLmedian$mad
GLmedian$lowerMAD = GLmedian$median - GLmedian$mad

# plot averaged spectrum
ggplot(GLmedian, aes(x=ppm, y=median)) + geom_line(size=1, alpha=0.8) + geom_ribbon(aes(ymin=GLmedian$lowerMAD, ymax=GLmedian$upperMAD), fill="black", alpha="0.3") + xlab("frequency (ppm)") + ylab("signal (arbituary unit)") + scale_x_reverse(breaks=seq(-2.5, 7.0, 0.5)) + scale_y_continuous(limits=c(-6, 35), position="right") + theme_classic(base_size=16) + ggtitle("glioblastoma") + theme(plot.title=element_text(hjust=0.5, size=30), panel.grid.major=element_line(), panel.border=element_rect(fill=NA, size=2))

#####################################################

# MM - meningioma
MMspectra = subset(spectra,Type=='mm')
MMmedian = data.frame(ppm = seq(-2.7, 7.1, length.out=512), median = sapply(MMspectra[, 2:513], median))
MMmedian$mad = sapply(MMspectra[, 2:513], mad)
MMmedian$upperMAD = MMmedian$median + MMmedian$mad
MMmedian$lowerMAD = MMmedian$median - MMmedian$mad
ggplot(MMmedian, aes(x=ppm, y=median)) + geom_line(size=1, alpha=0.8) + geom_ribbon(aes(ymin=MMmedian$lowerMAD, ymax=MMmedian$upperMAD), fill="black", alpha="0.3") + xlab("frequency (ppm)") + ylab("signal (arbituary unit)") + scale_x_reverse(breaks=seq(-2.5, 7.0, 0.5)) + scale_y_continuous(limits=c(-6, 35), position="right") + theme_classic(base_size=16) + ggtitle("meningioma") + theme(plot.title=element_text(hjust=0.5, size=30), panel.grid.major=element_line(), panel.border=element_rect(fill=NA, size=2))

######################################################

# ME - metastasis
MEspectra = subset(spectra,Type=='me')
MEmedian = data.frame(ppm = seq(-2.7, 7.1, length.out=512), median = sapply(MEspectra[, 2:513], median))
MEmedian$mad = sapply(MEspectra[, 2:513], mad)
MEmedian$upperMAD = MEmedian$median + MEmedian$mad
MEmedian$lowerMAD = MEmedian$median - MEmedian$mad
ggplot(MEmedian, aes(x=ppm, y=median)) + geom_line(size=1, alpha=0.8) + geom_ribbon(aes(ymin=MEmedian$lowerMAD, ymax=MEmedian$upperMAD), fill="black", alpha="0.3") + xlab("frequency (ppm)") + ylab("signal (arbituary unit)") + scale_x_reverse(breaks=seq(-2.5, 7.0, 0.5)) + scale_y_continuous(limits=c(-6, 35), position="right") + theme_classic(base_size=16) + ggtitle("metastasis") + theme(plot.title=element_text(hjust=0.5, size=30), panel.grid.major=element_line(), panel.border=element_rect(fill=NA, size=2))

#######################################################

# NO - normal
NOspectra = subset(spectra,Type=='no')
NOmedian = data.frame(ppm = seq(-2.7, 7.1, length.out=512), median = sapply(NOspectra[, 2:513], median))
NOmedian$mad = sapply(NOspectra[, 2:513], mad)
NOmedian$upperMAD = NOmedian$median + NOmedian$mad
NOmedian$lowerMAD = NOmedian$median - NOmedian$mad
ggplot(NOmedian, aes(x=ppm, y=median)) + geom_line(size=1, alpha=0.8) + geom_ribbon(aes(ymin=NOmedian$lowerMAD, ymax=NOmedian$upperMAD), fill="black", alpha="0.3") + xlab("frequency (ppm)") + ylab("signal (arbituary unit)") + scale_x_reverse(breaks=seq(-2.5, 7.0, 0.5)) + scale_y_continuous(limits=c(-6, 35), position="right") + theme_classic(base_size=16) + ggtitle("normal") + theme(plot.title=element_text(hjust=0.5, size=30), panel.grid.major=element_line(), panel.border=element_rect(fill=NA, size=2))

########################################################

A2spectra = subset(spectra,Type=='a2')
A2median = data.frame(ppm = seq(-2.7, 7.1, length.out=512), median = sapply(A2spectra[, 2:513], median))
A2median$mad = sapply(A2spectra[, 2:513], mad)
A2median$upperMAD = A2median$median + A2median$mad
A2median$lowerMAD = A2median$median - A2median$mad

A3spectra = subset(spectra,Type=='a3')
A3median = data.frame(ppm = seq(-2.7, 7.1, length.out=512), median = sapply(A3spectra[, 2:513], median))
A3median$mad = sapply(A3spectra[, 2:513], mad)
A3median$upperMAD = A3median$median + A3median$mad
A3median$lowerMAD = A3median$median - A3median$mad

ABspectra = subset(spectra,Type=='ab')
ABmedian = data.frame(ppm = seq(-2.7, 7.1, length.out=512), median = sapply(ABspectra[, 2:513], median))
ABmedian$mad = sapply(ABspectra[, 2:513], mad)
ABmedian$upperMAD = ABmedian$median + ABmedian$mad
ABmedian$lowerMAD = ABmedian$median - ABmedian$mad

HBspectra = subset(spectra,Type=='hb')
HBmedian = data.frame(ppm = seq(-2.7, 7.1, length.out=512), median = sapply(HBspectra[, 2:513], median))
HBmedian$mad = sapply(HBspectra[, 2:513], mad)
HBmedian$upperMAD = HBmedian$median + HBmedian$mad
HBmedian$lowerMAD = HBmedian$median - HBmedian$mad

LYspectra = subset(spectra,Type=='ly')
LYmedian = data.frame(ppm = seq(-2.7, 7.1, length.out=512), median = sapply(LYspectra[, 2:513], median))
LYmedian$mad = sapply(LYspectra[, 2:513], mad)
LYmedian$upperMAD = LYmedian$median + LYmedian$mad
LYmedian$lowerMAD = LYmedian$median - LYmedian$mad

OAspectra = subset(spectra,Type=='oa')
OAmedian = data.frame(ppm = seq(-2.7, 7.1, length.out=512), median = sapply(OAspectra[, 2:513], median))
OAmedian$mad = sapply(OAspectra[, 2:513], mad)
OAmedian$upperMAD = OAmedian$median + OAmedian$mad
OAmedian$lowerMAD = OAmedian$median - OAmedian$mad

ODspectra = subset(spectra,Type=='od')
ODmedian = data.frame(ppm = seq(-2.7, 7.1, length.out=512), median = sapply(ODspectra[, 2:513], median))
ODmedian$mad = sapply(ODspectra[, 2:513], mad)
ODmedian$upperMAD = ODmedian$median + ODmedian$mad
ODmedian$lowerMAD = ODmedian$median - ODmedian$mad

PIspectra = subset(spectra,Type=='pi')
PImedian = data.frame(ppm = seq(-2.7, 7.1, length.out=512), median = sapply(PIspectra[, 2:513], median))
PImedian$mad = sapply(PIspectra[, 2:513], mad)
PImedian$upperMAD = PImedian$median + PImedian$mad
PImedian$lowerMAD = PImedian$median - PImedian$mad

PNspectra = subset(spectra,Type=='pn')
PNmedian = data.frame(ppm = seq(-2.7, 7.1, length.out=512), median = sapply(PNspectra[, 2:513], median))
PNmedian$mad = sapply(PNspectra[, 2:513], mad)
PNmedian$upperMAD = PNmedian$median + PNmedian$mad
PNmedian$lowerMAD = PNmedian$median - PNmedian$mad

RAspectra = subset(spectra,Type=='ra')
RAmedian = data.frame(ppm = seq(-2.7, 7.1, length.out=512), median = sapply(RAspectra[, 2:513], median))
RAmedian$mad = sapply(RAspectra[, 2:513], mad)
RAmedian$upperMAD = RAmedian$median + RAmedian$mad
RAmedian$lowerMAD = RAmedian$median - RAmedian$mad

SCspectra = subset(spectra,Type=='sc')
SCmedian = data.frame(ppm = seq(-2.7, 7.1, length.out=512), median = sapply(SCspectra[, 2:513], median))
SCmedian$mad = sapply(SCspectra[, 2:513], mad)
SCmedian$upperMAD = SCmedian$median + SCmedian$mad
SCmedian$lowerMAD = SCmedian$median - SCmedian$mad

##########

# plot remaining 11
ggplot(A2median, aes(x=ppm, y=median)) + geom_line(size=1, alpha=0.8) + geom_ribbon(aes(ymin=A2median$lowerMAD, ymax=A2median$upperMAD), fill="black", alpha="0.3") + xlab("frequency (ppm)") + ylab("signal (arbituary unit)") + scale_x_reverse(breaks=seq(-2.5, 7.0, 0.5)) + scale_y_continuous(limits=c(-6, 35), position="right") + theme_classic(base_size=16) + ggtitle("grade II astrocytoma") + theme(plot.title=element_text(hjust=0.5, size=30), panel.grid.major=element_line(), panel.border=element_rect(fill=NA, size=2))
ggplot(A3median, aes(x=ppm, y=median)) + geom_line(size=1, alpha=0.8) + geom_ribbon(aes(ymin=A3median$lowerMAD, ymax=A3median$upperMAD), fill="black", alpha="0.3") + xlab("frequency (ppm)") + ylab("signal (arbituary unit)") + scale_x_reverse(breaks=seq(-2.5, 7.0, 0.5)) + scale_y_continuous(limits=c(-6, 35), position="right") + theme_classic(base_size=16) + ggtitle("grade III astrocytoma") + theme(plot.title=element_text(hjust=0.5, size=30), panel.grid.major=element_line(), panel.border=element_rect(fill=NA, size=2))
ggplot(ABmedian, aes(x=ppm, y=median)) + geom_line(size=1, alpha=0.8) + geom_ribbon(aes(ymin=ABmedian$lowerMAD, ymax=ABmedian$upperMAD), fill="black", alpha="0.3") + xlab("frequency (ppm)") + ylab("signal (arbituary unit)") + scale_x_reverse(breaks=seq(-2.5, 7.0, 0.5)) + scale_y_continuous(limits=c(-6, 35), position="right") + theme_classic(base_size=16) + ggtitle("abscess") + theme(plot.title=element_text(hjust=0.5, size=30), panel.grid.major=element_line(), panel.border=element_rect(fill=NA, size=2))
ggplot(HBmedian, aes(x=ppm, y=median)) + geom_line(size=1, alpha=0.8) + geom_ribbon(aes(ymin=HBmedian$lowerMAD, ymax=HBmedian$upperMAD), fill="black", alpha="0.3") + xlab("frequency (ppm)") + ylab("signal (arbituary unit)") + scale_x_reverse(breaks=seq(-2.5, 7.0, 0.5)) + scale_y_continuous(limits=c(-6, 35), position="right") + theme_classic(base_size=16) + ggtitle("haemangioblastoma") + theme(plot.title=element_text(hjust=0.5, size=30), panel.grid.major=element_line(), panel.border=element_rect(fill=NA, size=2))
ggplot(LYmedian, aes(x=ppm, y=median)) + geom_line(size=1, alpha=0.8) + geom_ribbon(aes(ymin=LYmedian$lowerMAD, ymax=LYmedian$upperMAD), fill="black", alpha="0.3") + xlab("frequency (ppm)") + ylab("signal (arbituary unit)") + scale_x_reverse(breaks=seq(-2.5, 7.0, 0.5)) + scale_y_continuous(limits=c(-6, 35), position="right") + theme_classic(base_size=16) + ggtitle("lymphoma") + theme(plot.title=element_text(hjust=0.5, size=30), panel.grid.major=element_line(), panel.border=element_rect(fill=NA, size=2))
ggplot(OAmedian, aes(x=ppm, y=median)) + geom_line(size=1, alpha=0.8) + geom_ribbon(aes(ymin=OAmedian$lowerMAD, ymax=OAmedian$upperMAD), fill="black", alpha="0.3") + xlab("frequency (ppm)") + ylab("signal (arbituary unit)") + scale_x_reverse(breaks=seq(-2.5, 7.0, 0.5)) + scale_y_continuous(limits=c(-6, 35), position="right") + theme_classic(base_size=16) + ggtitle("oligoastrocytoma") + theme(plot.title=element_text(hjust=0.5, size=30), panel.grid.major=element_line(), panel.border=element_rect(fill=NA, size=2))
ggplot(ODmedian, aes(x=ppm, y=median)) + geom_line(size=1, alpha=0.8) + geom_ribbon(aes(ymin=ODmedian$lowerMAD, ymax=ODmedian$upperMAD), fill="black", alpha="0.3") + xlab("frequency (ppm)") + ylab("signal (arbituary unit)") + scale_x_reverse(breaks=seq(-2.5, 7.0, 0.5)) + scale_y_continuous(limits=c(-6, 35), position="right") + theme_classic(base_size=16) + ggtitle("oligodendroglioma") + theme(plot.title=element_text(hjust=0.5, size=30), panel.grid.major=element_line(), panel.border=element_rect(fill=NA, size=2))
ggplot(PImedian, aes(x=ppm, y=median)) + geom_line(size=1, alpha=0.8) + geom_ribbon(aes(ymin=PImedian$lowerMAD, ymax=PImedian$upperMAD), fill="black", alpha="0.3") + xlab("frequency (ppm)") + ylab("signal (arbituary unit)") + scale_x_reverse(breaks=seq(-2.5, 7.0, 0.5)) + scale_y_continuous(limits=c(-6, 35), position="right") + theme_classic(base_size=16) + ggtitle("pituitary tumour") + theme(plot.title=element_text(hjust=0.5, size=30), panel.grid.major=element_line(), panel.border=element_rect(fill=NA, size=2))
ggplot(PNmedian, aes(x=ppm, y=median)) + geom_line(size=1, alpha=0.8) + geom_ribbon(aes(ymin=PNmedian$lowerMAD, ymax=PNmedian$upperMAD), fill="black", alpha="0.3") + xlab("frequency (ppm)") + ylab("signal (arbituary unit)") + scale_x_reverse(breaks=seq(-2.5, 7.0, 0.5)) + scale_y_continuous(limits=c(-6, 35), position="right") + theme_classic(base_size=16) + ggtitle("primitive neuroectodermal tumour") + theme(plot.title=element_text(hjust=0.5, size=30), panel.grid.major=element_line(), panel.border=element_rect(fill=NA, size=2))
ggplot(RAmedian, aes(x=ppm, y=median)) + geom_line(size=1, alpha=0.8) + geom_ribbon(aes(ymin=RAmedian$lowerMAD, ymax=RAmedian$upperMAD), fill="black", alpha="0.3") + xlab("frequency (ppm)") + ylab("signal (arbituary unit)") + scale_x_reverse(breaks=seq(-2.5, 7.0, 0.5)) + scale_y_continuous(limits=c(-6, 35), position="right") + theme_classic(base_size=16) + ggtitle("rare type tumour") + theme(plot.title=element_text(hjust=0.5, size=30), panel.grid.major=element_line(), panel.border=element_rect(fill=NA, size=2))
ggplot(SCmedian, aes(x=ppm, y=median)) + geom_line(size=1, alpha=0.8) + geom_ribbon(aes(ymin=SCmedian$lowerMAD, ymax=SCmedian$upperMAD), fill="black", alpha="0.3") + xlab("frequency (ppm)") + ylab("signal (arbituary unit)") + scale_x_reverse(breaks=seq(-2.5, 7.0, 0.5)) + scale_y_continuous(limits=c(-6, 35), position="right") + theme_classic(base_size=16) + ggtitle("schwannoma") + theme(plot.title=element_text(hjust=0.5, size=30), panel.grid.major=element_line(), panel.border=element_rect(fill=NA, size=2))

