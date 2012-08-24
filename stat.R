num_drug = 31

for (i in 1:num_drug) {
	if(i!=26) {
		file = paste("/Users/junehawk/Downloads/DREAM_RPPA_REG/drug",as.character(i),"_RPPA_reg.txt",sep="")
		data = read.table(file, header=TRUE, stringsAsFactors=F)
		a = subset(data, corrected<=0.05)
		cat(dim(a)[1],", ")
	}
	else
		cat(0, ", ")
}