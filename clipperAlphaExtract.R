library(graphite)

drugName = c("Aclacinomycin A", "Blebbistatin", "Camptothecin", "Cycloheximide", "Doxorubicin hydrochloride", "Etoposide", "Geldanamycin", "H-7, Dihydrochloride", "Methotrexate", "Mitomycin C", "Monastrol", "Rapamycin", "Trichostatin A", "Vincristine")
alphaMean = matrix(data = 1, nrow=232, ncol=14)
alphaVar = matrix(data = 1, nrow=232, ncol=14)
colnames(alphaMean) = drugName
rownames(alphaMean) = names(kegg)
colnames(alphaVar) = drugName
rownames(alphaVar) = names(kegg)
for(i in 1:14) {
	alpha_addr = paste("/Users/Junehawk/Dream7/clipperOut/Drug",i,"_alpha.RData", sep="")
	load(alpha_addr)
	for(j in 1:length(kegg)) {
		alphaMean[j,i] = alphaMeanResult[[j]]$alphaMean
		alphaVar[j,i] = alphaMeanResult[[j]]$alphaVar
	}
}
write_addr = "/Users/Junehawk/Dream7/clipperOut/meanResult.txt"
write_addr2 = "/Users/Junehawk/Dream7/clipperOut/varResult.txt"
write.table(alphaMean,file=write_addr,sep="\t", quote=F)
write.table(alphaVar,file=write_addr2,sep="\t", quote=F)
