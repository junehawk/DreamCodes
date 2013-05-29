library(graphite)
library(clipper)
library(doMC)

classes = c(rep(2,3), rep(1,8))
len = length(kegg)
filename = c("Aclacinomycin A_24_IC20_gene_based.txt", "Blebbistatin_24_IC20_gene_based.txt", "Camptothecin_24_IC20_gene_based.txt", "Cycloheximide_24_IC20_gene_based.txt", "Doxorubicin hydrochloride_24_IC20_gene_based.txt", "Etoposide_24_IC20_gene_based.txt", "Geldanamycin_24_IC20_gene_based.txt", "H-7, Dihydrochloride_24_IC20_gene_based.txt", "Methotrexate_24_IC20_gene_based.txt", "Mitomycin C_24_IC20_gene_based.txt", "Monastrol_24_IC20_gene_based.txt", "Rapamycin_24_IC20_gene_based.txt", "Trichostatin A_24_IC20_gene_based.txt", "Vincristine_24_IC20_gene_based.txt")
for(i in 1:length(filename)) {
	cat("Analysis start for drug ", i, "\n")

	exp = read.table(paste("/Users/Junehawk/Downloads/gene_based/", filename[i], sep=""), skip=1, stringsAsFactors=F, row.names=1)
	
	alphaMeanResult <- foreach(j=1:len) %dopar% {
		graph = NULL
		a = list()
		a$alphaVar = 1
		a$alphaMean = 1
		try(graph <- convertIdentifiers(kegg[[j]], "symbol"), silent=T)
		if(is.null(graph)) {
			a
		}
		else {
			graph <- pathwayGraph(graph)
			genes <- nodes(graph)
			genes <- intersect(genes, row.names(exp))
			graph <- subGraph(genes, graph)
			exps <- as.matrix(exp[genes,,drop=FALSE])
			pathwayAnalysis = NULL
			cat("pathway ", j, "\n")
			try(pathwayAnalysis <- pathQ(exps, classes, graph, nperm=100, alphaV=0.05, b=100), silent=T)
			if(is.null(pathwayAnalysis)) {
				a
			}
			else {
				pathwayAnalysis
			}
		}
	}
	cat("alpha calculation done, clipper start\n")
	clippedResult <- foreach(j=1:len) %dopar% {
		if(alphaMeanResult[[j]]$alphaVar==1)
			"NOT AVAILABLE"
		else {
			graph = NULL
			try(graph <- convertIdentifiers(kegg[[j]], "symbol"), silent=T)
			if(is.null(graph)) {
				"NOT AVAILABLE"
			}
			else {
				graph <- pathwayGraph(graph)
				genes <- nodes(graph)
				genes <- intersect(genes, row.names(exp))
				graph <- subGraph(genes, graph)
				exps <- as.matrix(exp[genes,,drop=FALSE])
				pathwayAnalysis = NULL
				try(pathwayAnalysis <- pathQ(exps, classes, graph, nperm=100, alphaV=0.05, b=100), silent=T)
				if(is.null(pathwayAnalysis)) {
					"NOT AVAILABLE"
				}
				else {
					clipped <- clipper(exps, classes, graph, "var", trZero=0.01)
					clipped <- prunePaths(clipped, thr=0.2)
					clipped
				}
			}		
		}
	}
	cat("clipping done. Saving results\n")
	write_addr = paste("Dream7/clipperOut/Drug",i,"_clipped.RData", sep="")
	save(clippedResult, file=write_addr)
	write_addr = paste("Dream7/clipperOut/Drug",i,"_alphaMean.RData", sep="")
	save(alphaMeanResult, file=write_addr)
}