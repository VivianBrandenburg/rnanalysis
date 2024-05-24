# @author: vivian brandenburg  

# This is a GO enrichment using topGO (https://doi.org/doi:10.18129/B9.bioc.topGO). 



##### libraries #####
library(topGO)
library(dplyr)



##### files #####
GOuniverse <- "example_data/GOterms/GOuniverse.txt"

GOterms_translator <- "example_data/GOterms/Terms_MF.txt"
ontology <- 'MF' # put 'BP' / 'CC' for other GOenrichment category

path_to_input <- 'example_data/DESeq_tables_forEnrichment/'
path_to_output <- 'example_data/GOenrichments/'


##### prepare general data #####
# Gene universe file
exp_data= read.table(GOuniverse, header=FALSE, sep='')
bg_genes=as.character(exp_data[,1])

# read in GO terms 
geneID2GO <- readMappings(GOterms_translator)
geneNames <- names(geneID2GO)

# loop over all input files
filelist <- list.files(path=path_to_input,  pattern='.csv', all.files=FALSE, full.names=FALSE)




for (exp_file in filelist){

	  
	##### prepare experiment-specific data #####
	  
	# find right data
	infile <-  paste(path_to_input, exp_file, sep='')

	# Read in genes of interest
	dataOI <- read.csv(infile, header=TRUE, sep=',')
	dataOI <- filter(dataOI, abs(log2FoldChange)>1) 
		# the original dataset was already filterd for padj values.
		# If you did not do this before, you might want to add a filter here
	genesOI <- dataOI$X
	candidate_list = genesOI

	# remove any candidate genes without GO annotation
	keep = candidate_list %in% names(geneID2GO)
	keep =which(keep==TRUE)
	candidate_list=candidate_list[keep]

	# make named factor showing which genes are of interest
	geneList=factor(as.integer(bg_genes %in% candidate_list))
	names(geneList)= bg_genes





	##### topGO #####

	# make topGO data object 
	GOdata=new('topGOdata', ontology=ontology, allGenes = geneList, 
		   annot = annFUN.gene2GO, gene2GO = geneID2GO)
	
	
	# Test for significance
	weight01_fisher=runTest(GOdata, algorithm='weight01', statistic='fisher') 

	# generate a table of results
	allGO=usedGO(GOdata)
	all_res=GenTable(GOdata, weightFisher=weight01_fisher, 
			orderBy='weightFisher', topNodes=length(allGO))


	# make outfile
	outfile_name = paste(path_to_output, exp_file,  sep='') 
	write.table(all_res, file =outfile_name,  quote = F, col.names = T, row.names = F, sep='\t')

}


