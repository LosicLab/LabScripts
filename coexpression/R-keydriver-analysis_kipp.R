# the KDA script, but adjusted a few ways.  
# *_keydriver.xls now outputs all nodes, with an FDR value.  
# 

library( KDA )

args <- commandArgs(TRUE);
fcausalnet <- args[1]
finputlist <- args[2]
directed <- TRUE
layer <- 6
minDsCut <- -1  #not used?
fgeneinfo <- NULL
Pcutoff <- 0.05
Bonf <- TRUE


keydriverInSubnetwork<-
function (linkpairs, signature, background = NULL, directed = T, 
    nlayers = 6, enrichedNodes_percent_cut = -1, FET_pvalue_cut = 0.05, 
    boost_hubs = T, dynamic_search = T, bonferroni_correction = T, 
    expanded_network_as_signature = F) 
{
    allnodes = sort(union(linkpairs[, 1], linkpairs[, 2]))
    no.subnetsize = length(allnodes)
    if (no.subnetsize <= 4) {
        return(NULL)
    }
    network_as_signature = length(setdiff(allnodes, signature)) == 
        0
    network_as_signature = expanded_network_as_signature | network_as_signature
    overlapped = intersect(allnodes, signature)
    no.overlapped = length(overlapped)
    if (is.null(background)) {
        background2 = c(no.subnetsize, no.overlapped)
    }
    else {
        background2 = background
    }
    keydrivers = NULL
    kdMatrix = NULL
    kdIndex = NULL
    dsnodes_list = as.list(rep(0, no.subnetsize))
    no.dsnodes = rep(0, no.subnetsize)
    cnt = 1
    intv = as.integer(no.subnetsize/10)
    if (intv == 0) {
        intv = 1
    }
    print("find downstream genes")
    if (dynamic_search) {
        layers_4search = c(1:nlayers)
    }
    else {
        layers_4search = c(nlayers)
    }
    if (network_as_signature) {
        layers_4search = c(nlayers)
    }
    dn_matrix = matrix(1, no.subnetsize, nlayers)
    for (i in c(1:no.subnetsize)) {
        if (i%%intv == 0) {
            print(paste(i, "/", no.subnetsize))
        }
        minpv = 1
        min_nohits = 0
        min_noidn = 0
        min_layer = 0
        min_dn = 0
        min_fc = 0
        minpvW = 1
        min_fcW = 0
        for (y in layers_4search) {
            idn = downStreamGenes(netpairs = linkpairs, seednodes = allnodes[i], 
                N = y, directed = directed)
            idn = setdiff(idn, allnodes[i])
            no.idn = length(idn)
            dn_matrix[i, y] = no.idn
            if (!network_as_signature) {
                hits = intersect(idn, overlapped)
                no.hits = length(hits)
                if (no.hits == 0) {
                  next
                }
                foldchg = (no.hits/no.idn)/(no.overlapped/no.subnetsize)
                pv = phyper(no.hits - 1, no.idn, no.subnetsize - 
                  no.idn, no.overlapped, lower.tail = F)
                foldchgW = (no.hits/no.idn)/(background2[2]/background2[1])
                pvW = phyper(no.hits - 1, no.idn, background2[1] - 
                  no.idn, background2[2], lower.tail = F)
				#print(  rbind(kippdat, (cbind(i, y, foldchg, pv))))
				#kippdat<<-rbind(kippdat, (cbind(i, y, foldchg, pv)))
				
                if (pv < minpv) {
                  minpv = pv
                  min_nohits = no.hits
                  min_noidn = no.idn
                  min_layer = y
                  min_fc = foldchg
                  minpvW = pvW
                  min_fcW = foldchgW
                }
            }
            else {
                no.hits = no.idn
                minpv = 0
                min_nohits = no.idn
                min_noidn = no.idn
                min_layer = y
                min_fc = 1
            }
        }
        if (no.idn > 0) {
            dsnodes_list[[i]] = idn
            no.dsnodes[i] = no.idn
        }
        correct_minpv = minpv * no.subnetsize
        correct_minpv = ifelse(correct_minpv > 1, 1, correct_minpv)
        res = c(min_nohits, min_noidn, no.overlapped, no.subnetsize, 
            background2[2], background2[1], length(signature), 
            min_layer, min_fcW, minpvW, min_fc, minpv, correct_minpv)
        kdMatrix = rbind(kdMatrix, res)
		#kippdat<<-kdMatrix
    }
    colnames(kdMatrix) <- c("hits", "downstream", "signature_in_subnetwork", 
        "subnetwork_size", "signature_in_network", "network_size", 
        "signature", "optimal_layer", "fold_change_whole", "pvalue_whole", 
        "fold_change_subnet", "pvalue_subnet", "pvalue_corrected_subnet")
    optLidx = getMatchedIndexFast(colnames(kdMatrix), "optimal_layer")
    mymincut = enrichedNodes_percent_cut * no.overlapped
    if (enrichedNodes_percent_cut <= 0) {
        mincutLayers = apply(dn_matrix, 2, mean) + apply(dn_matrix, 
            2, sd)
        opL = ifelse(kdMatrix[, optLidx] == 0, 1, kdMatrix[, 
            optLidx])
        mymincut = mincutLayers[opL]
    }
    cutmatrix = c(mean(no.dsnodes), sd(no.dsnodes), mymincut)
    ncols = dim(kdMatrix)[2]
    if (bonferroni_correction) {
        kdSel = (kdMatrix[, ncols] < FET_pvalue_cut) & (kdMatrix[,2] >= mymincut)
    }
    else {
        kdSel = (kdMatrix[, ncols - 1] < FET_pvalue_cut) & (kdMatrix[, 
            2] >= mymincut)
    }
    keydrv = rep(0, no.subnetsize)
    if (sum(kdSel) > 0) {
        keydrivers = allnodes[kdSel]
        kdIndex = c(1:no.subnetsize)[kdSel]
        n.drivers = length(keydrivers)
        for (i in c(1:n.drivers)) {
            iselA = (kdMatrix[, 2] > kdMatrix[kdIndex[i], 2]) & 
                kdSel
            isel = c(1:no.subnetsize)[iselA]
            if (sum(isel) > 0) {
                if (directed) {
                  ilocal = setInSets(setC = allnodes[kdIndex[i]], 
                    setlist = dsnodes_list[isel])
                }
                else {
                  ilocal = setInSets(setC = dsnodes_list[[kdIndex[i]]], 
                    setlist = dsnodes_list[isel])
                }
                keydrv[kdIndex[i]] = !ilocal + 0
            }
            else {
                keydrv[kdIndex[i]] = TRUE
            }
        }
    }
    else {
        print("Warning: the downstream metric is not used as the specified minimumal downstream size is too big !!")
    }
    if (boost_hubs) {
        if (!network_as_signature) {
            kdSelB = rep(F, no.subnetsize)
            kdSelB[kdIndex] = T
            psel = kdMatrix[, ncols - 3] * no.subnetsize < 0.05
            kdPsd = 1
            kdPmean = 1
            kdpvalues = -log10(kdMatrix[, ncols - 3])
            kdpvalues = ifelse(is.na(kdpvalues), 0, kdpvalues)
            if (sum(psel) > 0) {
                kdPmean = mean(kdpvalues[psel])
                kdPsd = sd(kdpvalues[psel])
                print(as.numeric(signif(kdpvalues[psel], 2)))
                directSel = (kdpvalues > (kdPmean + kdPsd))
                directSel = ifelse(is.na(directSel), FALSE, directSel)
                if (sum(directSel) > 0) {
                  kdSel = kdSel | directSel
                  dIndex = c(1:no.subnetsize)[kdSel]
                  keydrv = rep(F, no.subnetsize)
                  keydrv[dIndex] = TRUE
                }
            }
            cutmatrix = rbind(c(mean(no.dsnodes), sd(no.dsnodes), 
                concatenate(mymincut, ";"), kdPmean, kdPsd, kdPmean + 
                  kdPsd))
            colnames(cutmatrix) <- c("mean_downstream", "sd_downstream", 
                "enrichedNodes_cut", "mean_logP", "sd_logP", 
                "cut_logP")
        }
        else {
            kdSelB = rep(TRUE, no.subnetsize)
            mydegree = degreeByLinkPairs(linkpairs = linkpairs, 
                directed = directed, cleangarbage = F)
            mIdx = getMatchedIndexFast(rownames(mydegree), allnodes)
            mydegree = mydegree[mIdx, ]
            if (directed) {
                directSel = mydegree[, 2] > mean(mydegree[, 2]) + 
                  2 * sd(mydegree[, 2])
                cutmatrix = rbind(c(mean(no.dsnodes), sd(no.dsnodes), 
                  concatenate(mymincut, ";"), mean(mydegree[, 
                    2]), sd(mydegree[, 2]), mean(mydegree[, 2]) + 
                    2 * sd(mydegree[, 2])))
            }
            else {
                directSel = mydegree[, 3] > mean(mydegree[, 3]) + 
                  2 * sd(mydegree[, 3])
                cutmatrix = rbind(c(mean(no.dsnodes), sd(no.dsnodes), 
                  concatenate(mymincut, ";"), mean(mydegree[, 
                    3]), sd(mydegree[, 3]), mean(mydegree[, 3]) + 
                    2 * sd(mydegree[, 3])))
            }
            directSel = directSel & kdSelB
            directeHub = rownames(mydegree)[directSel]
            isDirectHub = setElementInSet(allnodes, directeHub)
            keydrv[isDirectHub] = T
            kdSel = kdSel | isDirectHub
            colnames(cutmatrix) <- c("mean_downstream", "sd_downstream", 
                "cut_downstream", "mean_degree", "sd_degree", 
                "cut_degree")
        }
    }
    else {
        cutmatrix = rbind(c(mean(no.dsnodes), sd(no.dsnodes), 
            concatenate(mymincut, ";"), "F"))
        colnames(cutmatrix) <- c("mean_downstream", "sd_downstream", 
            "cut_downstream", "boost_directhubs")
    }
    if (sum(kdSel) == 0) {
        return(NULL)
    }
    is_signature = rep(0, no.subnetsize)
    names(is_signature) <- as.character(allnodes)
    is_signature[as.character(overlapped)] = 1
    fkd = cbind(allnodes, is_signature, kdMatrix, keydrv + 0)[kdSel,]
	kippdat<<-cbind(allnodes, is_signature, kdMatrix, keydrv + 0)
    if (sum(kdSel) > 1) {
        nf.cols = dim(fkd)[2]
        if (network_as_signature) {
            mo = order(-as.integer(fkd[, 3]))
        }
        else {
            mo = order(as.numeric(fkd[, nf.cols - 1]))
        }
        fkd = fkd[mo, ]
        mo = order(-as.integer(fkd[, nf.cols]))
        fkd = fkd[mo, ]
    }
    else {
        fkd = rbind(fkd)
    }
    colnames(fkd) <- c("keydrivers", "is_signature", "hits", 
        "downstream", "signature_in_subnetwork", "subnetwork_size", 
        "signature_in_network", "network_size", "signature", 
        "optimal_layer", "fold_change_whole", "pvalue_whole", 
        "fold_change_subnet", "pvalue_subnet", "pvalue_corrected_subnet", 
        "keydriver")
	colnames(kippdat)<<-colnames(fkd)
	kippdat<<-as.data.frame(kippdat)
	kippdat$FDR<<-p.adjust(as.numeric(as.character(kippdat$pvalue_whole)), method="BH")
    
	print(fkd)
    return(list(fkd, cutmatrix))
}


###################################################################################################
#
#  This file is the main program for Key Driver Analysis (KDA)
#
#  Function: to identify the key regulators for a list of genes with regard to a given network.
#
#  Author: Bin Zhang, PhD, Mount Sinai School of Medicine, new York, USA
#  Contact: bin.zhang@mssm.edu
#
#  Release Date: Oct 30, 2012
#  
#  Modified by ______________  Date ____________
#
#
#    1. Jun Zhu, Bin Zhang, Erin Smith, Becky Drees, Rachel Brem, Roger Bumgarner, Eric E. Schadt. 
#       (2008) Integrating Large-Scale Functional Genomic Data to Dissect the Complexity of Yeast 
#       Regulatory Networks. Nature Genetics 40: 854-861
#    2. I-Ming Wang*§, Bin Zhang*§, Xia Yang*§ et al. (2012) Systems Analysis of Eleven Rodent Disease 
#       Models Reveals an Inflammatome Signature and Key Drivers. Molecular Systems Biology 8:594
#    3. Linh Tran*, Bin Zhang*§, et al. (2011) Inferring Causal Genomic Alterations in Breast Cancer 
#       Using Gene Expression Data. BMC Systems Biology 5(121)
#    4. Xia Yang*, Bin Zhang*, et al. (2010) Genetic and Genomic Analysis of Cytochrome P450 Enzyme 
#       Activities in Human Liver. Genome Research 20: 1020-1036.
#    5. Bin Zhang, Linh M. Tran, Jonathan Derry, Jun Zhu, Stephen Friend. Breast Cancer Transcriptional
#       Networks: Pathways, Regulators, and Prediction. Manuscript
#

##----------------------------------------------------------------------------------------------------
# Input: 
#   1. gene list: finputlist ="input_gene_lists.xls"
#   2. causal network: fcausalnet="input_breastcancer_BN.xls"; 
#   3. directed or undirected network: directed = T
#   4. expand subnetwork based on L-layer neighbors: layer=0
#   5. gene annotation file(NULL if not available): fgeneinfo
#
# Output:
#    1. keydrivers for each subnetwork: "*_keydriver.xls"
#    2. Cytoscape network: *_cys.txt
#    3. Cytoscape node properties: *_cys-nodes.txt
#    4. combined keydrivers for all runs: "_KDx_combined.xls"
#    5. parameters used for key driver analysis: "_KDx_parameters.xls"
#
#
# ---------------------------------- Parameters to be changed -----------------------------------------

# example 1: key drivers for breast cancer gene modules
#
#data( breastcausalnet )
#data( breastinputlist )
#fcausalnet <- breastcausalnet
#finputlist <- breastinputlist
#directed <- TRUE
#layer <- 0
#minDsCut <- -1
#fgeneinfo <- NULL
#
## example 2: key drivers for yeast eQTL hotspots
##
#data( yeastcausalnet )
#data( yeastinputlist )
#fcausalnet <- yeastcausalnet
#finputlist <- yeastinputlist
#fgeneinfo <- "yeast-geneinfo.xls"
#directed <- TRUE
#layer <- 1
#minDsCut <- 5


# 2. specify the directory for holding analysis results
#
if ( directed )
{
  outputDir <- "KeyDrivers/"
} else
{
  outputDir <- "KeyDrivers-undirected/"
}

dir.create( outputDir )

#
# -----------------------------End of Parameters to be changed --------------------------------------

library( class )
library( cluster )
library( rpart )
#library( sma ) # this is needed for plot.mat below
library( lattice ) # require is design for use inside functions 
#memory.size( TRUE )   # check the maximum memory that can be allocated
#memory.limit( size = 3800 )   # increase the available memory

################################################################################################
#    1. read in network

cnet <- read.delim( fcausalnet , sep = "\t" , header = TRUE )
cnet <- as.matrix( cnet )
dim( cnet )

totalnodes <- union( cnet[,1] , cnet[,2] )

fname <- getFileName( fcausalnet )
fname <- paste( fname , "_L" , layer , sep = "" )

################################################################################################
# 2. read in gene lists

listMatrix <- read.delim( finputlist , sep="\t" , header = TRUE )
dim( listMatrix )
listMatrix <- as.matrix( listMatrix )
listMatrix[1:2,]
ncols <- dim( listMatrix )[2]

modules <- names( table( listMatrix[,ncols] ) )

xkdFall <- paste( outputDir , fname , "_KDx_combined.xls" , sep = "" )
xkdFpara <- paste( outputDir , fname , "_KDx_parameters.xls" , sep = "" )
xkdrMatrix <- NULL
paraMatrix <- NULL

################################################################################################
# 3. process each gene list
#

for ( em in modules )
{

#em="green"

  print( paste( "*****************" , em , "********************" ) )

  esel <- listMatrix[,ncols] == em

# remove abnormal gene names
#
  genes <- union( listMatrix[esel,1] , NULL )
  genes <- genes[genes != ""]
  genes <- genes[!is.na( genes )]
  no.genes <- length( genes )

  em2 <- replaceString( em , ":" , "" )
  em2 <- replaceString( em2 , " " , "-" )

  key2 <- paste( fname , "_KD_" , em2 , sep = "" )
  onetFname <- paste( outputDir , key2 , ".pair" , sep = "" )
  snpFname <- paste( outputDir , key2 , ".snp" , sep = "" )
  kdFname <- paste( outputDir , key2 , "_keydriver.xls" , sep = "" )

  if(layer >=1 )
  {
     # expand network by K-hop nearest neighbors layers
     expandNet <- findNLayerNeighborsLinkPairs( linkpairs = cnet , subnetNodes = genes ,
			                                     nlayers = layer , directed = FALSE )
  }  else
  {
   # no expansion
     expandNet <- getSubnetworkLinkPairs( linkpairs = cnet , subnetNodes = genes )
  }
  dim( expandNet )

  allnodes <- union( expandNet[,1] , expandNet[,2] )
#write.table(expandNet, onetFname, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)

################################################################################################
# 4. keydriver for a given network
#

  if (directed)
  {
    ret <- keydriverInSubnetwork( linkpairs = expandNet , signature = genes, background=NULL, directed = directed ,
			         nlayers = 6 , enrichedNodes_percent_cut=-1, FET_pvalue_cut=Pcutoff,
			         #nlayers = 6 , enrichedNodes_percent_cut=-1, FET_pvalue_cut=0.05,
			         boost_hubs=T, dynamic_search=T, bonferroni_correction=Bonf, expanded_network_as_signature =F)
			         #boost_hubs=T, dynamic_search=T, bonferroni_correction=T, expanded_network_as_signature =F)
  }  else
  {
    ret <- keydriverInSubnetwork( linkpairs = expandNet , signature = genes , directed = directed ,
			         nlayers = 2 , enrichedNodes_percent_cut=-1, FET_pvalue_cut=0.05,
				boost_hubs=T, dynamic_search=T, bonferroni_correction=T, expanded_network_as_signature =F)
  }

  if ( is.null( ret ) )
  {
	next
  }

  fkd <- ret[[1]]
  parameters <- ret[[2]]

  fkd2 <- cbind( rep( em , dim( fkd )[1] ) , fkd )
  xkdrMatrix <- rbind( xkdrMatrix , fkd2 )

  paraMatrix <- rbind( paraMatrix , c( key2 , parameters ) )
  
  write.table( kippdat , kdFname , sep = "\t" , quote = FALSE , col.names = TRUE , row.names = FALSE )

################################################################################################
# 4. output networks & key drivers for visualization
#
#     Cytoscape output: 1) network file - *_cys.txt 2) node property file: *_cys-nodes.txt
#

      nodeprop = configureNodeVisualization(allnodes=allnodes, signature=genes, kdaMatrix=fkd)

      hnList     = nodeprop[[1]] # node subcategpries
      listprop   = nodeprop[[2]] # visual properties for each subcategory
      legend     = nodeprop[[3]] # legend table for visual propertie

      resf = makeSNP(netpairsWtype   = expandNet, 
               edgecolorlevels = c("grey"),
               highlightNodes  = hnList,
               normColor="grey",   highColor=listprop[,1],
               normShape="circle", highShape=listprop[,2],
               normNodeSize ="40",  highNodeSize =listprop[,3],
               normFontSize ="12",  highFontSize =listprop[,4],
               legendtable=legend, snafile=snpFname )

} #for (em in modules) {


# save all key drivers
colnames( xkdrMatrix ) <- c( "module" , colnames( fkd ) )
write.table( xkdrMatrix , xkdFall , sep = "\t" , quote = FALSE , col.names = TRUE , 
		     row.names = FALSE )     

# save parameters used
#
colnames( paraMatrix ) <- c( "subnet" , colnames( parameters ) )
write.table( paraMatrix , xkdFpara , sep = "\t" , quote = FALSE , col.names = TRUE ,
		     row.names = FALSE )

if( !is.null( fgeneinfo ) )
{
   infoMatrix <- read.delim( fgeneinfo , sep = "\t" , header = TRUE )
   dim( infoMatrix )
   xkdrMatrix2 <- cbind( xkdrMatrix , c( 1:( dim( xkdrMatrix )[1] ) ) )

   ic <- dim( infoMatrix )[2] + 1
   merged <- merge( infoMatrix , xkdrMatrix2 , by.x = 1 , by.y = 2 , all.y = TRUE )
   merged <- as.matrix( merged )
   ic2 <- dim( merged )[2]
   xf <- merged[,c( ic , setdiff( 1:ic2 , ic ) )]
   ic3 <- dim( merged )[2]
   mo <- order( as.integer( merged[,ic3] ) )
   write.table( xf[mo,-ic3] , xkdFall , sep = "\t" , quote = FALSE , col.names = TRUE ,
		        row.names = FALSE )
}

## ------------------------------------- END ------------------------------------------------
