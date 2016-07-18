# CALL (example):
# source('robust.clustering.metagenomics.functions.r')

# With phyloseq object in a R data file:
# robust.clustering.all.steps("~/RobustClustering/David2014","~/RobustClustering/David2014/data.norm_David2014.RData",'David2014',"COLLECTION_DAY")

# Or with .biom file:
# robust.clustering.all.steps("~/RobustClustering/David2014","~/RobustClustering/David2014/David2014.biom",'David2014',"COLLECTION_DAY","~/RobustClustering/David2014/mapping_David2014.tsv")


#####################################
### Function write.silhouette.output
#####################################
# Writes a clustering assessment in a file, based on Silhouette values.
#
# Args:
#   ssi: numerical summary of the individual silhouette widths s(i).
#        Object of class summary.silhouette.
#   si: object of class silhouette, returned by silhouette() function.
#   fileName: file name where writting the results.
write.silhouette.output <- function(ssi=NULL,si=NULL,fileName=NULL){
  sink(fileName)
  #if(is.na(si)){
  #  print('Only 1 cluster!!!!')  
  #}else{
    print(ssi)
    print('Sd. among clusters: ')
    print(sd(ssi$clus.avg.widths))
    print('Sd. among samples: ')
    print(sd(si[,"sil_width"]))
  #} # end-if
  sink()
} #end-function write.silhouette.output



################################################################
### Function plot.robust.clustering.all.together.formatted   ###
################################################################
# Plots different summary graphs, with the results of the different clustering procedures,
# combined in one output file (robustClustering_allTogether_formatted.pdf).
#
# Args:
#   eval.array4d: 4-dimensions array with clustering results, returned by robust.clustering() function.
plot.robust.clustering.all.together.formatted <- function(eval.array4d){
  library(grid)  
  library(gridExtra)  
  #library(cowplot)  
  
  plots<- list()
  # To join in a data.frame the selected metrics for each graph:
  for (alg in c('pam','hclust')){
    for (score in c('SI','PS','Jaccard')){
      assFrame <- NULL
      assFrame <- data.frame(k=2:10, JSD=eval.array4d[,score,'jsd',alg], rJSD=eval.array4d[,score,'rjsd',alg], BrayCurtis=eval.array4d[,score,'bray',alg], MorisitaHorn=eval.array4d[,score,'horn',alg], Kulczynski=eval.array4d[,score,'kulczynski',alg])
      assFrame <- melt(assFrame, id.vars=c("k"), variable.name="metric", value.name="score")
      gp <-NULL
      if(score=='SI'){
        gp <- (ggplot(assFrame, aes(x=k, y=score, color=metric)) +
                 geom_point(aes(shape=metric),size=3) + geom_line(aes(linetype=metric),size=0.75) +
                 theme_bw() + # background grey to white, and re-initialize all associated parameters
                 scale_x_continuous(breaks=2:10, labels=2:10) +
                 coord_cartesian(ylim = c(0,1)) +
                 geom_hline(yintercept = 0.25, linetype=2, size=1) + 
                 geom_hline(yintercept = 0.50, linetype=2, size=1) + 
                 ggtitle(paste(toupper(alg),score,sep=' - ')) +
                 theme(text = element_text(size = 20)) # increase all font size
               # + geom_errorbar(aes(x = 3, y = 0.9, ymin = 0.9 - sd, ymax = 0.9 + sd), colour = 'red', width = 0.4)
        ) #end-gplot
      }else if(score=="PS"){
        gp <- (ggplot(assFrame, aes(x=k, y=score, color=metric)) +
                 theme_bw() +
                 geom_point(aes(shape=metric),size=3) + geom_line(aes(linetype=metric),size=0.75) +
                 scale_x_continuous(breaks=2:10, labels=2:10) +
                 coord_cartesian(ylim = c(0,1)) +
                 geom_hline(yintercept = 0.8, linetype=2, size=1) + 
                 ggtitle(paste(toupper(alg),score,sep=' - ')) +
                 theme(text = element_text(size = 20))
        ) #end-gplot
      }else{
        gp <- (ggplot(assFrame, aes(x=k, y=score, color=metric)) +
                 theme_bw() +
                 geom_point(aes(shape=metric),size=3) + geom_line(aes(linetype=metric),size=0.75) +
                 scale_x_continuous(breaks=2:10, labels=2:10) +
                 coord_cartesian(ylim = c(0,1)) +
                 geom_hline(yintercept = 0.75, linetype=2, size=1) + 
                 ggtitle(paste(toupper(alg),score,sep=' - ')) +
                 theme(text = element_text(size = 20))
        ) #end-gplot
      } # end-if
      plots=c(plots,list(gp))
    } #end-for-score
  } #end-for-alg
  
  # Print one combined legend at the bottom.
  pdf("robustClustering_allTogether_formatted.pdf",width=12,height=8)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl=c(gl,list(nrow=2))
  combined<-arrangeGrob(do.call(arrangeGrob, gl),
                        legend,
                        ncol = 1,
                        nrow = 2,
                        heights = unit.c(unit(1, "npc") - lheight, lheight))
  grid.draw(combined)
  dev.off()
  #graphics.off()
} #end-function-plot.robust.clustering.all.together.formatted



####################################
### Function robust.clustering   ###
####################################
# Main function of the algorithm Robust Clustering in Metagenomics. It computes
# and collect the results of all the different clustering procedures, including:
#     2 clustering algorithms
#   x 9 (k=2:10) potential number of clusters
#   x 5 distance measures
#   x 3 assessment scores
# 
# Args:
#   data.norm: phyloseq object, with otu_table with relative and normalized abundances.
#   var.color: string, variable within phyloseq object to use in the PCoA plots. Ex:"age".
#   label: string, used as label and suffix of files and folders. Often, it is the dataset name.
#
# Returns:
#   4-dimensions array with clustering results.
robust.clustering <- function(data.norm=NULL,var.color=NULL,label=NULL){
  # Load libraries
  library(lattice,quietly=TRUE)
  library(cluster,quietly=TRUE)
  library(fpc,quietly=TRUE)
  library(reshape2,quietly=TRUE) # Load the reshape2 package (for the melt() function)
  library(ggplot2,quietly=TRUE)
  
  data.norm.backup <- data.norm
  iniSeed <- 12345
  # 2 clustering algorithms x 9 (k=2:10) potential number of clusters x 5 distance measures x 3 assessment scores
  # Assigning values: eval.array4d['2','SI','jsd','pam']=0.2365
  # Assigning values: eval.array4d[as.character(k),score,distanceMethod,alg]
  eval.array4d <- array(0,dim=c(9,3,5,2),dimnames=list(list('2','3','4','5','6','7','8','9','10'), list('SI','PS','Jaccard'), list('jsd','rjsd','bray','horn','kulczynski'), list('pam','hclust')))
  
  # a.- different distance measures: JSD, rootJSD, Bray-Curtis, Morisita-Horn (Horn in Vegan) and Kulczynski.
  for (distanceMethod in c("jsd", "rjsd", "bray", "horn", "kulczynski")){
    dir.create(distanceMethod)
    setwd(distanceMethod)
    
    # Compute distance
    if(distanceMethod == "rjsd"){
      dist <- distance(otu_table(data.norm),method="jsd")
      distMat <- as.matrix(dist)
      rootDistMat <- apply(distMat,MARGIN=c(1,2),sqrt)
      dist <- as.dist(rootDistMat)
      rm(distMat,rootDistMat)
    }else{
      dist <- distance(otu_table(data.norm),method=distanceMethod)
    }
    save(dist,file=paste('dist_',distanceMethod,'.RData',sep=''))
    # Distance plot
    pdf(paste("sampleDistribution_",ntaxa(data.norm),"taxa_ordinateMDS_distance_",distanceMethod,".pdf",sep=''),width=10)
    print(plot_ordination(data.norm, ordinate(data.norm, "MDS", distance=dist), color=var.color))
    dev.off()
    
    # Applying clustering: to fill in the eval.array4d
    for (alg in c("pam","hclust")){
      ### To compute robustness/assestment measures:
      ## A.-Metric SI: Silhouette
      score='SI'  
      switch(alg,
             pam={
               for(k in 2:10){
                 fit <- pam(dist,k,diss=TRUE)
                 eval.array4d[as.character(k),score,distanceMethod,alg] <- summary(silhouette(fit))$avg.width
               } # end-for-k
             }, # end-Pam
             
             hclust={
               fit <- hclust(dist,method="average") 
               for(k in 2:10){
                 eval.array4d[as.character(k),score,distanceMethod,alg] <- summary(silhouette(cutree(fit,k=k),dist))$avg.width
               } # end-for-k
             }, # end-Hclust
             { print('ERROR: algorithm should be "pam" or "hclust"!!') }
      ) #end-switch-alg
      
      ## B.-Metric PS: Prediction strength
      score='PS'
      set.seed(iniSeed)
      switch(alg,
             # claraCBI is the interface for pam and clara. In our case, claraCBI should apply pam, since clara is not available with distance matrix as input, and we indicate clearly distance=true.
             pam={ out.pred.str <- prediction.strength(dist, Gmin=2, Gmax=10, M=100, clustermethod=claraCBI, classification="centroid", distances=TRUE) },
             hclust={ out.pred.str <- prediction.strength(dist, Gmin=2, Gmax=10, M=100, clustermethod=disthclustCBI,method="average", classification="averagedist", distances=TRUE) }
      ) #end-switch-alg
      eval.array4d[,score,distanceMethod,alg] <- out.pred.str$mean.pred[2:10]
      
      ## C.-Metric Jaccard: Jaccard similarity after resampling
      # With 'boot' (multiple points or not?) or with 'subset'? They are the only two bootmethods working with dissimilarity matrix. --> boot is selected, as the default method, with the default options (omitting multiple points in the computation of the Jaccard index).
      score='Jaccard'
      switch(alg,
             pam={
               for(k in 2:10){
                 cf <- clusterboot(dist,B=100,distances=TRUE,bootmethod="boot",clustermethod=claraCBI,k=k,seed=iniSeed,count=FALSE)
                 eval.array4d[as.character(k),score,distanceMethod,alg] <- mean(cf$bootmean)
                 ## sd(cf$bootmean)
               } # end-for-k
             }, # end-Pam
             hclust={ 
               for(k in 2:10){
                 cf <- clusterboot(dist,B=100,distances=TRUE,bootmethod="boot",clustermethod=disthclustCBI,method='average',k=k,seed=iniSeed,count=FALSE)
                 eval.array4d[as.character(k),score,distanceMethod,alg] <- mean(cf$bootmean)
               } # end-for-k
             } # end-Hclust
      ) #end-switch-alg
      # Alternative with subset: cf <- clusterboot(dist,B=100,distances=TRUE,bootmethod="subset",subtunning=(2/3)*nsamples(data.norm),clustermethod=disthclustCBI,method="average",k=k,seed=iniSeed)

      # Computes the results with the best SI, with both algorithm,  without comparison between different distances, just as additional information in output files.
      switch(alg, # Switch computes 'k' and 'si' for different algorithms
             pam={
               fit <- pamk(dist,krange=2:10,diss=TRUE)
               k <- fit$nc
               # Once the k value is selected, to compute si with silhouette() method to be comparable with the remainder runs, to avoid the re-ordering of the sample columns, with the corresponding wrong assignments of cluster number to each sample.
               si <- silhouette(fit$pamobject$clustering,dist)
             }, # end-pam
             hclust={ 
               fit <- hclust(dist,method="average")
               my.k.choices <- 2:10
               current.avgWidth <- 0
               for (ii in (1:length(my.k.choices)) ){
                 new.avgWidth <- summary(silhouette(cutree(fit,k=my.k.choices[ii]),dist))$avg.width
                 #print(paste('k=',my.k.choices[ii],': ',new.avgWidth,sep=''))
                 if(new.avgWidth > current.avgWidth){
                   si <- silhouette(cutree(fit,k=my.k.choices[ii]),dist)
                   current.avgWidth <- new.avgWidth
                   k <- my.k.choices[ii]
                 } #end-if better SI
               } #end-for k values
             } # end-hclust
      ) #end-switch-alg
      ssi <- summary(si)
      write.silhouette.output(ssi,si,paste('output_silhouette_kbestSI_',alg,'_',distanceMethod,'.txt',sep=''))
    } # end-for clustering algorithm    
    # Go back to the root pathway before next loop iteration
    setwd("..")
  } # end-for distance measure  
  save(eval.array4d,file=paste('eval.arrays_',label,'.RData',sep=''))
  
  rm(data.norm)
  return("array4d" = eval.array4d)
} #end-function robust.clustering


###########################################
### Function robust.clustering.decision ###
###########################################
# Automatically decides the grouping in states/clusters. According the following
# criteria:
#   1.- to select the “k” with the  highest SI (in 2-10) in 2 measures
#   2.- to check if this “k” value pass the PS limits to be robust
#   3.- to check if those “k” clusters are stable according Jaccard threshold from bootstrapping process
#
# Args:
#  data.now: phyloseq object, with all samples and OTU information.
#  eval.array4d: 4-dimensions array with clustering results, returned by robust.clustering() function.
#  var.color: string, variable within phyloseq object to use in the PCoA plots. Ex:"age".
#  label: string, used as label and suffix of files and folders. Often, it is the dataset name.
#
# Returns:
#   Phyloseq object with the cluster assigned to each sample, in a new variable.
robust.clustering.decision <- function(data.now=NULL,eval.array4d=NULL,var.color=NULL,label=NULL){
  library(modeest,quietly=TRUE) # For compute the mode with mlv()
  
  thresh.SI <- 0.25
  thresh.PS <- 0.80
  thresh.Jac <- 0.75
  
  no.agree.methods <- 2
  
  # 1.- to select the “k” with the  highest SI (in 2-10) in 2 measures
  #kbestSI <- array(0,dim=c(1,5),dimnames=list(listlist('jsd','rjsd','bray','horn','kulczynski'))
  SI.ok <- TRUE
  kbestSI <- NULL
  for (distanceMethod in c("jsd", "rjsd", "bray", "horn", "kulczynski")){
    kbestSI[distanceMethod] <- as.integer(names(which.max(eval.array4d[,'SI',distanceMethod,'pam'])))
  } # end-for distanceMethod
  kbestSI.mode <- mlv(kbestSI)[1]
  kbestSI.mode.char <- as.character(kbestSI.mode)
  kbestSI.mode.methods <- names(which(kbestSI==kbestSI.mode))
  # To check if the "k" mode value is higher than the SI threshold in at least two measures
  kbestSI.mode.methods <- names(eval.array4d[kbestSI.mode.char,'SI',kbestSI.mode.methods,'pam'] >= thresh.SI)
  if(length(kbestSI.mode.methods) < no.agree.methods){
    print(paste('ERROR: Not shared k value greater than SI threshold (',thresh.SI,') in at least ',no.agree.methods,' distance methods!!',sep=''))
    SI.ok <- FALSE
  } # end-if k > SI.thresh in >= no.agree.methods (i.e. 2) measures
  
  # 2.- to check if this “k” value pass the PS limits to be robust (in the selected 2 measures)
  PS.ok <- TRUE
  if(SI.ok){
    # To check if the "k" mode value is higher than the PS threshold in at least two measures
    kbestSI.mode.methods <- names(eval.array4d[kbestSI.mode.char,'PS',kbestSI.mode.methods,'pam'] >= thresh.PS)
    if(length(kbestSI.mode.methods) < no.agree.methods){
      print(paste('ERROR: Not shared k value greater than PS threshold (',thresh.PS,') in at least ',no.agree.methods,' distance methods!!',sep=''))
      PS.ok <- FALSE
    } # end-if k > Ps.thresh in >= no.agree.methods (i.e. 2) measures
  } # end-if SI.ok
  
  # 3.- to check if those “k” clusters are stable according Jaccard threshold from bootstrapping process
  Jac.ok <- TRUE
  if((SI.ok)&&(PS.ok)){
    kbestSI.mode.methods <- names(eval.array4d[kbestSI.mode.char,'Jaccard',kbestSI.mode.methods,'pam'] >= thresh.Jac)
    if(length(kbestSI.mode.methods) < no.agree.methods){
      print(paste('ERROR: Not shared k value greater than Jaccard threshold (',thresh.Jac,') in at least ',no.agree.methods,' distance methods!!',sep=''))
      Jac.ok <- FALSE
    } # end-if k > Ps.thresh in >= no.agree.methods (i.e. 2) measures
  } # end-if SI.ok
  
  # 4.- When all the criterion are satisfied, the distance measure with the highest SI is selected to apply and return the clustering results in a phyloseq object.
  if((SI.ok)&&(PS.ok)&&(Jac.ok)){
    distanceMethod <- names(which.max(eval.array4d[kbestSI.mode.char,'SI',kbestSI.mode.methods,'pam']))
    
    # Print definitive numerical results in a text file
    sink(paste('definitiveClusteringResults_',label,'.txt',sep=''))
    cat(paste('No.clusters: ',kbestSI.mode,'\n',sep=''))
    cat(paste('Distance Method: ',distanceMethod,'\n',sep=''))
    cat('Clustering algorithm: PAM\n')
    cat(paste('Assessment score: Average width silhouette: ',eval.array4d[kbestSI.mode.char,'SI',distanceMethod,'pam'],'\n',sep=''))
    cat(paste('Assessment score: Prediction strength: ',eval.array4d[kbestSI.mode.char,'PS',distanceMethod,'pam'],'\n',sep=''))
    cat(paste('Assessment score: Jaccard similarity: ',eval.array4d[kbestSI.mode.char,'Jaccard',distanceMethod,'pam'],'\n',sep=''))
    cat('--------------------------\n\n')
    
    # To load pre-computed distances rather than to recompute them
    #dist <- distance(otu_table(data.now),method=distanceMethod)
    load(paste(distanceMethod,'/dist_',distanceMethod,'.RData',sep=''))
    fit <- pam(dist,kbestSI.mode,diss=TRUE)
    cat('SILHOUETTE RESULTS:\n')
    # IMP: Not silhouette(fit,dist), because the order of the instance is WRONG, and the assignation of cluster_id is WRONG in si[,1], which is assigned to sample_data(data.now)$cluster!!!
    # Thus, row.names(sample_data(data.now)[1:10]) == labels(dist)[1:10] == labels(fitPam$clustering[1:10]) != labels(fitPam[1:10])
    si <- silhouette(fit$clustering,dist) 
    ssi <- summary(si)
    print(ssi)
    cat('Sd. among clusters: ')
    print(sd(ssi$clus.avg.widths))
    cat('Sd. among samples: ')
    print(sd(si[,"sil_width"]))
    cat('--------------------------\n\n') # cat() for newline instead of print()
    
    # Include clusters assignment in phyloseq object
    sample_data(data.now)$cluster<- as.factor(si[1:nsamples(data.now),1])
    sink()
  } # end-if all threshold are satisfied
  
  # Distance plot with definitive clusters
  pdf(paste("pcoa_definitiveClustering_",distanceMethod,"_k",kbestSI.mode,"_colorByCluster.pdf",sep=''),width=8)
  p <- plot_ordination(data.now, ordinate(data.now, "MDS", distance=dist), color='cluster') +
    geom_point(size=2.5) +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    theme(legend.position="bottom",legend.box="horizontal")
  print(p)
  dev.off()
  
  pdf(paste("pcoa_definitiveClustering_",distanceMethod,"_k",kbestSI.mode,".pdf",sep=''),width=8)
  p <- plot_ordination(data.now, ordinate(data.now, "MDS", distance=dist), shape='cluster', color=var.color) +
    geom_point(size=2.5) +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    theme(legend.position="bottom",legend.box="horizontal")
  print(p)
  dev.off()
  
  # Write output files
  data.norm <- data.now
  #   phyloseq object
  save(data.norm,dist,file=paste('data.normAndDist_definitiveClustering_',label,'.RData',sep=''))
  #   <sampleID,cluster> pairs
  sample_data(data.now)$sampleID <- row.names(sample_data(data.now))
  df.allVariables <- get_variable(data.now)
  df.clusterPerSample <- df.allVariables[,c('sampleID','cluster')]
  write.table(df.clusterPerSample,paste('sampleId-cluster_pairs_',label,'.txt',sep=''),quote=FALSE,sep=',',row.names=FALSE)
  
  return(data.norm)
} # end-function robust.clustering.decision



############################################
### Function robust.clustering.all.steps ###
############################################
# Runs all the required functions for getting a phyloseq object with a cluster
# assigned to each sample.
#
# Args:
#   path: where the phyloseq object RData file is. 
#   RDataOrBiomFile: the phyloseq object or the biom file.
#     If it is a RData file, the phyloseq object must call 'data.norm' and it must have its OTU counts already normalized!!
#     If it is a .biom file, the 5th argument is also required, because the .biom file will have only the OTU counts, which must be already normalized too!!
#   dataset.name: used as label and suffix of files. Example: 'David2014', 'Chicks', etc.
#   variable.in.PCoAgraph: name of a variable from sample_variables() in the phyloseq object, for the color of the samples in the PCoAgraph.
#   mapBiomFile: Mappping file instructions (only if biom format, else mapping info is included in the phyloseq object): comma separated values file, with samples in rows, being the first column the sampleID, and the remainder with their corresponding headers in the first row. The name to color the PCoA graphs must be one of these column headers within this mapping file.
# 
# Returns:
#   phyloseq object with a new variable in the phyloseq object ($cluster) with the cluster identifier per sample. This object also is saved in 'data.normAndDist_definitiveClustering_<dataset.label>.RData'. It could be used as input of other R scripts with posterior steps of microbiome dynamics analysis.
#   It also returns several text and graph files with the results.
robust.clustering.all.steps <- function(path,RDataOrBiomFile,dataset.label,variable.in.PCoAgraph,mapBiomFile){
  # Load libraries
  library(ggplot2,quietly=TRUE)
  library(phyloseq,quietly=TRUE)
  library(lattice,quietly=TRUE)
  library(cluster,quietly=TRUE)
  library(fpc,quietly=TRUE)
  library(reshape2,quietly=TRUE) # Load the reshape2 package (for the melt() function)
  
  # A) ROBUST CLUSTERING COMPUTATION
  # Test if phyloseq R object or biom file, depending on whether the 5th argument exists.
  if(missing(mapBiomFile)){ # RData phyloseq object
    load(RDataOrBiomFile)  
  }else{ # Biom + mapping file
    biom.otu<-import_biom(RDataOrBiomFile)
    map<-import_qiime_sample_data(mapBiomFile)
    data.norm <- merge_phyloseq(biom.otu,map)
  } # end-if RData or biom
  setwd(path)
  array4d.results <- robust.clustering(data.norm,variable.in.PCoAgraph,dataset.label)

  # B) PLOT GRAPHS
  load(paste('eval.arrays_',dataset.label,'.RData',sep=''))
  plot.robust.clustering.all.together.formatted(eval.array4d)

  # C) AUTOMATIC CLUSTERING DECISION
  data.norm <- robust.clustering.decision(data.norm,eval.array4d,variable.in.PCoAgraph,dataset.label)

  setwd("..")
  
  return(data.norm) 
} # end-function robust.clustering.all.steps

