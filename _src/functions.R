# All functions used across the Rmd's in this dir. 

# copied from: https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, 
                   attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}


chroms <- c(as.character(1:19), "X")
interp_bp <- function(df) {
  
  df <- arrange(df, peak_chr, peak_cM)
  peak_gpos <- select(df, peak_chr, peak_cM)
  chr <- peak_gpos$peak_chr
  f <- factor(chr, chroms)
  peak_gcoord_list <- split(peak_gpos$peak_cM, f)
  peak_pcoord_list <- qtl2::interp_map(peak_gcoord_list, gmap, pmap)
  df$interp_bp_peak <- unsplit(peak_pcoord_list, f)
  df
}


rankZ <- function (x) {
  x <- rank(x, na.last = "keep", ties.method = "average")/(sum(!is.na(x)) + 1)
  qnorm(x)
}

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

do_global_mediation <- function(gene.name,p.chr, peak.info, probs, expr.target, expr.mediator, covar, med.type){
  
  gene <- filter(peak.info, mgi_symbol==gene.name & peak_chr==p.chr)[1,]
  #if(dim(gene)[1] >1){ print("check the peaks!")}
  
  if(med.type=="p-p"){ # protein to protein
    meds <- all.prots %>% left_join(all.genes)
    meds <- meds[meds$protein_id %in% colnames(expr.mediator),]
    mediator <- expr.mediator[ ,meds$protein_id, drop=FALSE ]
    target   <- expr.target[ ,gene$protein_id,drop=FALSE ]
    marker   <- map_dat2 %>% filter(pos_cM == gene$peak_cM.esc_prot)
  }
  
  if(med.type=="r-p"){ # rna to protein
    meds <- all.genes 
    meds <- meds[meds$ensembl_gene_id %in% colnames(expr.mediator),]    
    mediator <- expr.mediator[ ,meds$ensembl_gene_id, drop=FALSE ]
    target   <- expr.target[ ,gene$protein_id,drop=FALSE ]
    marker   <- map_dat2 %>% filter(pos_cM == gene$peak_cM.esc_prot)
  }
  if(med.type=="r-r"){ # rna to rna
    meds <- all.genes
    meds <- meds[meds$ensembl_gene_id %in% colnames(expr.mediator),]
    mediator <- expr.mediator[ ,meds$ensembl_gene_id, drop=FALSE ]
    target   <- expr.target[ ,gene$ensembl_gene_id,drop=FALSE ]
    marker   <- map_dat2 %>% filter(pos_cM == gene$peak_cM.esc_rna)
  }
  if(med.type=="p-r"){ # protein to rna
    meds <- all.prots %>% left_join(all.genes)
    meds <- meds[meds$protein_id %in% colnames(expr.mediator),]
    mediator <- expr.mediator[ ,meds$protein_id, drop=FALSE ]
    target   <- expr.target[ ,gene$ensembl_gene_id,drop=FALSE ]
    marker   <- map_dat2 %>% filter(pos_cM == gene$peak_cM.esc_rna)
  }
  
  
  annot    <- meds %>% mutate(chr=gene_chr,pos=abs(gene_end+gene_start)/2)
  geno     <- pull_genoprobpos(probs,marker$marker)
  
  med.scan <- mediation.scan(target= target,
                             mediator = mediator,
                             annotation = annot,
                             covar =  covar,
                             qtl.geno = geno, verbose=FALSE) 
  med.scan <- med.scan %>% select(-chr) %>% mutate(target=gene.name) %>% 
    left_join(.,all.genes2) %>% rename("mediator"="mgi_symbol")
  return(med.scan)
}

plot_mediation <-function(med.scan, p.chr){
  xax <- list(tickmode="array",tickvals=chrom_lens_midpt, ticktext= names(chrom_lens), title="Chr")
  med.scan <- med.scan %>% mutate(chr=ifelse(gene_chr=="X", 20,chrom))
  p<-plot_ly(data=med.scan,x=~cumsum_bp_gene,y=~LOD,hoverinfo="text",type="scatter",mode="markers",
             text=med.scan$mediator,color =~(as.integer(chr) %% 2 == 0),
             alpha = 0.6) %>% layout(showlegend = FALSE, title = paste0(med.scan$target[1]," mediation for peak on chr ",p.chr ), 
                                     xaxis=xax) 
  return(p)
}

plot_mediation2 <-function(med.scan, p.chr){
  xax <- list(tickmode="array",tickvals=chrom_lens_midpt, ticktext= names(chrom_lens), title="Chr")
  med.scan <- med.scan %>% mutate(chr=ifelse(mediator.chr=="X", 20,mediator.chr))
  p<-plot_ly(data=med.scan,x=~cumsum_bp_gene,y=~mediation.z,hoverinfo="text",type="scatter",mode="markers",
             text=~mediator.symbol,color =~(as.integer(chr) %% 2 == 0),
             alpha = 0.6) %>% layout(showlegend = FALSE, title = paste0(~target.symbol[1]," mediation for peak on chr ",~p.chr ), 
                                     xaxis=xax) 
  return(p)
}


getplots <-function(g.name,probs,exprZ,kinship_loco,covar,p.chr, dist=FALSE,thr=7.5, prot=FALSE, ...){
  if(prot==FALSE){
    gene <- filter(all.genes2, mgi_symbol==g.name)
    gene.id <- gene$ensembl_gene_id_v85
  }
  if(prot ==TRUE){ 
    prot <- filter(all.prots, protein_id_v98==g.name)
    gene <- filter(all.genes2, ensembl_gene_id_v85 == prot$ensembl_gene_id_v85)
    gene.id <- g.name
    
  }
  out <- scan1(probs, exprZ[,gene.id , drop=FALSE],
               kinship_loco, addcovar=covar )
  peak <- find_peaks(out, gmap, drop=1.5,
                     threshold=thr, thresholdX=thr) 
  if(p.chr %in% c(1:19) & dim(peak) > 0){
    peak <- peak %>% filter(chr == p.chr)
  }
  for(i in 1:dim(peak)[1]){
    chrom <- peak$gene_chr[i]
    eff <- scan1blup(probs[,chrom], exprZ[,gene.id , drop=FALSE],
                     kinship_loco[[chrom]], addcovar=covar )
    # Let's find the genomic position - marker that is closest to my gene to abline
    gene.pos <- map_dat2[map_dat2$chr==gene$chrom,][which.min(abs(map_dat2[map_dat2$chr==gene$chrom,]$pos_bp-gene$midpoint)),]
    peak.pos <- peak$pos[i]
    par(mfrow=c(2,1))
    plot(out, map = gmap, chr = chrom , gap=0, col="dark green", xlab="cM", bgcolor="gray95",
         main=paste0(gene$mgi_symbol," (",gene.id,") "," plot"))
    if(gene$chrom == chrom){ abline(v = gene.pos$pos_cM, col="black", lwd=2,lty=2)}
    abline(h=thr,lwd=1)
    abline(v=peak.pos,col="red",lwd=2,lty=3)
    if(peak.pos >30){plot_coefCC(eff, map = gmap[chrom], main=paste0(gene$mgi_symbol," (",gene.id,") "," plot"),bgcolor="gray95",legend="bottomleft")}
    if(peak.pos <30){plot_coefCC(eff, map = gmap[chrom], main=paste0(gene$mgi_symbol," (",gene.id,") "," plot"),bgcolor="gray95",legend="bottomright")}
    if(gene$chrom == chrom){ abline(v = gene.pos$pos_cM, col="black", lwd=2,lty=2)}
    abline(v=peak.pos,col="red",lwd=2,lty=3)
  }
}

gene_plot <-  function(g.obj,source,term){
  
  term.name <- strsplit(strsplit(term, "Factor: ")[[1]][2],";")[[1]][1]
  genes <- tibble(ensembl_gene_id = unlist(str_split((g.obj$result %>% filter(term_name ==term))$intersection,","))) %>%
    left_join(., all.genes) %>% mutate(category = term.name) %>% select(mgi_symbol,category)
  p.value <- formatC((g.obj$result %>% filter(term_name ==term))$p_value, digits=2)
  
  
  g.graph <- igraph::graph_from_data_frame(genes , directed=FALSE)
  
  
  g1 <- ggraph(g.graph)+
    geom_edge_link0(alpha=0.1)+
    geom_node_point(size=5)+
    geom_node_text(aes(label = name), repel = TRUE )+
    theme_classic()+
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())+
    ggtitle(paste0(source," database, p-value < ", p.value))
  return(g1)
  
}

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
# from http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software


# Making downloadable data tables
# https://www.r-bloggers.com/vignette-downloadable-tables-in-rmarkdown-with-the-dt-package/
create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                rownames = FALSE, 
                filter="top",
                options = list(dom = 'Blfrtip',
                               buttons = c('copy', 'csv', 'excel'),
                               pageLength = 5, 
                               scrollX= TRUE
                               ))
  
}

create_dt_alt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                rownames = FALSE, 
                filter="top",
                options = list(dom = "<Blf<\"datatables-scroll\"t>ipr>",
                               buttons = c('copy', 'csv', 'excel'),
                               pageLength = 5, 
                               scrollX = TRUE,
                               scrollY = 400,
                               autoWidth = TRUE,
                               columnDefs = list(list(width = '800px', targets = c(4)) )
                ))
  
}


get_peaks <- function(pc, probs, kinship,covar.mat){
  #rownames(pc) <- rownames(atac.pca$x)
  scan1 <- scan1(genoprobs = probs, pheno = pc, kinship = kinship, addcovar = covar.mat)
  peaks.lod5 <- find_peaks(scan1, threshold = 5,map = gmap)
  return(as_tibble(peaks.lod5))
}

get_meds <- function(target, meds, mediator, peaks, probs,covar, z_thres = -4,  pos_thres = 10 ){
  med.scan.all <- c()
  samples <- rownames(mediator)
  for( i in 1:nrow(peaks)){
    #print( paste0("i is ", i))
    
    marker    <- map_dat2 %>% filter(pos_cM == peaks$pos_cM[i])
    annot     <- meds %>% mutate(chr=gene_chr,pos=abs(gene_end+gene_start)/2)
    geno      <- pull_genoprobpos(probs,marker$marker)
    geno      <- geno[samples,]
    
    if( !is.null(covar)){ 
      covar <- covar[samples,,drop=FALSE] 
    }
    
    med.scan <- mediation_scan(target= target[samples,peaks$phenotype[i], drop=FALSE],
                               mediator = mediator,
                               annotation = annot,
                               covar =  covar,
                               driver = geno,
                               verbose=FALSE,
                               method     = "double-lod-diff") 
    
    med.scan <- med.scan %>% 
      mutate(phenotype=peaks$phenotype[i], 
             peak_chr = peaks$peak_chr[i], 
             peak_lod = peaks$lod[i],
             med_chr = chr) %>% 
      mutate(scaled_LOD = scale(lod), 
             middle = (gene_end+gene_start)/2e06) %>%
      filter( (scaled_LOD < z_thres & 
                 peak_chr ==  med_chr & 
                 abs(middle - pos/1e06) <= pos_thres) |
                (lod < peak_lod*0.5))
    
    med.scan.all <- rbind(med.scan.all,med.scan)
  }
  return(med.scan.all)
}


# got from: https://github.com/federicogiorgi/aracne.networks/blob/master/R/write.regulon.R
write.regulon<-function(
  regulon,
  file="",
  sep="\t",
  header=TRUE,
  n=Inf,
  regulator=NULL
){
  if(header){
    cat(paste0("Regulator\tTarget\tMoA\tlikelihood\n"),file=file)
  }
  nn<-0
  if(is.null(regulator)){
    for(tf in names(regulon)){
      x<-regulon[[tf]]
      targets<-names(x$tfmode)
      moas<-x$tfmode
      likelihoods<-x$likelihood
      tab<-cbind(rep(tf,length(targets)),targets,moas,likelihoods)
      i<-1
      while(nn<n&i<=nrow(tab)){
        cat(
          tab[i,],
          file=file,
          append=TRUE,
          sep=sep
        )
        cat("\n",file=file,append=TRUE)
        nn<-nn+1
        i<-i+1
      }
    }
  } else {
    tf<-regulator
    x<-regulon[[tf]]
    targets<-names(x$tfmode)
    moas<-x$tfmode
    likelihoods<-x$likelihood
    tab<-cbind(rep(tf,length(targets)),targets,moas,likelihoods)
    i<-1
    while(nn<n&i<=nrow(tab)){
      cat(
        tab[i,],
        file=file,
        append=TRUE,
        sep=sep
      )
      cat("\n",file=file,append=TRUE)
      nn<-nn+1
      i<-i+1
    }
  }
  
}

subset_probs <- function(this_probs, this_chrom, this_markers) {
  att <- attributes(this_probs)
  att$names <- this_chrom
  att$is_x_chr <- setNames(FALSE, this_chrom)
  #assert_that(all(this_markers %in% dimnames(this_probs[[this_chrom]])[[3]]))
  newprobs <- list(this_probs[[this_chrom]][, , this_markers, drop=FALSE])
  names(newprobs) <- this_chrom
  attributes(newprobs) <- att
  newprobs
}
