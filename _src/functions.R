# All functions used across the Rmd's in this dir. 

chroms <- c(as.character(1:19), "X")
interp_bp <- function(df) {
  chroms <- c(as.character(1:19), "X")
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
