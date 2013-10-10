library(Kmisc) ## devtools::install_github("Kmisc", "kevinushey", subdir="Kmisc")
library(data.table.extras) ## devtools::install_github("data.table.extras", "kevinushey")
library(reshape2)
library(testthat)
library(plyr)

source("common_functions.R")

## debug
dat_file <- "data/res_cd8.rds"
meta_file <- "data/meta.rds"
counts_file <- "data/cd8_counts.rds"
name <- "cd8_spice"

make_spice_file <- function(dat_file, meta_file, counts_file, name) {
  
  dat <- readRDS(dat_file)
  meta <- readRDS(meta_file)
  counts <- readRDS(counts_file)
  cytokines <- unclass( unname( colnames(dat[[1]] ) ) )
  
  ## remove the NULL cases
  non_null_cases <- which( sapply( dat, function(x) {
    !is.null(x)
  }))
  
  if( length(non_null_cases) > 0 ) {
    warning("The following .fcs files are NULL:\n\t", 
      paste( names(dat)[non_null_cases], collapse=", ")
    )
  }
  
  dat <- dat[ non_null_cases ]
  meta <- droplevels(meta[ non_null_cases, ])
  
  ## ensure the matrices are all in the same order
  nm <- colnames(dat[[1]])
  for (i in seq_along(dat)) {
    stopifnot( all( colnames(dat[[i]]) %in% nm ) )
  }
  dat[] <- lapply(dat, function(x) {
    x[,nm, drop=FALSE]
  })
  
  ## assign informative row names
  for (i in 1:length(dat)) {
    rownames( dat[[i]] ) <- paste( names(dat)[i], 1:nrow(dat[[i]]), sep=" : " )
    colnames( dat[[i]] ) <- unclass( colnames( dat[[i]] ) )
  }
  
  long_dat <- do.call( rbind, dat )
  long_dat_bin <- colApply(long_dat, function(x) {
    return(x>0)
  })
  
  long_dat_bin <- colApply(long_dat_bin, function(x) {
    swap(x, c(FALSE, TRUE), c("-", "+"))
  })
  
  rownames(long_dat_bin) <- rownames(long_dat)
  colnames(long_dat_bin) <- unname( unclass( colnames( dat[[1]] ) ) )
  attr(long_dat_bin, "names") <- NULL
  
  head(long_dat_bin)
  head(long_dat)
  
  tmp <- data.frame(
    Subject=rownames(long_dat_bin),
    stringsAsFactors=FALSE
  )
  
  tmp <- factor_to_char(cbind(tmp, long_dat_bin))
  tmp$Sample <- gsub(" :.*", "", rownames(tmp))
  tmp <- data.frame(
    Sample=tmp$Sample,
    Cytokines=do.call( function(...) paste(..., sep=','), tmp[cytokines] ),
    stringsAsFactors=FALSE
  )
  tmp <- data.table(tmp)
  tmp[, counts := counts(Cytokines), by=Sample]
  out <- group_subset(tmp, by=list(Sample, Cytokines))
  
  totals <- as.data.frame(counts)
  totals$Sample <- rownames(totals)
  names(totals) <- c("total", "Sample")
  out <- merge(out, totals, by="Sample")
  
  out$prop <- out$counts / out$total
  keep(out, Sample, Cytokines, prop)
  
  ## trying tabs
  out$Cytokines <- gsub(",", "\t", out$Cytokines)
  
  ## write it out
  dir.create(file.path("data", name), showWarnings=FALSE)
  file <- paste0("data/", name, "/", name, ".txt")
  header <- paste( "Sample", paste(cytokines, collapse="\t"), "value", sep="\t" )
  cat( paste0(header, "\n"), file=file )
  write.table(out,
    file=file,
    quote=FALSE,
    sep="\t",
    row.names=FALSE,
    col.names=FALSE,
    append=TRUE
  )
  
}

make_spice_file(
  "data/res_cd4.rds",
  "data/meta.rds",
  "data/cd4_counts.rds",
  "cd4_spice"
)

make_spice_file(
  "data/res_cd8.rds",
  "data/meta.rds",
  "data/cd8_counts.rds",
  "cd8_spice"
)
