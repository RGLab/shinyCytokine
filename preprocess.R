## preprocessing of data for cytokine visualization script

library(Kmisc) ## devtools::install_github("Kmisc", "kevinushey", subdir="Kmisc")
library(data.table.extras) ## devtools::install_github("data.table.extras", "kevinushey")
library(reshape2)
library(testthat)
library(plyr)

source("common_functions.R")

## debug
dat_file <- "data/res_cd4.rds"
meta_file <- "data/meta.rds"
name <- "cd4"

preprocess <- function(dat_file, meta_file, name) {
  
  dat <- readRDS(dat_file)
  meta <- readRDS(meta_file)
  
  stopifnot( all( names(dat) == meta$name ) )
  
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
  
  ## ensure each entry is a matrix
  stopifnot( identical( names( counts( sapply( dat, class ) ) ), "matrix" ) )
  
  ## find negatives
  stopifnot(
    all( sapply(dat, function(x) {
      all( x >= 0 )
    }))
  )
  
  ## generate all possible permutations
  cytokines <- unname( unclass( colnames( dat[[1]] ) ) )
  saveRDS(cytokines, "data/cytokines.rds")
  combos <- combinations( 1:length(cytokines) )
  ncol <- sum( choose( length(cytokines), 1:length(cytokines) ) )
  colnames <- combinations(cytokines)
  
  ## re-order each matrix to be the same column order
  dat[] <- lapply(dat, function(x) {
    x[, cytokines, drop=FALSE]
  })
  
  ## test all the names
  tmp <- unlist( colnames[1:6] )
  for( i in 1:length(dat) ) {
    stopifnot( all( tmp == colnames(dat[[i]]) ) )
  }
  
  ## make sure they're identical across 'combos', 'colnames'
  for( i in 1:length(combos) ) {
    for( j in 1:length(dat) ) {
      stopifnot( all( colnames(dat[[j]])[combos[[i]]] == colnames[[i]] ) )  
    }
  }
  
  ## collapse the colnames
  colnames <- sapply(colnames, paste, collapse=" x ")
  
  ## generate proportions for all cytokine combinations
  ## first, generate raw TRUE/FALSE entries for each cell
  ## in each sample; TRUE if all cytokines of a particular combo are TRUE,
  ## FALSE if all cytokines of a particular combo are FALSE
  d <- vector("list", length(dat))
  for( i in 1:length(dat) ) {
    x <- dat[[i]]
    output <- matrix(FALSE, nrow=nrow(x), ncol=ncol)
    for( j in 1:ncol(output) ) {
      output[, j] <- as.logical(apply( x[, combos[[j]], drop=FALSE], 1, prod))
    }
    colnames(output) <- colnames
    d[[i]] <- output
  }
  
  ## compute proportions
  d_prop <- do.call( rbind, 
    lapply(d, function(x) {
      colApply(x, mean)
    })
  )
  
  overall_sample_prop <- rowMeans(d_prop)
  
  ## compute the 'degree of functionality'
  ## for each cell, we compute the number of cytokines
  ## that have been expressed simultaneously
  ## we then collapse over all cells in a sample
  dof <- sapply(dat, function(x) {
    rowSums(x > 0)
  })
  
  stopifnot( all( names(dof) == meta$name ) )
  
  ## merge the data sets together
  merged <- cbind(d_prop, meta)
  
  ## melt for ggplot2 usage
  dat_prop_melt <- melt(
    merged, 
    id.vars=names(meta),
    value.name="Proportion", 
    variable.name="Cytokine"
  )
  
  ## cytokine order
  dat_prop_melt$Cytokine_order <- 
    stringr::str_count( dat_prop_melt$Cytokine, " x " ) + 1L
  
  ## more informative name for the samples
  names(dat_prop_melt)[ names(dat_prop_melt) == "name" ] <- "Sample"
  
  ## perform background subtraction on the proportions
  ## ie, for a given PTID x VISIT combo, subtract the
  ## neg. controls
  d_bg <- ddply( dat_prop_melt, .(PTID, VISITNO), function(DF) {
    if( length( grep("^negctrl", DF$Stim) ) == 0 ) {
      warning("No negative controls for this sample! PTID = ", DF$PTID[1], ", VISITNO = ", DF$VISITNO[1])
      DF$Proportion_bg <- DF$Proportion
    } else {
      DF$Proportion_bg <-
        DF$Proportion - mean( DF$Proportion[ grep("^negctrl", DF$Stim) ] )
    }
    return(DF)
  })
  
  ## FIXME: Add / remove negctrl?
  ## remove sebctrl, negctrl
#   d2 <- droplevels( d_bg[
#     d_bg$Stim %nin% grep("ctrl", unique(d_bg$Stim), value=TRUE),
#   ])
  d2 <- d_bg
  
  ## an individual-level data set
  d <- data.table(dat_prop_melt)
  del(d, VISITNO, Stim, Sample)
  setkeyv(d, c("PTID", "Cytokine"))
  d <- d[, Proportion := mean(Proportion), by=list(PTID, Cytokine)]
  d <- group_subset(d, by=list(PTID, Cytokine))
  
  ## write the files out
  dir.create( paste0("data/", name), showWarnings=FALSE )
  saveRDS( dat_prop_melt, file=paste0("data/", name, "/sample_proportions_postProcess_marginal.rds") )
  saveRDS( d, file=paste0("data/", name, "/indiv_proportions_postProcess_marginal.rds") )
  saveRDS( d2, file=paste0("data/", name, "/sample_proportions_bg_subtracted_marginal.rds") )
  saveRDS( dof, file=paste0("data/", name, "/cell_dof_marginal.rds") )
  dat <- dat_prop_melt[ dat_prop_melt$Cytokine %in% cytokines, ]
  saveRDS( dat, file=paste0("data/", name, "/sample_proportions_postProcess_noCombo_marginal.rds") )
  
}

preprocess("data/res_cd4.rds", "data/meta.rds", "cd4")
preprocess("data/res_cd8.rds", "data/meta.rds", "cd8")