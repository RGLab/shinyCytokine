## preprocessing of data for cytokine visualization script

library(car)
library(Kmisc) ## devtools::install_github("Kmisc", "kevinushey", subdir="Kmisc")
library(data.table.extras) ## devtools::install_github("data.table.extras", "kevinushey")
library(reshape2)
library(testthat)
library(plyr)
library(flowCore)

source("common_functions.R")

## debug
dat_file <- "data/RV144 Pilot Study/res_cd4.rds"
meta_file <- "data/RV144 Pilot Study/meta.rds"
counts_file <- "data/RV144 Pilot Study/cd4_counts.rds"
name <- "cd4_joint"
generate_proportions <- marginal

preprocess <- function(dat_file, meta_file, counts_file, name, generate_proportions) {
  
  generate_proportions <- match.fun(generate_proportions)
  
  dat <- readRDS(dat_file)
  meta <- readRDS(meta_file)
  counts <- readRDS(counts_file)
  
  ## remove FCS files with abnormally low cell counts
  low_cell_counts <- names(counts)[ counts < 1E4 ]
  
  dat <- dat[ names(dat) %nin% low_cell_counts ]
  meta <- droplevels(meta[ meta$name %nin% low_cell_counts, ])
  counts <- counts[ names(counts) %nin% low_cell_counts ]
  
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
  counts <- counts[ non_null_cases ]
  
  stopifnot( all( names(dat) == names(counts) ) )
  
  ## ensure each entry is a matrix
  stopifnot( identical( names( Kmisc::counts( sapply( dat, class ) ) ), "matrix" ) )
  
  ## find negatives
  stopifnot(
    all( sapply(dat, function(x) {
      all( x >= 0 )
    }))
  )
  
  ## generate all possible permutations
  cytokines <- unname( unclass( colnames( dat[[1]] ) ) )
  saveRDS(cytokines, "data/RV144 Pilot Study/cytokines.rds")
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
  d <- generate_proportions(dat, combos, ncol, colnames)
  
  ## compute proportions
  ## use the total counts to compute proportions
  d_prop <- do.call( rbind, 
    lapply(d, function(x) {
      colApply(x, sum)
    })
  )
  
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
  d <- melt(
    merged, 
    id.vars=names(meta),
    value.name="Counts", 
    variable.name="Cytokine"
  )
  
  ## cytokine order
  d$Cytokine_order <- 
    stringr::str_count( d$Cytokine, " x " ) + 1L
  
  ## more informative name for the samples
  names(d)[ names(d) == "name" ] <- "Sample"
  
  ## merge in cell counts
  tmp <- data.frame(
    Sample=names(counts),
    TotalCellCount=counts
  )
  d <- merge(d, tmp, by="Sample", all.x=TRUE)
  rm(tmp)
  
  ## compute 'number of cells expressing at least one cytokine'
  tmp <- sapply(dat, function(x) {
    sum( apply(x, 1, function(xx) {
      any( xx > 0 )
    } ) )
  })
  tmp <- data.frame(
    Sample=names(tmp),
    TotalActivated=tmp
  )
  d <- merge( d, tmp, by="Sample", all.x=TRUE)
  rm(tmp)
  
  d <- data.table(d)
  
  ## Compute the proportions
  d[, PropTotal := Counts / TotalCellCount]
  d[, PropActivated := Counts / TotalActivated]
  
  ## BG Correction
  d[, 
    PropTotalBG := PropTotal - mean(PropTotal[ grep("^negctrl", Stim) ]), 
    by=list(PTID, VISITNO, Cytokine)
  ]
  
  d[, 
    PropActivatedBG := PropActivated - mean(PropActivated[ grep("^negctrl", Stim) ]), 
    by=list(PTID, VISITNO, Cytokine)
  ]
  
  ## correct for missingness
  negmean <- mean( d[ 
    Stim %in% c("negctrl 1", "negctrl 2"), 
    mean(PropTotal), 
    by=list(PTID, VISITNO, Cytokine) 
  ]$V1 )
  
  d[ 
    is.na(PropTotalBG), 
    PropTotalBG := PropTotal - negmean, 
    by=list(PTID, VISITNO, Cytokine) 
  ]
  
  negmean <- mean( 
    d[ 
      Stim %in% c("negctrl 1", "negctrl 2"), 
      mean(PropActivated), 
      by=list(PTID, VISITNO, Cytokine) 
    ]$V1 
  )
  
  d[ 
    is.na(PropActivatedBG), 
    PropActivatedBG := PropActivated - negmean, 
    by=list(PTID, VISITNO, Cytokine) 
  ]
  
  ## Log scale
  ## Compute log fold changes instead of the subtraction
  ## log(w_s / w_u), s=stimulated, u=unstimulated)
  d[,
    LogFoldChange := log2(
      (Counts+1) / mean(Counts[ Stim == "negctrl 1" ] + 1)
    ),
    by=list(PTID, VISITNO, Cytokine)
  ]
  
  ## write the files out
  dir.create( paste0("data/RV144 Pilot Study/", name), showWarnings=FALSE )
  ## saveRDS( d, file=paste0("data/", name, "/sample_proportions_postProcess.rds") )
  ## saveRDS( d, file=paste0("data/", name, "/indiv_proportions_postProcess.rds") )
  saveRDS( d, file=paste0("data/RV144 Pilot Study/", name, "/", name, "_sample_proportions_bg_subtracted.rds") )
  ## saveRDS( dof, file=paste0("data/", name, "/cell_dof_marginal.rds") )
  ## dat <- d[ d$Cytokine %in% cytokines, ]
  ## saveRDS( dat, file=paste0("data/", name, "/sample_proportions_postProcess_noCombo_marginal.rds") )
  
}

## CD4
preprocess("data/RV144 Pilot Study/res_cd4.rds", "data/RV144 Pilot Study/meta.rds", "data/RV144 Pilot Study/cd4_counts.rds", "cd4_marginal", marginal)
preprocess("data/RV144 Pilot Study/res_cd4.rds", "data/RV144 Pilot Study/meta.rds", "data/RV144 Pilot Study/cd4_counts.rds", "cd4_joint", joint)

## CD8
preprocess("data/RV144 Pilot Study/res_cd8.rds", "data/RV144 Pilot Study/meta.rds", "data/RV144 Pilot Study/cd8_counts.rds", "cd8_marginal", marginal)
preprocess("data/RV144 Pilot Study/res_cd8.rds", "data/RV144 Pilot Study/meta.rds", "data/RV144 Pilot Study/cd8_counts.rds", "cd8_joint", joint)