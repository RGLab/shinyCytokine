## preprocessing of data for cytokine visualization script

suppressPackageStartupMessages({
  library(Rcpp)
  library(stringr)
  library(car)
  library(Kmisc) ## devtools::install_github("Kmisc", "kevinushey", subdir="Kmisc")
  library(data.table.extras) ## devtools::install_github("data.table.extras", "kevinushey")
  library(reshape2)
  library(testthat)
  library(plyr)
  library(flowCore)
  library(COMPASS)
})

source("common_functions.R")
sourceCpp("preprocess_scripts/generate_proportions.cpp")

## debug
dat_file <- "data/RV144 CaseControl/res.rds"
meta_file <- "data/RV144 CaseControl/meta.rds"
counts_file <- "data/RV144 CaseControl/cd4_counts.rds"
name <- "cd4"

preprocess <- function(dat_file, meta_file, counts_file, name) {
  
  dat <- readRDS(dat_file)
  meta <- readRDS(meta_file)
  counts <- readRDS(counts_file)
  
  ## fit a mixture model to the counts to determine which are too 'low'
  ## and which are okay
  
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
  
  null_cases <- which( sapply(dat, is.null) )
  
  if( length(null_cases) > 0 ) {
    warning("The following .fcs files are NULL:\n\t", 
      paste( names(dat)[null_cases], collapse=", ")
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
  saveRDS(cytokines, "data/RV144 CaseControl/cytokines.rds")
  combos <- combinations( c(1:length(cytokines), -1:-length(cytokines)) )
  combos <- combos[ sapply(combos, function(x) !any( duplicated( abs(x) ) ) ) ]
  ncol <- length(combos)
  colnames <- sapply(combos, function(x) {
    paste( collapse="",
      swap(x, c(1:6, -1:-6), c(paste0(cytokines, "+"), paste0(cytokines, "-")))
    )
  })
  
  ## re-order each matrix to be the same column order
  dat[] <- lapply(dat, function(x) {
    x[, cytokines, drop=FALSE]
  })
  
  ## test all the names
  for( i in 1:length(dat) ) {
    stopifnot( all( cytokines == colnames(dat[[i]]) ) )
  }
  
  ## convert to logical (on/off)
  dat <- lapply(dat, function(x) {
    colApply(x, drop=FALSE, function(x) {
      as.logical(x)
    })
  })
  
  ## generate proportions for all cytokine combinations
  d_counts <- generate_proportions(dat, combos)
  colnames(d_counts) <- colnames
  rownames(d_counts) <- names(dat)
  
  ## a test
  i <- 700
  c_combo <- combos[i]
  c_d <- d_counts[,i]
  j <- which.max(c_d)
  tmp <- rowApply( dat[[j]], function(xx) {
    if (xx[1] == FALSE && xx[2] == TRUE && xx[3] == TRUE && xx[4] == TRUE && xx[5] == FALSE && xx[6] == FALSE) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  })
  stopifnot( identical( sum(tmp), c_d[[j]] ) )
  
  ## compute the 'degree of functionality'
  ## for each cell, we compute the number of cytokines
  ## that have been expressed simultaneously
  ## we then collapse over all cells in a sample
  dof <- sapply(dat, function(x) {
    rowSums(x > 0)
  })
  
  stopifnot( all( names(dof) == meta$name ) )
  
  ## merge the data sets together
  merged <- cbind(d_counts, meta)
  
  ## melt for ggplot2 usage
  d <- melt(
    merged, 
    id.vars=names(meta),
    value.name="Counts", 
    variable.name="Cytokine"
  )
  
  ## more informative name for the samples
  names(d)[ names(d) == "name" ] <- "Sample"
  
  ## cytokine order
  d$CytokineOrder <- stringr::str_count(d$Cytokine, "[-\\+]")
  
  ## merge in cell counts
  tmp <- data.frame(
    Sample=names(counts),
    TotalCellCount=counts
  )
  d <- merge(d, tmp, by="Sample", all.x=TRUE)
  rm(tmp)
  
  ## compute 'number of cells expressing at least one cytokine'
  tmp <- sapply(dat, nrow)
  tmp <- data.frame(
    Sample=names(tmp),
    TotalActivated=tmp
  )
  d <- merge(d, tmp, by="Sample", all.x=TRUE)
  rm(tmp)
  
  ## denote whether a particular cytokine is on or off for a particular
  ## combination
  
  d <- data.table(d)
  for (i in seq_along(cytokines)) {
    cytokine <- cytokines[i]
    d[[ cytokine ]] <- 0L
    d[[ cytokine ]][ str_detect(d$Cytokine, paste0(cytokine, "+")) ] <- 2
    d[[ cytokine ]][ str_detect(d$Cytokine, paste0(cytokine, "-")) ] <- 1
  }
  
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
  negmean <- mean( 
    d[ 
      Stim %in% c("negctrl 1", "negctrl 2"), 
      mean(PropTotal), 
      by=list(PTID, VISITNO, Cytokine) 
    ]$V1 
  )
  
  d[ 
    is.na(PropTotalBG), 
    PropTotalBG := PropTotal - negmean, by=list(PTID, VISITNO) 
  ]
  
  negmean <- mean( d[ Stim %in% c("negctrl 1", "negctrl 2"), mean(PropActivated), by=list(PTID, VISITNO) ]$V1 )
  d[ is.na(PropActivatedBG), PropActivatedBG := PropActivated - negmean, by=list(PTID, VISITNO) ]
  
  ## Log scale
  ## Compute log fold changes instead of the subtraction
  ## log(w_s / w_u), s=stimulated, u=unstimulated)
  d[,
    LogFoldChange := log2(
      (Counts+1) / (Counts[ Stim == "negctrl 1" ] + 1)
    ),
    by=list(PTID, VISITNO, Cytokine)
  ]
  
  ## merging in new data from Raphael
  ## from Raphael -- there are some 'bad ptids' to exclude.
  bad_ptid_old <- c(523845, 525881, 600929, 616848, 631338, 631772, 644891, 646477, 704809, 719186, 849938, 104173, 125948, 132409, 134551, 210617, 211497, 217402, 223868, 233547, 323142, 330034, 330272, 339553)
  master <- read.csv("data/RV144 CaseControl/rv144_master_wk26.csv")
  correlates <- read.csv("data/RV144 CaseControl/rv144_data_wk26_correlates_long_mapped.csv")
  
  bad_ptid <- as.character(unique(correlates$PTID[ correlates$PTID_orig %in% bad_ptid_old ]))
  
  ## merge in the infection status from master
  master <- data.table(master)
  master$PTID <- swap(master$pin, correlates$PTID_orig, as.character(correlates$PTID))
  master$PTID[ master$PTID %nin% correlates$PTID ] <- NA
  master <- master[ !is.na(PTID), ]
  keep(master, PTID, site, vaccno, infect, time_to_infect, perprot, ittmod, BRA_risk)
  
  setkey(d, PTID)
  setkey(master, PTID)
  d <- master[d]
  
  d <- d[ PTID %nin% bad_ptid, ]
  
  ## write the files out
  dir.create( paste0("data/RV144 CaseControl/", name), showWarnings=FALSE )
  ## saveRDS( d, file=paste0("data/", name, "/sample_proportions_postProcess.rds") )
  ## saveRDS( d, file=paste0("data/", name, "/indiv_proportions_postProcess.rds") )
  saveRDS( d, file=paste0("data/RV144 CaseControl/", name, "/", name, "_sample_proportions_bg_subtracted.rds") )
  ## saveRDS( dof, file=paste0("data/", name, "/cell_dof_marginal.rds") )
  ## dat <- d[ d$Cytokine %in% cytokines, ]
  ## saveRDS( dat, file=paste0("data/", name, "/sample_proportions_postProcess_noCombo_marginal.rds") )
  
  saveRDS(d_counts, file=paste0("data/RV144 CaseControl/", name, "/", name, "_cytokine_combos_preprocessed.rds"))
  
}

## CD4
preprocess("data/RV144 CaseControl/res.rds", "data/RV144 CaseControl/meta.rds", "data/RV144 CaseControl/cd4_counts.rds", "cd4")

## A function for merging a set of results (as retrieved from Lynn)

## debug
data_loc <- "data/RV144 CaseControl/cd4/cd4_sample_proportions_bg_subtracted.rds"
results_loc <- "data/RV144 CaseControl/Paperdata_RV144_CC.RData"
map_loc <- "data/RV144 CaseControl/PTID_map_case_control.csv"

merge_results <- function(data_loc, results_loc, map_loc) {
  
  ## get the things we need
  data <- readRDS(data_loc)
  env <- new.env()
  load(results_loc, envir=env)
  
  ## we need to make sure we're using cytokine names in the same order
  ## as the previous preprocessing script
  cyto <- unlist(strsplit( as.character(data$Cytokine)[ data$CytokineOrder == 6 ][1], "[-\\+]"))
    
  Mgamma <- env$Mgamma
  map <- read.csv(map_loc, header=TRUE, sep=",", stringsAsFactors=FALSE)
  rownames(Mgamma) <- swap( rownames(Mgamma), map$PTID_orig, map$PTID )
  
  ## we want to fill Mgamma with zeroes for the cytokines not included
  combos <- do.call( paste0, (do.call( expand.grid, replicate(6, c(0, 1), simplify=FALSE) ) ) )
  for (combo in combos) {
    if (combo %nin% colnames(Mgamma)) {
      Mgamma[combo] <- 0
    }
  }
  
  Mgamma <- Mgamma[ order(colnames(Mgamma)) ]
  Mgamma <- as.matrix(Mgamma)
  
  ## convert the column names of Mgamma
  colnames(Mgamma) <- unlist( lapply(strsplit(colnames(Mgamma), "", fixed=TRUE), function(x) {
    paste0(paste0(cyto, swap(x, c("0", "1"), c("-", "+"))), collapse="")
  }))
  
  melted <- melt_( as.matrix(Mgamma) )
  melted <- data.table(melted)
  setnames(melted, c("PTID", "Cytokine", "MeanGamma"))
  melted[, Cytokine2 := gsub("[[:alnum:]]", "", Cytokine)]
  
  cyto_plus <- paste0(cyto, "+")
  
  ## we need to pre-process and compute all possible collapses of the data
  melted[, (cyto) := {
    tmp <- strsplit(Cytokine2, "", fixed=TRUE)
    tmp <- lapply(tmp, function(x) {
      as.numeric( swap(x, c("-", "+"), c("0", "1")) )
    })
    tmp <- as.data.frame( t( as.matrix( as.data.frame(tmp) ) ) )
    setNames( tmp, cyto )
  }]
  
  combos <- combinations(cyto)
  marginals <- vector("list", length(combos))
  
  for (i in seq_along(combos)) {
    
    message("Combo ", i, " of ", length(combos), ".")
    combo <- combos[[i]]
    
    ## first, generate the tables
    m <- melted[, list(MeanGamma=mean(MeanGamma)), 
      by=eval(paste("PTID", paste(combo, collapse=","), sep=","))
    ]
    
    ## now, add the cytokine column
    m[, Cytokine := {
      tmp <- character( nrow(m) )
      for (i in 1:length(tmp)) {
        tmp[[i]] <- paste( sep="", collapse="",
          names(.SD), 
          swap(
            unlist(.SD[i]), 
            c("0", "1"), 
            c("-", "+")
          ) )
      }
      list(tmp)
    }, .SDcols=combo]
    
    del(m, list=combo)
    
    ## assign it to the list
    marginals[[i]] <- m
  }
  
  m <- rbindlist(marginals)
  data$Cytokine <- as.character(data$Cytokine)
  
  ## do some incredibly terrible things to ensure the two datasets
  ## can be merged correctly
  data[, Barf := {
    cytos <- strsplit(Cytokine, "[-\\+]")
    signs <- strsplit( gsub("[[:alnum:]]", "", Cytokine), "", fixed=TRUE )
    orders <- lapply(cytos, order)
    
    stopifnot( identical( length(cytos), length(orders) ) )
    stopifnot( identical( length(cytos), length(signs) ) )
    
    garbage <- vector("list", length(cytos))
    for (i in seq_along(garbage)) {
      if (i %% 10000 == 0) message("Iteration ", i, " of ", length(garbage), ".")
      garbage[[i]] <- paste( cytos[[i]][ orders[[i]] ], signs[[i]][ orders[[i]] ], sep='' )
    }
    unlist( lapply( garbage, function(x) paste(x, collapse="") ) )
  }]
  
  m[, Barf := {
    cytos <- strsplit(Cytokine, "[-\\+]")
    signs <- strsplit( gsub("[[:alnum:]]", "", Cytokine), "", fixed=TRUE )
    orders <- lapply(cytos, order)
    
    stopifnot( identical( length(cytos), length(orders) ) )
    stopifnot( identical( length(cytos), length(signs) ) )
    
    garbage <- vector("list", length(cytos))
    for (i in seq_along(garbage)) {
      if (i %% 10000 == 0) message("Iteration ", i, " of ", length(garbage), ".")
      garbage[[i]] <- paste( cytos[[i]][ orders[[i]] ], signs[[i]][ orders[[i]] ], sep='' )
    }
    unlist( lapply( garbage, function(x) paste(x, collapse="") ) )
  }]
  
  m[, Cytokine := Barf]
  data[, Cytokine := Barf]
  
  del(m, Barf)
  del(data, Barf)
  
  setkeyv(m, c("PTID", "Cytokine"))
  setkeyv(data, c("PTID", "Cytokine"))
  
  stopifnot( identical(
    sort(unique(m$Cytokine)),
    sort(unique(data$Cytokine))
  ) )
  
  data <- m[data]
  
  saveRDS(data, file="data/RV144 CaseControl/cd4/cd4_data.rds")
  
}

