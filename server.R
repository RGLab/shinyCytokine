library(hexbin)
library(scales)
library(gridExtra)
library(ggplot2)
library(Kmisc) ## devtools::install_github("Kmisc", "kevinushey")
library(shiny)
library(reshape2)
library(data.table)
library(data.table.extras) ## devtools::install_github("data.table.extras", "kevinushey")
library(gtools)
library(stringr)
library(pheatmap) ## devtools::install_github("pheatmap", "RGLab", "polyfunctionality")
library(rCharts)

CELL_DATA <- "data/RV144 CaseControl/res.rds"
COMBO_DATA <- "data/RV144 CaseControl/cd4/cd4_cytokine_combos_preprocessed.rds"
COUNTS <- "data/RV144 CaseControl/cd4_counts.rds"
DATA <- "data/RV144 CaseControl/cd4/cd4_data.rds"
META <- "data/RV144 CaseControl/meta.rds"

source("common_functions.R")
source("pheatmap.R")

## Used for the d3 splom
renderSplom <- function(expr, env=parent.frame(), quoted=FALSE) {
  func <- exprToFunction(expr, env, quoted)
  function() {
    val <- func()
    df <- as.data.frame(val[[1]])
    id <- as.character(val[[2]])
    #str(df)
    #str(id)
    id_ind <- match(id, names(df))
    cbind( df[-id_ind], df[id_ind] )
  }
}

## read in the data
## FIXME: only allocate one reference for the data
## need to properly handle subsetting and filtering in the plot commands
counts <- readRDS(COUNTS)
cell_data <- readRDS(CELL_DATA)
combo_data <- readRDS(COMBO_DATA)
d <- dat <- droplevels(readRDS(DATA))
l <- levels( d$Cytokine )
orig_cytokines <- readRDS("data/RV144 CaseControl/cytokines.rds")
meta <- readRDS(META)
counts <- readRDS(COUNTS)
cell_counts <- data.frame(
  Sample=names(counts),
  counts=counts
)

## FIXME: remove sebctrl, negctrl?
d <- dat <- dat[ Stim %nin% c("sebctrl", "negctrl 1", "negctrl 2"), ]
# data_marginal <- data_marginal[ data_marginal$Stim %nin% c("sebctrl", "negctrl 1", "negctrl 2"), ]

## translates phenotype to label for different plots
phenoToLabel <- function(x) {
  return( switch(x,
    `MeanGamma`="Posterior Probability of Expression",
    `LogFoldChange`="Log Fold Change",
    `PropTotal`="Proportion (rel. Total)",
    `PropActivated`="Proportion (rel. Total Activated)",
    `PropTotalBG`="Proportion (rel. Total, BG Corrected)",
    `PropActivatedBG`="Proportion (rel. Total Activated, BG Corrected)"
  ) )
}

## A function used to filter a dataset based on a vector of
## cytokines
filter_cytokines <- function(d, cytokines, 
  cytokines_to_marginalize_over, 
  p, order, phenotype, max_combos) {
  
  ## first, only keep cytokines that have been selected
  for (cytokine in cytokines) {
    d <- droplevels(d[ fgrep(cytokine, Cytokine), ])
  }
  
  ## then, only keep marginals that have been explicitly specified
  ## everything else must go
  remove <- orig_cytokines %knin% cytokines_to_marginalize_over
  get_out <- as.logical(Reduce(`+`, lapply(remove, function(x) {
    return( d[[x]] == 0 )
  })))
  d <- d[!get_out, ]
  
  ## then, remove cytokine combos that have an overall proportion < p
  props <- sort( decreasing=TRUE,
    with(d, tapply_(d[[phenotype]], Cytokine, function(x) median(x, na.rm=TRUE)))
  )
  
  keep <- names(props)[props >= p]
  
  ## keep only terms of a certain order
  keep_order <- str_count(keep, "\\+")
  keep <- keep[ keep_order >= order[1] & keep_order <= order[2] ]
  
  ## and only keep the top 'n' of these combinations
  if (!is.na(max_combos) && max_combos > 0 && length(keep) > max_combos) {
    d <- d[ Cytokine %in% keep[1:max_combos], ]
  } else {
    d <- d[ Cytokine %in% keep, ]
  }
  
  ## and finally, return
  return(d)
}

## a function to filter the data, as based on some custom input from
## the user
customFilter <- function(dat, expr) {
  if( missing(expr) || expr == '' ) {
    return(dat)
  } else {
    return( dat[ eval( parse( text=expr ), envir=dat ), ] )
  }
}

## filter a function based on levels of a factor
filter1 <- function(dat, var, levels) {
  if( var == "None" ) {
    return(dat)
  } else {
    return( dat[ dat[[var]] %in% levels, ] )
  }
}

shinyServer( function(input, output, session) {
  
  ## global variables for the data state
  ## NULL variables are updated by other functions
  d <- dat
  m <- NULL
  dof <- NULL
  row_order <- NULL
  col_order <- NULL
  
  ## define getters
  
  getPhenotype <- reactive({
    return( input$phenotype )
  })
  
  getFacet1 <- reactive({
    return( input$facet1 )
  })
  
  getFacet2 <- reactive({
    return( input$facet2 )
  })
  
  getFacet3 <- reactive({
    return( input$facet3 )
  })
  
  getSample <- reactive({
    return( input$sample )
  })
  
  #   getCytokine <- reactive({
  #     return( input$cytokine )
  #   })
  
  getCytokineFilter <- reactive({
    return( input$cytokine_filter )
  })
  
  getProportionFilter <- reactive({
    return( input$proportion_filter )
  })
  
  getSamplesToKeep <- reactive({
    return( names(overall_sample_prop)[ overall_sample_prop > input$sample_proportion_filter ] )
  })
  
  getPlotType <- reactive({
    return( input$plot_type )
  })
  
  getIndividual <- reactive({
    return( input$individual )
  })
  
  getHeatmapLevel <- reactive({
    return( input$heatmap_level )
  })
  
  getCytokines <- reactive({
    return( input$cytokines )
  })
  
  getCytokineOrderMin <- reactive({
    return( input$cytokine_dof_min )
  })
  
  getCytokineOrderMax <- reactive({
    return( input$cytokine_dof_max )
  })
  
  getCytokineOrder <- reactive({
    return( c(input$cytokine_dof_min, input$cytokine_dof_max) )
  })
  
  getCustomFilter <- reactive({
    return( input$custom_filter )
  })
  
  getFilter1 <- reactive({
    return( input$filter1 )
  })
  
  getFilter1Cb <- reactive({
    return( input$filter1_cb )
  })
  
  #   isMarginal <- reactive({
  #     return( as.logical(input$marginal) )
  #   })
  
  getBoxplotOrientation <- reactive({
    return( input$boxplot_by_cytokine_orientation )
  })
  
  getFlipHeatmap <- reactive({
    return( input$flip_heatmap )
  })
  
  getFlipLinechart <- reactive({
    return( input$flip_linechart )
  })
  
  getBoxplotCoordFlip <- reactive({
    return( input$boxplot_coord_flip )
  })
  
  getBoxplotManualLimits <- reactive({
    return( input$boxplot_manual_limits )
  })
  
  getBoxplotLowerLimit <- reactive({
    return( input$boxplot_lower_limit )
  })
  
  getBoxplotUpperLimit <- reactive({
    return( input$boxplot_upper_limit )
  })
  
  getMaxCombosToShow <- reactive({
    return( as.integer(input$max_combos_to_show) )
  })
  
  getCytokinesToMarginalizeOver <- reactive({
    return( input$cytokines_to_marginalize_over )
  })
  
  ## observers
  observe({
    x <- input$individual
    samples <- meta$name[ meta$PTID == x ]
    updateSelectInput(session, "sample",
      label="Sample",
      choices=samples
    )
  })
  
  observe({
    x <- input$boxplot_manual_limits
    updateNumericInput(session, "boxplot_lower_limit", value=0)
    updateNumericInput(session, "boxplot_upper_limit", value=0.005)
  })
  
  ## if the phenotype is changed, make sure the 'set limits manually'
  ## thing gets unchecked
  observe({
    x <- input$phenotype
    updateCheckboxInput(session, "boxplot_manual_limits", value=FALSE)
  })
  
  #   observe({
  #     x <- input$heatmap_level
  #     updateSliderInput(session, 
  #       "proportion_filter", 
  #       label=paste("Remove", x, "with proportion < x")
  #     )
  #   })
  
  ## an observer for conditionally updating the filter widget,
  ## depending on the type of variable being used
  observe({
    x <- input$filter1
    if( x != "None" ) {
      m <- d[[x]]
      if( is.factor(m) || is.character(m) ) {
        
        updateCheckboxGroupInput(
          session, 
          "filter1_cb", 
          choices=levels(d[[x]]),
          selected=levels(d[[x]])
        )
        
      } else {
        
        ## TODO: numeric filter
        
      }
    }
  })
  
  ## wrap data updates in an observer
  observe({
    
    ## The variables we need to observe for updating the data
    cytokines <- getCytokines()
    cytokine_filter <- getCytokineFilter()
    cytokine_dof <- getCytokineOrder()
    phenotype <- getPhenotype()
    filter1 <- getFilter1()
    filter1_cb <- getFilter1Cb()
    custom_filter <- getCustomFilter()
    max_combos_to_show <- getMaxCombosToShow()
    cytokines_to_marginalize_over <- getCytokinesToMarginalizeOver()
    
    ## filter the cytokines
    tmp <- filter_cytokines(
      d, 
      cytokines, 
      cytokines_to_marginalize_over,
      cytokine_filter, 
      cytokine_dof, 
      phenotype, 
      max_combos_to_show
    )
    
    tmp <- customFilter(tmp, custom_filter)
    tmp <- filter1(tmp, filter1, filter1_cb)
    
    ## Cast the cytokine data down into a matrix
    m <<- acast(tmp, Sample ~ Cytokine, value.var=phenotype)
    
    ## Rows are samples, columns are cytokine combinations
    
    ## Generate the cytokine / column annotations
    annot <<- data.frame(do.call(cbind, lapply(orig_cytokines, function(x) {
      
      ## the positives
      pos <- as.integer(grepl( paste0(x, "\\+"), colnames(m)))
      
      ## the negatives
      neg <- -as.integer(grepl( paste0(x, "-"), colnames(m)))
      
      ## merge them
      output <- rep(0, length(pos))
      output[pos == 1] <- 1
      output[neg == -1] <- -1
      return(output)
      
    })))
    
    rownames(annot) <- colnames(m)
    colnames(annot) <- orig_cytokines
    
    ## Order the columns by degree of functionality, number of
    ## marginals, and then actual mean phenotype.
    dof <<- apply(annot, 1, function(x) {
      sum(x == 1)
    })
    
    nom <<- apply(annot, 1, function(x) {
      sum(x == 0)
    })
    
    mean_pheno <<- apply(m, 2, function(x) {
      mean(x, na.rm=TRUE)
    })
    
    col_order <<- do.call(order, list(dof, nom, mean_pheno))
    
    d <<- tmp
    
  })
  
  ## heatmap plot
  output$heatmap <- renderPlot({
    
    facet1 <- getFacet1()
    facet2 <- getFacet2()
    facet3 <- getFacet3()
    heatmap_level <- getHeatmapLevel()
    phenotype <- getPhenotype()
    
    ## A quick function used for getting a nice color palette
    heatmap_color <- function(phenotype) {
      switch(phenotype,
        `MeanGamma`=colorRampPalette(brewer.pal(n=9, name="Purples")[2:7])(100),
        colorRampPalette(brewer.pal(n=7, name="PuOr"))(100)
      )
    }
    
    if (facet1 == "Original Ordering") facet1 <- NULL
    if (facet2 == "None") facet2 <- NULL
    if (facet3 == "None") facet3 <- NULL
    
    ## Generate the row annotations; ie, annotate each row (sample) by
    ## facet
    if (any( !is.null( c(facet1, facet2, facet3) ) )) {
      
      row_annot <- d[, c("Sample", facet1, facet2, facet3), with=FALSE]
      row_annot <- as.data.frame(collapse(row_annot, by=Sample))
      row_annot$Sample <- as.character(row_annot$Sample)
      row_annot <- row_annot[ match( rownames(m), row_annot$Sample ), ]
      rownames(row_annot) <- row_annot$Sample
      del(row_annot, Sample)
      
      row_order <- do.call(order, row_annot[ c(facet1, facet2, facet3) ])
      
      pheatmap(m[row_order, col_order, drop=FALSE], 
        color=heatmap_color(phenotype),
        cluster_rows=FALSE,
        cluster_cols=FALSE,
        show_rownames=FALSE,
        show_colnames=FALSE,
        cytokine_annotation=annot[col_order, , drop=FALSE],
        row_annotation=row_annot[row_order, , drop=FALSE],
        main=paste(phenoToLabel(phenotype), "by Sample")
      )
      
    } else {
      
      pheatmap(m[, col_order, drop=FALSE], 
        color=heatmap_color(phenotype),
        cluster_rows=FALSE,
        cluster_cols=FALSE,
        show_rownames=FALSE,
        show_colnames=FALSE,
        cytokine_annotation=annot[col_order, , drop=FALSE],
        main=paste(phenoToLabel(phenotype), "by Sample")
      )
      
    }
    
  })
  
  ## linechart
  output$linechart <- renderPlot({
    
    individual <- getIndividual()
    facet1 <- getFacet1()
    facet2 <- getFacet2()
    heatmap_level <- getHeatmapLevel()
    phenotype <- getPhenotype()
    flip_linechart <- getFlipLinechart()
    
    d_sub <- droplevels( d[ d$PTID == individual, ] )
    p <- ggplot( d_sub, aes_string(x="Cytokine", y=phenotype, group="Sample")) +
      geom_point() +
      geom_line( aes(x=as.integer( factor(Cytokine) ))) +
      xlab(paste("Individual:", individual)) +
      ggtitle(paste(phenoToLabel(phenotype), "by Cytokine")) +
      #guides(color=guide_legend(nrow=10)) +
      theme( 
        axis.text.x = element_text(angle=45, hjust=1),
        title = element_text(family="Gill Sans")
      )
    
    if( facet1 != "Original Ordering" ) {
      p <- p + aes_string(x="Cytokine", y=phenotype, group="Sample", color=facet1)
      if( facet2 != "None" ) {
        p <- p + facet_wrap(facet2)
      }
    }
    
    ## y label -- for some reason, it has to occur after the facetting code
    p <- p +
      ylab( phenoToLabel(phenotype) )
    
    if (flip_linechart) {
      p <- p +
        coord_flip() +
        theme(
          title=element_text(family="Gill Sans")
        )
    }
    
    plot(p)
    
  })
  
  output$dofplot <- renderPlot({
    
    individual <- getIndividual()
    facet1 <- getFacet1()
    facet2 <- getFacet2()
    heatmap_level <- getHeatmapLevel()
    phenotype <- getPhenotype()
    flip_dofplot <- input$flip_dofplot
    
    ## get the samples belonging to the currently selected individual
    samples <- cell_data[ names(cell_data) %in% d$Sample ]
    
    ## computing the degree of functionality
    ## for each cell, we compute the number of cytokines that were
    ## expressed (expression value > 0). then, we generate 
    ## counts of these values. for example, given 3 cells in a sample like
    
    ## IL1 TNFa INFg
    ## 0 1 0
    ## 0 0 1
    ## 1 1 0
    
    ## we should get output like
    
    ## sample num_cytokines_expressed count
    ## 1 0 0
    ## 1 1 2
    ## 1 2 1
    ## 1 3 0
    
    ## the zero counts can be left out
    
    ## compute the degree of functionality:
    ## for each cell, the DOF == the number of cytokines expressed
    samples_dof <- stack( lapply( samples, function(x) {
      rowSums(x > 0)
    }) )
    
    samples_dof <- data.table(samples_dof)
    samples_dof <- samples_dof[, {
      counts <- Kmisc::counts(values)
      list(counts=counts, order=names(counts))
    }, by=ind]
    ## merge in cell counts
    dof <- merge.data.frame( samples_dof, cell_counts, by.x="ind", by.y="Sample", all.x=TRUE )
    dof$prop <- dof$counts.x / dof$counts.y
    
    ## merge in meta-info
    dof <- unique( merge( 
      dof, 
      d[ names(d) %in% c("Sample", names(meta)) ], 
      by.x="ind", 
      by.y="Sample",
      all.x=TRUE
    ) )
    
    p <- ggplot(dof, aes(x=factor(order), y=prop)) +
      theme(
        title=element_text(family="Gill Sans")
      ) +
      xlab("Degree of Functionality") +
      ylab("Proportion of Cells") +
      ggtitle(paste("Degree of Functionality\nIndividual:", individual))
    
    if( facet1 != "Original Ordering" ) {
      if( facet2 == "None" ) {
        p <- p + 
          geom_bar( aes_string(fill=facet1), stat="identity", position="dodge" )
        #facet_grid( paste( facet1, "~ ." ) )
      } else {
        p <- p + 
          geom_bar( aes_string(fill=facet1), stat="identity", position="dodge") +
          facet_grid( paste( facet2, "~ ." ) )
      }
    } else {
      p <- p +
        geom_bar(stat="identity", position="dodge")
      
    }
    
    if (flip_dofplot) {
      p <- p + coord_flip()
    }
    
    plot(p)
    
  })
  
  ## boxplot by cytokine
  output$boxplot_by_cytokine <- renderPlot({
    
    facet1 <- getFacet1()
    facet2 <- getFacet2()
    facet3 <- getFacet3()
    heatmap_level <- getHeatmapLevel()
    cytokines <- getCytokines()
    cytokine_filter <- getCytokineFilter()
    cytokine_dof <- getCytokineOrder()
    phenotype <- getPhenotype()
    filter1 <- getFilter1()
    filter1_cb <- getFilter1Cb()
    custom_filter <- getCustomFilter()
    boxplot_orientation <- getBoxplotOrientation()
    boxplot_coord_flip <- getBoxplotCoordFlip()
    boxplot_manual_limits <- getBoxplotManualLimits()
    boxplot_lower_limit <- getBoxplotLowerLimit()
    boxplot_upper_limit <- getBoxplotUpperLimit()
    max_combos_to_show <- getMaxCombosToShow()
    cytokines_to_marginalize_over <- getCytokinesToMarginalizeOver()
    plot_type <- getPlotType()
    
    if( facet1 == "Original Ordering" ) {
      plot( 1 ~ 1, type='n', axes=FALSE, frame.plot=FALSE, ann=FALSE)
      text( x=1, y=1, labels="Select a Facet to view Boxplots of Proportions by Cytokine")
      return(NULL)
    }
    
    p <- ggplot(d) +
      theme( 
        #axis.text.x = element_text(angle=45, hjust=1),
        title = element_text(family="Gill Sans")
      )
    
    ## what kind of plot to generate?
    if( facet2 == "None" ) {
      p <- p + switch(plot_type,
        boxplot=geom_boxplot( aes_string(x=facet1, y=phenotype), fill="#ABCDEF"),
        histogram=geom_histogram( aes_string(x=phenotype, fill=facet1), alpha=0.5, position="identity" ),
        density=geom_density( aes_string(x=phenotype, fill=facet1), alpha=0.5 )
      )
    } else {
      p <- p + switch(plot_type,
        boxplot=geom_boxplot( aes_string(x=facet2, y=phenotype, fill=facet1)),
        histogram=geom_histogram( aes_string(x=phenotype, fill=facet1), alpha=0.5, position="identity" ),
        density=geom_density( aes_string(x=phenotype, fill=facet1), alpha=0.5 )
      )
    }
    
    if (boxplot_coord_flip) {
      if (facet3 != "None") {
        p <- p +
          ylab( phenoToLabel(phenotype) ) +
          facet_grid(paste("Cytokine ~", facet3))
      } else {
        p <- p +
          ylab( phenoToLabel(phenotype) ) +
          facet_grid("Cytokine ~ .")
      }
    } else {
      if (facet3 != "None") {
        p <- p +
          ylab( phenoToLabel(phenotype) ) +
          facet_grid(paste(facet3, "~ Cytokine"))
      } else {
        p <- p +
          ylab( phenoToLabel(phenotype) ) +
          facet_grid(". ~ Cytokine")
      }
    }
    
    if (boxplot_orientation == "Horizontal") {
      if (boxplot_manual_limits) {
        p <- p +
          coord_flip(ylim=c(boxplot_lower_limit, boxplot_upper_limit))
      } else {
        p <- p +
          coord_flip()
      }
    } else {
      if (boxplot_manual_limits) {
        p <- p +
          coord_cartesian(ylim=c(boxplot_lower_limit, boxplot_upper_limit))
      } 
    }
    
    plot(p)
    
  })
  
  output$stats <- renderChart2({
    
    individual <- getIndividual()
    facet1 <- getFacet1()
    facet2 <- getFacet2()
    facet3 <- getFacet3()
    heatmap_level <- getHeatmapLevel()
    cytokines <- getCytokines()
    cytokine_filter <- getCytokineFilter()
    cytokine_dof <- getCytokineOrder()
    phenotype <- getPhenotype()
    filter1 <- getFilter1()
    filter1_cb <- getFilter1Cb()
    custom_filter <- getCustomFilter()
    cytokines_to_marginalize_over <- getCytokinesToMarginalizeOver()
    max_combos_to_show <- getMaxCombosToShow()
    
    if (nrow(d) == 0) {
      return("Empty data")
    }
    
    ## construct 'by' from the included facets
    by <- "Cytokine"
    
    if (facet1 != "Original Ordering") {
      by <- paste(by, facet1, sep=",")
    }
    
    if (facet2 != "None") {
      by <- paste(by, facet2, sep=",")
    }
    
    if (facet3 != "None") {
      by <- paste(by, facet3, sep=",")
    }
    
    m <- d[, list(
      Mean=mean( get(phenotype), na.rm=TRUE ),
      Median=median( get(phenotype), na.rm=TRUE ),
      SD=sd( get(phenotype), na.rm=TRUE ),
      N=.N
    ), keyby=by]
    
    dTable(m)
    
  })
  
  output$splom <- renderPlot({
    
    individual <- getIndividual()
    facet1 <- getFacet1()
    facet2 <- getFacet2()
    heatmap_level <- getHeatmapLevel()
    cytokines <- getCytokines()
    cytokine_filter <- getCytokineFilter()
    cytokine_dof <- getCytokineOrder()
    phenotype <- getPhenotype()
    filter1 <- getFilter1()
    filter1_cb <- getFilter1Cb()
    custom_filter <- getCustomFilter()
    cytokines_to_marginalize_over <- getCytokinesToMarginalizeOver()
    max_combos_to_show <- getMaxCombosToShow()
    
    m <- combo_data[, colnames(combo_data) %in% d$Cytokine ]
    gp <- factor(rownames(m))
    nm <- colnames(m)
    nm <- gsub("\\+", "\\+\n", nm)
    nm <- gsub("-", "-\n", nm)
    
    m[ m < 200 ] <- NA
    
    p1 <- splom(m,
      pscales=3,
      type=c('g', 'p'),
      varnames=nm,
      varname.cex=0.8,
      panel=function(x, y, i, j, ...) {
        cor <- sprintf("%.3f", cor(x, y, use="pairwise"))
        expr <- as.expression( substitute( paste(
          rho == x
        ), list(x=cor) ) )
        panel.hexplom(x, y, ...)
        panel.loess(x, y, col='red', ...)
        grid.text2(expr, x=0.5, y=0.85)
      }
    )
    
    print(p1)
    
  })
  
  output$splom_cell <- renderPlot({
    
    individual <- getIndividual()
    facet1 <- getFacet1()
    facet2 <- getFacet2()
    heatmap_level <- getHeatmapLevel()
    cytokines <- getCytokines()
    cytokine_filter <- getCytokineFilter()
    cytokine_dof <- getCytokineOrder()
    phenotype <- getPhenotype()
    filter1 <- getFilter1()
    filter1_cb <- getFilter1Cb()
    custom_filter <- getCustomFilter()
    cytokines_to_marginalize_over <- getCytokinesToMarginalizeOver()
    max_combos_to_show <- getMaxCombosToShow()
    
    samples <- unique(d$Sample)
    dat <- cell_data[ names(cell_data) %in% samples ]
    m <- do.call(rbind, dat)
    if (nrow(m) > 1E4) {
      m <- m[1:1E4, ]
    }
    
    p1 <- hexplom(m,
      xbins=8,
      pscales=3,
      type=c('g', 'p'),
      varname.cex=0.8,
      panel=function(x, y, ...) {
        tx <- x[x > 0 & y > 0]
        ty <- y[x > 0 & y > 0]
        if (length(tx) == 0 || length(ty) == 0) {
          panel.text("No data")
        }
        panel.hexbinplot(tx, ty, ...)
      }
    )
    
    print(p1)
    
  })
  
  output$splom <- renderSplom({
    
    individual <- getIndividual()
    facet1 <- getFacet1()
    facet2 <- getFacet2()
    heatmap_level <- getHeatmapLevel()
    phenotype <- getPhenotype()
    filter1 <- getFilter1()
    filter1_cb <- getFilter1Cb()
    custom_filter <- getCustomFilter()
    cytokines_to_marginalize_over <- getCytokinesToMarginalizeOver()
    max_combos_to_show <- getMaxCombosToShow()
    
    ## take the selected individual
    
    samples <- as.character(unique(d$Sample[ d$PTID == individual ]))
    dat <- cell_data[ names(cell_data) %in% samples ]
    dat <- lapply(dat, as.data.frame)
    m <- rbindlistn(dat, TRUE)
    setnames(m, ".Names", "Sample")
    list(m, "Sample")
    
  })
  
})
