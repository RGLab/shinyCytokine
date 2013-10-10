library(hexbin)
library(scales)
library(gridExtra)
library(ggplot2)
library(Kmisc) ## devtools::install_github("Kmisc", "kevinushey")
library(shiny)
library(reshape2)
library(data.table)
library(data.table.extras) ## install_github("data.table.extras", "kevinushey")
library(gtools)
library(stringr)
#library(rCharts)

CELL_DATA <- "data/RV144 CaseControl/res.rds"
COMBO_DATA <- "data/RV144 CaseControl/cd4/cd4_cytokine_combos_preprocessed.rds"
COUNTS <- "data/RV144 CaseControl/cd4_counts.rds"
DATA <- "data/RV144 CaseControl/cd4/cd4_sample_proportions_bg_subtracted.rds"
META <- "data/RV144 CaseControl/meta.rds"

source("common_functions.R")

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
    `LogFoldChange`="Log Fold Change",
    `PropTotal`="Proportion (rel. Total)",
    `PropActivated`="Proportion (rel. Total Activated)",
    `PropTotalBG`="Proportion (rel. Total, BG Corrected)",
    `PropActivatedBG`="Proportion (rel. Total Activated, BG Corrected)"
  ) )
}

## A function used to filter a dataset based on a vector of
## cytokines
filter_cytokines <- function(d, cytokines, cytokines_to_exclude, 
  p, order, phenotype, max_combos) {
  
  ## first, only keep cytokines that have been selected
  for (cytokine in cytokines) {
    d <- droplevels(d[ fgrep(cytokine, Cytokine), ])
  }
  
  ## then, remove cytokines we wish to explicitly exclude
  for (cytokine in cytokines_to_exclude) {
    d <- droplevels(d[ ngrep(cytokine, Cytokine), ])
  }
  
  ## then, remove cytokine combos that have an overall proportion < p
  props <- sort( decreasing=TRUE,
    with(d, tapply_(d[[phenotype]], Cytokine, function(x) median(x, na.rm=TRUE)))
  )
  
  keep <- names(props)[props >= p]
  
  ## keep only terms of a certain order
  keep_order <- str_count(keep, "[-\\+]")
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
    return( input$cytokine_order_min )
  })
  
  getCytokineOrderMax <- reactive({
    return( input$cytokine_order_max )
  })
  
  getCytokineOrder <- reactive({
    return( c(input$cytokine_order_min, input$cytokine_order_max) )
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
  
  getCytokinesToExclude <- reactive({
    return( input$cytokines_to_exclude )
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
  
  # wrap all the plots in an observer -- this is done so that the data is
  ## updated / filtered only once per new parameter update, rather than
  ## once for each plot generated
  
  ## heatmap plot
  output$heatmap <- renderPlot({
    
    facet1 <- getFacet1()
    facet2 <- getFacet2()
    heatmap_level <- getHeatmapLevel()
    cytokines <- getCytokines()
    cytokine_filter <- getCytokineFilter()
    cytokine_order <- getCytokineOrder()
    phenotype <- getPhenotype()
    filter1 <- getFilter1()
    filter1_cb <- getFilter1Cb()
    custom_filter <- getCustomFilter()
    flip_heatmap <- getFlipHeatmap()
    max_combos_to_show <- getMaxCombosToShow()
    cytokines_to_exclude <- getCytokinesToExclude()
    
    ## filter the cytokines
    d <- filter_cytokines(
      d, 
      cytokines, 
      cytokines_to_exclude,
      cytokine_filter, 
      cytokine_order, 
      phenotype, 
      max_combos_to_show
    )
    
    d <- customFilter(d, custom_filter)
    d <- filter1(d, filter1, filter1_cb)
    
    ## the actual heatmap
    p1 <- ggplot( d, aes_string(x="Cytokine", y="Sample", fill=phenotype)) +
      geom_raster() +
      coord_flip() +
      scale_fill_gradient2(
        low="steelblue",
        mid="grey97",
        high="darkorange",
        name=paste(strwrap(phenoToLabel(phenotype), 10), collapse="\n"),
        na.value="maroon"
      ) +
      scale_x_discrete('', expand=c(0,0), breaks=NULL) +
      scale_y_discrete('', expand=c(0,0)) +
      xlab('') +
      ylab(heatmap_level) +
      theme(
        axis.text.x=element_text(color="white"),
        axis.ticks.x=element_line(color="white"),
        title=element_text(family="Gill Sans"),
        plot.margin=unit( c(0, 0, 0.9, 0), "cm" )
      )
    
    ## the color codes
    d_sub <- d[, c("Cytokine", orig_cytokines), with=FALSE]
    d_sub <- data.table(melt_(d_sub, id.vars="Cytokine"))
    setnames(d_sub, c("Combination", "Cytokine", "Used"))
    d_sub$Cytokine <- factor_(d_sub$Cytokine, levels=orig_cytokines)
    d_sub$Used <- factor(d_sub$Used, levels=c(0, 1, 2), labels=c("NA", "Off", "On"))
    
    ## we repeat 'd_sub' to be the number of unique levels generated by
    ## facet 1 with the data set
    if (facet1 != "Original Ordering") {
      n <- length( unique( d[[facet1]] ) )
    } else {
      n <- 1
    }
    
    old_nrow <- nrow(d_sub)
    d_sub <- rbindlist( replicate(n, d_sub, simplify=FALSE) )
    d_sub[, Facet := rep(1:n, each=old_nrow)]
    
    p2 <- ggplot(d_sub, aes(x=Cytokine, y=Combination, fill=Used)) +
      geom_raster() +
      scale_x_discrete('', expand=c(0, 0)) +
      scale_y_discrete('', breaks=NULL, expand=c(0, 0)) +
      scale_fill_manual(
        values=c(
          `NA`="grey95", 
          `Off`="blue",
          `On`="orange"
        )
      ) +
      theme( 
        legend.position='none',
        strip.background=element_blank(),
        strip.text=element_blank(),
        plot.margin=unit( c(0, 0, 0, 0), "cm" ),
        axis.text.x=element_text(angle=90)
      ) +
      facet_grid(Facet ~ .)
    
    ## reorder by phenotype
    if (facet1 != "Original Ordering") {
      
      ## double facetting doesn't work for the heatmap (yet? it might be
      ## too cluttered / ugly anyhow)
      #       if( facet2 == "None" ) {
      p1 <- p1 + facet_grid( paste( facet1, "~ ." ), space="free", scales="free" )
      #       } else {
      #        p1 <- p1 + facet_grid( paste( facet1, "~", facet2 ), space="free", scales="free" )
      #       }
      
    }
    
    grid.arrange( p2, p1, ncol=2, widths=c(15, 85))
    
  })
  
  ## linechart
  output$linechart <- renderPlot({
    
    individual <- getIndividual()
    facet1 <- getFacet1()
    facet2 <- getFacet2()
    heatmap_level <- getHeatmapLevel()
    cytokines <- getCytokines()
    cytokine_filter <- getCytokineFilter()
    cytokine_order <- getCytokineOrder()
    phenotype <- getPhenotype()
    filter1 <- getFilter1()
    filter1_cb <- getFilter1Cb()
    custom_filter <- getCustomFilter()
    flip_linechart <- getFlipLinechart()
    max_combos_to_show <- getMaxCombosToShow()
    cytokines_to_exclude <- getCytokinesToExclude()
    
    ## filter the cytokines
    d <- filter_cytokines(
      d, 
      cytokines, 
      cytokines_to_exclude,
      cytokine_filter, 
      cytokine_order, 
      phenotype, 
      max_combos_to_show
    )
    
    d <- customFilter(d, custom_filter)
    d <- filter1(d, filter1, filter1_cb)
    
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
    cytokines <- getCytokines()
    cytokine_filter <- getCytokineFilter()
    cytokine_order <- getCytokineOrder()
    phenotype <- getPhenotype()
    filter1 <- getFilter1()
    filter1_cb <- getFilter1Cb()
    custom_filter <- getCustomFilter()
    flip_dofplot <- input$flip_dofplot
    
    ## filter the cytokines
    d <- customFilter(d, custom_filter)
    d <- filter1(d, filter1, filter1_cb)
    d <- d[ d$PTID == individual, ]
    
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
    cytokine_order <- getCytokineOrder()
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
    cytokines_to_exclude <- getCytokinesToExclude()
    plot_type <- getPlotType()
    
    ## filter the cytokines
    d <- filter_cytokines(
      d, 
      cytokines, 
      cytokines_to_exclude,
      cytokine_filter, 
      cytokine_order, 
      phenotype, 
      max_combos_to_show
    )
    
    d <- customFilter(d, custom_filter)
    d <- filter1(d, filter1, filter1_cb)
    
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
    cytokine_order <- getCytokineOrder()
    phenotype <- getPhenotype()
    filter1 <- getFilter1()
    filter1_cb <- getFilter1Cb()
    custom_filter <- getCustomFilter()
    cytokines_to_exclude <- getCytokinesToExclude()
    max_combos_to_show <- getMaxCombosToShow()
    
    ## filter the cytokines
    d <- filter_cytokines(
      d, 
      cytokines,
      cytokines_to_exclude,
      cytokine_filter, 
      cytokine_order, 
      phenotype, 
      max_combos_to_show
    )
    
    d <- customFilter(d, custom_filter)
    d <- filter1(d, filter1, filter1_cb)
    d <- droplevels(d)
    
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
    cytokine_order <- getCytokineOrder()
    phenotype <- getPhenotype()
    filter1 <- getFilter1()
    filter1_cb <- getFilter1Cb()
    custom_filter <- getCustomFilter()
    cytokines_to_exclude <- getCytokinesToExclude()
    max_combos_to_show <- getMaxCombosToShow()
    
    ## filter the cytokines
    d <- filter_cytokines(
      d, 
      cytokines,
      cytokines_to_exclude,
      cytokine_filter, 
      cytokine_order, 
      phenotype, 
      max_combos_to_show
    )
    
    d <- customFilter(d, custom_filter)
    d <- filter1(d, filter1, filter1_cb)
    d <- droplevels(d)
    
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
    cytokine_order <- getCytokineOrder()
    phenotype <- getPhenotype()
    filter1 <- getFilter1()
    filter1_cb <- getFilter1Cb()
    custom_filter <- getCustomFilter()
    cytokines_to_exclude <- getCytokinesToExclude()
    max_combos_to_show <- getMaxCombosToShow()
    
    ## filter the cytokines
    d <- filter_cytokines(
      d, 
      cytokines,
      cytokines_to_exclude,
      cytokine_filter, 
      cytokine_order, 
      phenotype, 
      max_combos_to_show
    )
    
    d <- customFilter(d, custom_filter)
    d <- filter1(d, filter1, filter1_cb)
    d <- droplevels(d)
    
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
    cytokines <- getCytokines()
    cytokine_filter <- getCytokineFilter()
    cytokine_order <- getCytokineOrder()
    phenotype <- getPhenotype()
    filter1 <- getFilter1()
    filter1_cb <- getFilter1Cb()
    custom_filter <- getCustomFilter()
    cytokines_to_exclude <- getCytokinesToExclude()
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
