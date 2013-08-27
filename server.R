library(RColorBrewer)
library(ggplot2)
library(Kmisc) ## devtools::install_github("Kmisc", "kevinushey", subdir="Kmisc")
library(shiny)
library(reshape2)
library(data.table)
library(data.table.extras) ## install_github("data.table.extras", "kevinushey")
library(gtools)
library(stringr)
library(scales)
#library(rCharts)

CELL_DATA <- "data/RV144 CaseControl/res.rds"
COUNTS <- "data/RV144 CaseControl/cd4_counts.rds"
DATA <- "data/RV144 CaseControl/cd4_joint/cd4_joint_sample_proportions_bg_subtracted.rds"
DATA_MARGINAL <- "data/RV144 CaseControl/cd4_marginal/cd4_marginal_sample_proportions_bg_subtracted.rds"
META <- "data/RV144 CaseControl/meta.rds"

source("common_functions.R")

## read in the data
## FIXME: only allocate one reference for the data
## need to properly handle subsetting and filtering in the plot commands
cell_data <- readRDS(CELL_DATA)
d <- dat <- readRDS(DATA)
data_marginal <- readRDS(DATA_MARGINAL)
l <- levels( d$Cytokine )
orig_cytokines <- strsplit( l[ length(l) ], " x ", fixed=TRUE)[[1]]
meta <- readRDS(META)
counts <- readRDS(COUNTS)
cell_counts <- data.frame(
  Sample=names(counts),
  counts=counts
)

## FIXME: remove sebctrl, negctrl?
d <- dat 
# d <- dat <- dat[ dat$Stim %nin% c("sebctrl", "negctrl 1", "negctrl 2"), ]
# data_marginal <- data_marginal[ data_marginal$Stim %nin% c("sebctrl", "negctrl 1", "negctrl 2"), ]

## TODO: allow other plot types
plot_type <- "boxplot"

## translates phenotype to label for different plots
phenoToLabel <- function(x) {
  return( switch(x,
    `LogFoldChange`="Log Fold Change",
    `Proportion`="Proportion",
    `Proportion_bg`="Proportion\n(BG Corrected)",
    `Proportion_log10`="Proportion\n(log10 transformed)",
    `Proportion_bg_log10`="Proportion\n(BG corrected, log10 transformed)"
  ) )
}

y_var <- "Sample"

## A function used to filter a dataset based on a vector of
## cytokines
filter_cytokines <- function(d, cytokines, p, order, phenotype) {
  
  ## first, only keep cytokines that have been selected
  if( length(cytokines) == 0 ) {
    d <- d[ d$Cytokine %in% orig_cytokines, ]
  } else {
    n <- length(cytokines)
    combos <- unlist( combinations(cytokines, function(x) {
      paste(x, collapse=" x ")
    }) )
    d <- d[ d$Cytokine %in% combos, ]
  }
  
  ## then, remove cytokine combos that have an overall proportion < p
  props <- with(d, tapply_(d[[phenotype]], Cytokine, function(x) mean(x, na.rm=TRUE)))
  keep <- names(props)[props >= p]
  
  ## keep only terms of a certain order
  keep_order <- str_count(keep, " x ") + 1
  keep <- keep[ keep_order >= order[1] & keep_order <= order[2] ]
  
  d <- d[ d$Cytokine %in% keep, ]
  
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
  
  isMarginal <- reactive({
    return( as.logical(input$marginal) )
  })
  
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
  
  ## FIXME: plots no longer grey out when being updated inside an `observe`
  ## call
  
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
    marginal <- isMarginal()
    
    if (marginal) {
      d <- data_marginal
    } else {
      d <- dat
    }
    
    ## filter the cytokines
    d <- filter_cytokines(d, cytokines, cytokine_filter, cytokine_order, phenotype)
    d <- customFilter(d, custom_filter)
    d <- filter1(d, filter1, filter1_cb)
    
    p1 <- ggplot( d, aes_string(x="Cytokine", y=y_var, fill=phenotype)) +
      geom_tile() +
      scale_fill_gradient2(
        low="steelblue",
        mid="white",
        high="darkorange",
        name=phenoToLabel(phenotype)
      ) +
      scale_x_discrete(expand=c(0,0)) +
      scale_y_discrete(expand=c(0,0), breaks=NULL) +
      xlab("Cytokine") +
      ylab(heatmap_level) +
      theme( 
        axis.text.x = element_text(angle=45, hjust=1),
        title = element_text(family="Gill Sans")
      )
    
    ## reorder by phenotype
    if( facet1 != "Original Ordering" ) {
      
      ## double facetting doesn't work for the heatmap (yet? it might be
      ## too cluttered / ugly anyhow)
      #       if( facet2 == "None" ) {
      p1 <- p1 + facet_grid( paste( facet1, "~ ." ), space="free", scales="free" )
      #       } else {
      #        p1 <- p1 + facet_grid( paste( facet1, "~", facet2 ), space="free", scales="free" )
      #       }
      
    }
    
    ## flip the axes?
    if (flip_heatmap) {
      p1 <- p1 +
        coord_flip() +
        theme(
          title=element_text(family="Gill Sans")
        )
    }
    
    plot(p1)
    
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
    marginal <- isMarginal()
    flip_linechart <- getFlipLinechart()
    
    if(marginal) {
      d <- data_marginal
    } else {
      d <- dat
    }
    
    ## filter the cytokines
    d <- filter_cytokines(d, cytokines, cytokine_filter, cytokine_order, phenotype)
    d <- customFilter(d, custom_filter)
    d <- filter1(d, filter1, filter1_cb)
    
    d_sub <- droplevels( d[ d$PTID == individual, ] )
    p <- ggplot( d_sub, aes_string(x="Cytokine", y=phenotype, group="Sample")) +
      geom_point() +
      geom_line( aes(x=as.integer( factor(Cytokine) ))) +
      xlab(paste("Individual:", individual)) +
      ggtitle("Proportion by Cytokine") +
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
    marginal <- isMarginal()
    flip_dofplot <- input$flip_dofplot
    
    if(marginal) {
      d <- data_marginal
    } else {
      d <- dat
    }
    
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
      ggtitle(paste("DOF Plot\nIndividual:", individual))
    
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
    marginal <- isMarginal()
    boxplot_orientation <- getBoxplotOrientation()
    boxplot_coord_flip <- getBoxplotCoordFlip()
    boxplot_manual_limits <- getBoxplotManualLimits()
    boxplot_lower_limit <- getBoxplotLowerLimit()
    boxplot_upper_limit <- getBoxplotUpperLimit()
    
    if(marginal) {
      d <- data_marginal
    } else {
      d <- dat
    }
    
    ## filter the cytokines
    d <- filter_cytokines(d, cytokines, cytokine_filter, cytokine_order, phenotype)
    d <- customFilter(d, custom_filter)
    d <- filter1(d, filter1, filter1_cb)
    
    if( facet1 == "Original Ordering" ) {
      plot( 1 ~ 1, type='n', axes=FALSE, frame.plot=FALSE, ann=FALSE)
      text( x=1, y=1, labels="Select a Facet to view Boxplots of Proportions by Cytokine")
      return(NULL)
    }
    
    p <- ggplot(d) +
      theme( 
        axis.text.x = element_text(angle=45, hjust=1),
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
    
    if( plot_type == "boxplot" ) {
      
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
      
    }
    
    
    
    plot(p)
    
  })
  
  output$stats <- renderTable({
    
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
    marginal <- isMarginal()
    
    if(marginal) {
      d <- data_marginal
    } else {
      d <- dat
    }
    
    ## filter the cytokines
    d <- filter_cytokines(d, cytokines, cytokine_filter, cytokine_order, phenotype)
    d <- customFilter(d, custom_filter)
    d <- filter1(d, filter1, filter1_cb)
    d <- droplevels(d)
    
    if( facet1 == "Original Ordering" ) {
      split_var <- d$Cytokine
    } else {
      split_var <- paste( d$Cytokine, d[[facet1]], sep=' x ' )
    }
    
    m <- do.call( rbind, lapply( split(d, split_var), function(DF) {
      return( summary( DF[[phenotype]] ) )
    } ) )
    
    return( as.data.frame(m) )
    
  })
  
  #   output$rchart <- renderChart({
  #     
  #     p <- rPlot( "Cytokine", phenotype, data=d, color=facet1, facet=facet2, type='point' )
  #     p$addParams(dom="rchart", width=430, height=300)
  #     return(p)
  #     
  #   })
  
})