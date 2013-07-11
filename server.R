library(ggplot2)
library(Kmisc) ## devtools::install_github("Kmisc", "kevinushey", subdir="Kmisc")
library(shiny)
library(reshape2)
library(data.table)
library(data.table.extras) ## install_github("data.table.extras", "kevinushey")
library(gtools)
library(stringr)
#library(rCharts)

DATA <- "data/cd4/sample_proportions_bg_substracted.rds"
META <- "data/meta.rds"

source("common_functions.R")

## TODO: allow other plot types
plot_type <- "boxplot"

## translates phenotype to label for different plots
phenoToLabel <- function(x) {
  return( switch(x,
    `Proportion`="Proportion",
    `Proportion_bg`="Proportion\n(BG Corrected)"
  ) )
}

shinyServer( function(input, output, session) {
  
  y_var <- "Sample"
  
  ## A function used to filter a dataset based on a vector of
  ## cytokines
  filter_cytokines <- function(d, cytokines, p, order) {
    
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
    props <- with(d, tapply_(Proportion, Cytokine, mean))
    keep <- names(props)[props >= p]
    
    ## keep only terms of a certain order
    keep_order <- str_count(keep, " x ") + 1
    keep <- keep[ keep_order >= order[1] ]
    keep <- keep[ keep_order <= order[2] ]
    
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
  
  ## read in the data
  ## FIXME: only allocate one reference for the data
  ## need to properly handle subsetting and filtering in the plot commands
  d <- dat <- readRDS(DATA)
  l <- levels( d$Cytokine )
  orig_cytokines <- strsplit( l[ length(l) ], " x ", fixed=TRUE)[[1]]
  meta <- readRDS(META)
  
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
    x <- input$heatmap_level
    updateSliderInput(session, 
      "proportion_filter", 
      label=paste("Remove", x, "with proportion < x")
    )
  })
  
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
  
  ## wrap all the plots in an observer -- this is done so that the data is
  ## updated / filtered only once per new parameter update, rather than
  ## once for each plot generated
  
  ## FIXME: plots no longer grey out when being updated inside an `observe`
  ## call
  
  ## heatmap plot
  output$heatmap <- renderPlot({
    
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
    
    ## filter the cytokines
    d <- filter_cytokines(dat, cytokines, cytokine_filter, cytokine_order)
    d <- customFilter(d, custom_filter)
    d <- filter1(d, filter1, filter1_cb)
    
    p1 <- ggplot( d, aes_string(x="Cytokine", y=y_var, fill=phenotype)) +
      geom_tile() +
      scale_fill_gradient(low="white", high="steelblue", name=phenoToLabel(phenotype)) +
      scale_x_discrete(expand=c(0,0)) +
      scale_y_discrete(expand=c(0,0), breaks=NULL) +
      xlab("Cytokine") +
      ylab(heatmap_level) +
      theme( 
        axis.text.x = element_text(angle=45, hjust=1),
        title = element_text(family="Gill Sans")
      ) +
      ggtitle("Proportion of Cells Expressing Cytokines by Sample")
    
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
    
    ## filter the cytokines
    d <- filter_cytokines(dat, cytokines, cytokine_filter, cytokine_order)
    d <- customFilter(d, custom_filter)
    d <- filter1(d, filter1, filter1_cb)
    
    d_sub <- droplevels( d[ d$PTID == individual, ] )
    p <- ggplot( d_sub, aes_string(x="Cytokine", y=phenotype, group="Sample")) +
      geom_point() +
      geom_line( aes(x=as.integer( factor(Cytokine) ))) +
      xlab(paste("Individual:", individual)) +
      ggtitle("Proportion by Cytokine") +
      guides(color=guide_legend(nrow=10)) +
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
    
    ## filter the cytokines
    d <- filter_cytokines(dat, cytokines, cytokine_filter, cytokine_order)
    d <- customFilter(d, custom_filter)
    d <- filter1(d, filter1, filter1_cb)
    
    p <- ggplot(d) +
      theme( 
        axis.text.x = element_text(angle=45, hjust=1),
        title = element_text(family="Gill Sans")
      )
    
    if( facet2 != "None" ) {
      p <- p + facet_grid( paste(". ~", facet2) )
    }
    
    plot(p)
    
  })
  
  ## boxplot by cytokine
  output$boxplot_by_cytokine <- renderPlot({
    
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
    
    ## filter the cytokines
    d <- filter_cytokines(dat, cytokines, cytokine_filter, cytokine_order)
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
        boxplot=geom_boxplot( aes_string(x=facet1, y=phenotype) ),
        histogram=geom_histogram( aes_string(x=phenotype, fill=facet1), alpha=0.5, position="identity" ),
        density=geom_density( aes_string(x=phenotype, fill=facet1), alpha=0.5 )
      )
    } else {
      p <- p + switch(plot_type,
        boxplot=geom_boxplot( aes_string(x=facet2, y=phenotype, fill=facet1) ),
        histogram=geom_histogram( aes_string(x=phenotype, fill=facet1), alpha=0.5, position="identity" ),
        density=geom_density( aes_string(x=phenotype, fill=facet1), alpha=0.5 )
      )
    }
    
    if( plot_type == "boxplot" ) {
      p <- p +
        ylab( phenoToLabel(phenotype) ) +
        facet_grid(". ~ Cytokine") +
        coord_flip()
    } else {
      p <- p +
        xlab( phenoToLabel(phenotype) ) +
        facet_grid( paste(facet2, "~", "Cytokine") )
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
    
    ## filter the cytokines
    d <- filter_cytokines(dat, cytokines, cytokine_filter, cytokine_order)
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