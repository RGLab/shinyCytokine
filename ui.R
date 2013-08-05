#library(rCharts)
library(shinyGridster)

meta <- readRDS("data/meta.rds")
dat <- readRDS("data/res_cd4.rds")

## width, height for gridster + plot elements
width <- 430
height <- 300

## svg output
svgOutput <- function(outputId, width, height) {
  div(
    tag("svg", list(id=outputId, width=width, height=height, class="html-shiny-output"))
  )
}

## ensure that each matrix has the same column names
stopifnot( length( table( table( unlist( lapply( dat, names ) ) ) ) ) != 1 )

shinyUI( bootstrapPage(
  
  includeCSS("www/css/styles.css"),
  includeScript("www/js/fancyboxify.js"),
  
  singleton( tags$body( style="background-color: #789;" ) ),
  
  h1(style="text-align: center; color: white", "Cytokine Visualization"),
  
  gridster( width=width, height=height,
    
    gridsterItem( row=1, col=1, sizex=1, sizey=2,
      
      tags$div( style="overflow: auto; width: 100%;",
        
        tags$div( style="float: left; width: 45%;",
          selectInput("phenotype", label="Phenotype", choices=list(
            `Proportion`="Proportion",
            `Proportion (BG Corrected)`="Proportion_bg"
          ))
        ),
        
        tags$div( style="float: left; margin-left: 10px; width: 45%;",
          selectInput("marginal", label="Distribution", choices=list(
            Marginal=TRUE,
            Joint=FALSE
          ))
        )
        
      ),
      
      #       selectInput("sample",
      #         label="Sample",
      #         choices=names(dat)
      #       ),
      
      #       selectInput("cytokine",
      #         label="Cytokine",
      #         choices=unname( colnames(dat[[1]]) )
      #       ),
      
      #       h3("Filters"),
      #       sliderInput("proportion_filter",
      #         label="Remove Individuals with Overall Proportion < x",
      #         min=0,
      #         max=1,
      #         value=0.05,
      #         step=0.01
      #       ),
      
      ## multiple selectInput doesn't work with gridster :( (yet?)
      #h3("Cytokine Combinations"),
      checkboxGroupInput("cytokines",
        label="Cytokine Combinations",
        choices=unname( colnames( dat[[1]] ) ),
        selected=unname( colnames( dat[[1]] ) )
      ),
      
      sliderInput("cytokine_filter",
        label="Remove Cytokine Combinations with p < x",
        min=0,
        max=0.2,
        value=0.05,
        step=0.001
      ),
      
      ## overflow: auto keeps div from collapsing to zero height
      ## see: http://stackoverflow.com/questions/218760/how-do-you-keep-parents-of-floated-elements-from-collapsing
      tags$div( style="overflow: auto;",
        tags$div( style="width: 45%; float: left;",
          tags$label( `for`="cytokine_order_min", "Minimum Cytokine Order Combination"),
          tags$input( style="width: 80%;", id="cytokine_order_min", type="number", value="1", min="1", max=ncol( dat[[1]] ), step="1" )
        ),
        tags$div( style="width: 45%; float: right;", 
          tags$label( `for`="cytokine_order_max", "Maximum Cytokine Order Combination"),
          tags$input( style="width: 80%;", id="cytokine_order_max", type="number", value=ncol( dat[[1]] ), min="1", max=ncol( dat[[1]] ), step="1" )
        )
      ),
      
      #h3("Facets"),
      tags$div( style="overflow: auto;",
        
        tags$div( style="width: 45%; float: left;",
          selectInput("facet1",
            label="Facet 1",
            choices=c("Original Ordering", names(meta))
          )
        ),
        
        tags$div( style="width: 45%; float: right;",
          selectInput("facet2",
            label="Facet 2",
            choices=c("None", names(meta))
          )
        )
        
      ),
      
      #h3("Plot Type"),
      #       selectInput("plot_type",
      #         label="Plot Type",
      #         choices=list(
      #           `Boxplots`="boxplot",
      #           `Histograms`="histogram",
      #           `Density Plots`="density"
      #         )
      #       ),
      
      selectInput("filter1",
        label="Filter 1",
        choices=c("None", names(meta))
      ),
      
      ## this panel will be updated by server.R -- displays available
      ## levels for a factor
      conditionalPanel("input.filter1 != 'None'",
        checkboxGroupInput("filter1_cb", label='', choices='')
      ),
      
      textInput("custom_filter",
        label="Custom Filter",
        value=""
      )
      
      #       h3("Order Rows by Proportion?"),
      #       checkboxInput("orderByProportion", label="Yes / No")
      
    ),
    
    gridsterItem(row=1, col=2, sizex=1, sizey=2,
      h4( style="text-align: center;",
        "Cytokine Proportion by Sample"
      ),
      plotOutput("heatmap", width=width, height=height*2-25),
      checkboxInput("flip_heatmap", "Flip Axes?", value=FALSE)
    ),
    
    gridsterItem(row=1, col=3, sizex=1, sizey=1,
      tags$div(
        selectInput("individual",
          label="Individual",
          choices=sort(unique(as.character(meta$PTID)))
        ),
        plotOutput("linechart", width=width, height=height-85),
        checkboxInput("flip_linechart", "Flip Axes?", value=FALSE)
      )
    ),
    
    gridsterItem(row=2, col=3, sizex=1, sizey=1,
      #       selectInput("sample",
      #         label="Sample",
      #         choices=unique(as.character(meta$name))
      #       ),
      plotOutput("dofplot", width=width, height=height)
    ),
    
    gridsterItem(row=3, col=1, sizex=3, sizey=2,
      ## custom plot output -- set style manually
      tags$div( style=paste0(
        "width: ", width*3, "px; ",
        "height: ", height*2, "px; "
      ),
        tags$div( id="boxplot_by_cytokine", class="shiny-plot-output",
          style=paste0(
            "width: ", width*3, "px; ",
            "height: ", height*2-30, "px; ",
            "margin: 0 auto;"
          )
        ),
        tags$div( style="overflow: auto;",
          tags$div( style="float: left; display: inline-block;", 
            selectInput("boxplot_by_cytokine_orientation", 
              label="Boxplot Orientation",
              choices=c("Horizontal", "Vertical")
            )
          ),
          tags$div( style="float: right; display: inline-block; margin-top: 20px;",
            checkboxInput("boxplot_coord_flip",
              label="Flip Axes?"
            )
          )
        )
      )
    ),
    
    gridsterItem(row=4, col=1, sizex=3, sizey=2,
      tags$div( style="overflow: auto; width: 1290px; height: 600px;",
        h2("Summary Statistics"),
        tableOutput("stats")
      )
    )
    
    #     gridsterItem(row=5, col=1, sizex=1, sizey=1,
    #       tags$div( style=paste("width:", width, "; height:", height),
    #         showOutput("rchart", "polycharts")
    #       )
    #     )
    
    #     gridsterItem(row=4, col=1, sizex=3, sizey=1,
    #       verbatimTextOutput("debug")
    #     )
    
  )
  
))
