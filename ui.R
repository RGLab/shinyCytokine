#library(rCharts)
library(shinyGridster)

meta <- readRDS("data/meta.rds")
dat <- readRDS("data/res_cd8.rds")

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

## custom gridster item
gridItem <- function(..., col=NULL, row=NULL, sizex=NULL, sizey=NULL) {
  return( gridsterItem(...,
    style="border-radius: 10px;",
    row=row,
    col=col,
    sizex=sizex,
    sizey=sizey
  ) )
}

shinyUI( bootstrapPage(
  
  ## note: 'overflow: auto;' is used as a cheap way to get the
  ## checkboxes in a div line up side by side
  tags$head( tags$style( type="text/css","

    .gridster {
      width: 1400px;
      margin: 0 auto;
    }

    .gridster > * {
      margin: 0 auto;
    }

    #cytokines {
      overflow: auto;
    }

    #cytokines .checkbox {
      float: left;
      margin: 5px;
    }

    #facet1 {
      width: 100%;
    }

    #facet2 {
      width: 100%;
    }

    #phenotype {
      width: 100%;
    }

    #individual {
      width: 100%;
    }
    ")),
  
  singleton( tags$body( style="background-color: #789;" ) ),
  
  h1(style="text-align: center; color: white", "Cytokine Visualization"),
  
  gridster( width=width, height=height,
    
    gridItem( row=1, col=1, sizex=1, sizey=2,
      
      tags$div( style="overflow: auto;",
        
        tags$div( style="width: 45%; float: left;",
          selectInput("phenotype", label="Phenotype", choices=list(
            `Proportion`="Proportion",
            `Proportion (BG Corrected)`="Proportion_bg"
          ))
        ),
        
        tags$div( style="width: 45%; float: right;", 
          selectInput("individual",
            label="Individual",
            choices=sort(unique(as.character(meta$PTID)))
          )
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
    
    gridItem(row=1, col=2, sizex=1, sizey=2,
      plotOutput("heatmap", width=width, height=height*2, clickId="heatmap_click", hoverId="heatmap_hover")
    ),
    
    gridItem(row=1, col=3, sizex=1, sizey=1,
      plotOutput("linechart", width=width, height=height)
    ),
    
    gridItem(row=2, col=3, sizex=1, sizey=1,
      plotOutput("dofplot", width=width, height=height)
    ),
    
    gridItem(row=3, col=1, sizex=3, sizey=2,
      plotOutput("boxplot_by_cytokine", width=width*3, height=height*2)
    ),
    
    gridItem(row=4, col=1, sizex=3, sizey=2,
      tags$div( style="overflow: auto; width: 1290px; height: 600px;",
        h2("Summary Statistics"),
        tableOutput("stats")
      )
    )
    
    #     gridItem(row=5, col=1, sizex=1, sizey=1,
    #       tags$div( style=paste("width:", width, "; height:", height),
    #         showOutput("rchart", "polycharts")
    #       )
    #     )
    
    #     gridItem(row=4, col=1, sizex=3, sizey=1,
    #       verbatimTextOutput("debug")
    #     )
    
  )
  
))
