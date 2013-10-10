library(lattice)
library(rCharts) ## for DataTable
library(data.table.extras)
library(shinyGridster)

zoomButton <- function(inputId) {
  tags$i(class="icon-zoom-in", style="position: absolute; right: 10px; bottom: 10px")
}

splomOutput <- function(outputId) {
  tags$div( id="splom-container",
    tags$div(id=outputId, class="splom")
  )
}

#meta <- readRDS("data/RV144 CaseControl/meta.rds")
dat <- readRDS("data/RV144 CaseControl/res.rds")
dat_post <- readRDS("data/RV144 CaseControl/cd4/cd4_sample_proportions_bg_subtracted.rds")

## hacky fix for the new data for the case/control RV144 data
## ie, the meta-data is pulled from dat_post, rather than the separately defined
## 'meta' file as we have merged new information into it
meta <- collapse(dat_post, by=PTID)[, 1:grep("Stim", names(dat_post)), with=FALSE]

d <- dat_post[ Stim %nin% c("sebctrl", "negctrl 1", "negctrl 2"), ]

## width, height for gridster + plot elements
width <- 430
height <- 300

## svg output
svgOutput <- function(outputId, width, height) {
  tags$div(
    tag("svg", list(id=outputId, width=width, height=height, class="html-shiny-output"))
  )
}

## ensure that each matrix has the same column names
stopifnot( length( table( table( unlist( lapply( dat, names ) ) ) ) ) != 1 )

## cytokines
cytokines <- unname( colnames(dat[[1]]) )
cytokines_positive <- paste0( colnames(dat[[1]]), "+" )
cytokines_negative <- paste0( colnames(dat[[1]]), "-" )

shinyUI( bootstrapPage(
  
  includeScript("www/js/d3.js"),
  includeCSS("www/css/styles.css"),
  includeScript("www/js/fancyboxify.js"),
  
  includeScript("www/jquery-ui/js/jquery-ui-1.10.3.custom.min.js"),
  includeCSS("www/jquery-ui/css/ui-lightness//jquery-ui-1.10.3.custom.min.css"),
  
  includeScript("www/multiselect/multiselect.js"),
  includeCSS("www/multiselect/multiselect.css"),
  
  includeCSS("www/css/shinySplom.css"),
  includeScript("www/js/shinySplom.js"),
  
  includeScript("www/js/gridsterExtras.js"),
  includeScript("scripts.js"),
  
  singleton( tags$body( style="background-color: #789;" ) ),
  
  h1(style="text-align: center; color: white", "Cytokine Visualization"),
  
  ## Controls exist separate of gridster layout
  tags$div( id="gridster-control-container",
    
    tags$button( type="button", id="gridster-control-hide",
      "Show Controls"
    ),
    
    tags$button( type="button", id="gridster-control",
      "Gridster Enabled"
    ),
    
    tags$div(
      id='controls-container', 
      selectInput("phenotype", label="Phenotype", choices=list(
        `Log Fold Change`="LogFoldChange",
        `Proportion (rel. Total)`="PropTotal",
        `Proportion (rel. Total Activated`="PropActivated",
        `Proportion (BG Corrected, rel. Total)`="PropTotalBG",
        `Proportion (BG Corrected, rel. Total Activated)`="PropActivatedBG"
      )),
      
      ## multiselect requires the attribute 'multiple' to be set; can't set
      ## this thru regular shiny HTML functions
      HTML("<select id='cytokines' multiple='multiple'>"),
      HTML(
        paste0("<option value='", cytokines_positive, "'> ",
        cytokines, "</option>")
      ),
      HTML("</select>"),
      
#       selectInput("cytokines",
#         label="Cytokine combinations must contain...",
#         choices=c(
#           cytokines_positive,
#           cytokines_negative
#         ),
#         multiple=TRUE
#       ),
      
      tags$div( class="overflow-auto",
        checkboxGroupInput("cytokines_to_exclude",
          label="Do not include cytokine combinations containing...",
          choices=cytokines
        )
      ),
      
      tags$div( class="overflow-auto",
        tags$div( style="float: left; width: 50%;",
          numericInput("cytokine_filter",
            label='Remove cytokine combinations with average phenotype < x',
            value=0
          )
        ),
        
        tags$div( style="float: right; width: 50%;",
          numericInput("max_combos_to_show",
            label=HTML("Show <span style='font-family: monospace;'>n</span> most highly expressed cytokine combinations"),
            value=5
          )
        )
      ),
      
      ## overflow: auto keeps div from collapsing to zero height
      ## see: http://stackoverflow.com/questions/218760/how-do-you-keep-parents-of-floated-elements-from-collapsing
      tags$div(
        tags$div( style="width: 50%; float: left;",
          tags$label( `for`="cytokine_order_min", "Minimum Cytokine Order Combination"),
          tags$input( id="cytokine_order_min", type="number", value="1", min="1", max=ncol( dat[[1]] ), step="1" )
        ),
        tags$div( style="width: 50%; float: right;", 
          tags$label( `for`="cytokine_order_max", "Maximum Cytokine Order Combination"),
          tags$input( id="cytokine_order_max", type="number", value="2", min="1", max=ncol( dat[[1]] ), step="1" )
        )
      ),
      
      #h3("Facets"),
      tags$div(
        
        tags$div( style="width: 33%; float: left;",
          selectInput("facet1",
            label="Facet 1",
            choices=c("Original Ordering", gsub("name", "Sample", names(meta)))
          )
        ),
        
        tags$div( style="width: 33%; float: left;",
          selectInput("facet2",
            label="Facet 2",
            choices=c("None", gsub("name", "Sample", names(meta)))
          )
        ),
        
        tags$div( style="width: 33%; float: left;",
          selectInput("facet3",
            label="Facet 3",
            choices=c("None", gsub("name", "Sample", names(meta)))
          )
        )
        
      ),
      
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
      
    )
  ),
  
  ## Actual gridster object
  gridster( width=width, height=height,
    
    gridsterItem(row=1, col=2, sizex=2, sizey=1,
      h4( style="text-align: center;",
        "Phenotype by Sample"
      ),
      plotOutput("heatmap", width=width*2, height=height-60),
      checkboxInput("flip_heatmap", "Flip Axes?", value=FALSE),
      zoomButton("zoom-heatmap")
    ),
    
    gridsterItem(row=1, col=1, sizex=1, sizey=1,
      tags$div(
        selectInput("individual",
          label="Individual",
          choices=sort(unique(as.character(meta$PTID)))
        ),
        plotOutput("linechart", width=width, height=height-85),
        checkboxInput("flip_linechart", "Flip Axes?", value=FALSE)
      ),
      zoomButton("zoom-linechart")
    ),
    
    gridsterItem(row=2, col=1, sizex=1, sizey=1,
      #       selectInput("sample",
      #         label="Sample",
      #         choices=unique(as.character(meta$name))
      #       ),
      plotOutput("dofplot", width=width, height=height-20),
      checkboxInput("flip_dofplot", "Flip Axes?", value=FALSE),
      zoomButton("zoom-dofplot")
    ),
    
    ## NOTE: we set the associated CSS for the DataTable generated
    gridsterItem(row=2, col=2, sizex=2, sizey=1,
      chartOutput("stats", "datatables"),
      zoomButton("zoom-stats")
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
            "height: ", height*2-30, "px; "
          )
        ),
        tags$div( style="overflow: auto;",
          
          tags$div( style="float: left; display: inline-block; margin-right: 20px;",
            selectInput("plot_type",
              label="Type of Plot To View",
              choices=list(
                `Box Plot`="boxplot",
                `Histogram`="histogram",
                `Density Plot`="density"
              )
            )
          ),
          
          tags$div( style="float: left; display: inline-block;", 
            selectInput("boxplot_by_cytokine_orientation", 
              label="Boxplot Orientation",
              choices=c("Vertical", "Horizontal")
            )
          ),
          tags$div( style="float: left; display: inline-block; margin-left: 20px;",
            checkboxInput("boxplot_manual_limits", "Set Limits Manually?", FALSE)
          ),
          conditionalPanel("input.boxplot_manual_limits === true",
            tags$div( style="float: left; display: inline-block; margin-left: 20px;",
              numericInput("boxplot_lower_limit", "Lower Limit", 0)
            ),
            tags$div( style="float: left; display: inline-block; margin-left: 20px;",
              numericInput("boxplot_upper_limit", "Upper Limit", 1)
            )
          ),
          tags$div( style="float: right; display: inline-block; margin-top: 20px;",
            checkboxInput("boxplot_coord_flip",
              label="Flip Axes?"
            )
          )
          
        ),
        zoomButton("zoom-boxplot")
      )
    ),
    
    gridsterItem(row=3, col=1, sizex=2, sizey=2,
      splomOutput("splom"),
      zoomButton("zoom-splom")
    )
    
  )
  
))
