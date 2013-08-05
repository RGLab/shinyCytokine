Comparing SPICE and shinyCytokine
========================================================
author: Kevin Ushey
date: July 26, 2013
font-import: http://fonts.googleapis.com/css?family=Risque
font-family: 'Risque'

Description of shinyCytokine
========================================================

Our Shiny visualization application allows us to
visualize cytokine profiles.

By cytokine profile, we mean 'the proportion of cells
expressing a (particular combination of) cytokines.'

Multiple kinds of plots are generated and these can be
colored and facetted by different meta-data available.

We can also focus on particular 'order's of cytokine
combinations; e.g. focus on 2nd-order combinations.

We can view joint, marginal distributions as well.

Different Plots Available
========================================================

We have a number of different plots available.

- Heatmap of cytokine cell proportions
- Individual-level cytokine profile
- Degree of Functionality
- Facetted boxplots of Proportions

Heatmap
========================================================

Presents proportion of cells expressing a particular
(combination of) cytokines by sample. Rows are samples, 
columns are cytokines.

Allows us to identify highly-expressed combinations of
cytokines. Facetting allows us to see visually whether
these proportions differ by different meta-data.

Heatmap
========================================================
![shinyCytokine Heatmap Example](SPICE_Comparison-img/heatmap.png)

Heatmap
========================================================

Facetting reveals more 'obvious' patterns.
![shinyCytokine Heatmap Example 2](SPICE_Comparison-img/heatmap-facetted.png)

Cytokine Profile
=======================================================
Allows us to view the variation in cytokine proportion,
by sample, for a particular individual.

These plots can be facetted as well.

Cytokine Profile
========================================================
![shinyCytokine Cytokine Profile](SPICE_Comparison-img/cytokine-profile.png)

Degree of Functionality
=======================================================
TODO

Boxplots
=======================================================
Allows us to more readily compare distributions of
proportions by different meta-data variables.

Gives an alternative way of visualizing the distribution
of proportions; similar to the heatmap.

Can facet by 2 variables, colour by 1 other variable.

Boxplots
========================================================
![shinyCytokine Boxplots](SPICE_Comparison-img/boxplots.png)

Description of SPICE
========================================================

[SPICE](http://exon.niaid.nih.gov/SPICE/) is software
for visualizing and performing simple analyses of 
cytokine proportions.

Only runs on Mac OS X. No source code released.

Import data as a single tab/comma delimited file. 

* Each row is
a unique combination of variables; typically a unique 
cytokine profile, + other cell/sample-specific info.

SPICE
========================================================

![SPICE](SPICE_Comparison-img/spice.png)

Main Features
========================================================

* Interactively and easily select different 'roles' for
variables
  * Categorical (use variable in determining unique 'groups' of observations)
  * Overlay (facet pie charts by this variable)
  * Sum, average (collapse proportions over this variable, using sum or average)
  * Group (essentially, 'ignore' or 'collapse' over this variable)
  
Variable Roles
========================================================

![SPICE Variable Roles](SPICE_Comparison-img/spice-variable-roles.png)

***

Similar to the facetting controls in shinyCytokine, but
with more options for collapsing variables.
  
Main Features
========================================================

### Pie charts

Each slice is a [mean/median] of proportions.

Allows one to visualize how 'the [mean/median] for each 
unique combination of categorical variables' differs.

Arcs show marginals.

***

![SPICE Piecharts](SPICE_Comparison-img/spice-piecharts.png)

Criticisms
======================================================

- Cannot get a legend for each pie 'slice'
  - although this information is contained in the other
  bar chart
- No notion of statistical variance in chart
- Potentially misleading visual comparisons of means
  
Main Features
=======================================================

- Bar Chart: Plot bars that extend to the (mean, median, range, IQR)
  of the data for each unique combination
  
- Can augment with points, whiskers, error bars -- however, scale
is not easily modifiable.
  
***

![SPICE Barcharts](SPICE_Comparison-img/spice-barcharts.png)


Main Features
========================================================

- Can modify barchart settings to get something
like a box and whisker plot.

- Current settings: bar at the middle is the mean,
whiskers are SD. Median gets squashed to the bottom.

- It's ugly.

***

![SPICE Dotplots](SPICE_Comparison-img/spice-dotplot.png)

Main Features
=======================================================

### NPlot

Similar to barchart, but plotted as linechart.

- Automatically splits over 'overlay' variables.

- Each line is generated by the 'group'ed variables.

- I average over PTID, Timepoint, Group to get a
sensible image.

***

![SPICE NPlot](SPICE_Comparison-img/spice-nplot.png)

Main Features
=======================================================

### 'CoolPlot' **

Heatmap of 'overlay' variables vs 'categorical' variables.

Similar to the shinyCytokine heatmap, as though we 'overlay' all
'non-cytokine' variables. I think?

Selecting all non-cytokine variables as overlay kills SPICE.

<p style="font-size: 16px;">** Note: the name is 'CoolPlot', not CoolPlot</p>

***

![SPICE Heatmap](SPICE_Comparison-img/spice-heatmap.png)

Pros and Cons
========================================================

### shinyCytokine
\+ Open source

\+ Runs in the browser

\+/- Controls + Plots part of one view

\- Can't save 'state' / open multiple 'views'

\+/- Facetting options limited (less flexible, less complex)


***

### SPICE
\- Closed source

\- Only runs on Mac OS

\+/- Controls and plots somewhat separated

\+ Multiple views / save state

\+/- Can specify different 'roles' for each variable (more flexible, more complex)

Pros and Cons
========================================================

### shinyCytokine
\+ Easily zoom in on figure with double-click

\+ View updates automatically upon parameter change

\- (Currently) no option for plot export

\+ Default plot settings are usually informative

***

### SPICE
\+ Control a global workspace zoom level

\+ Can set automatic updates at time intervals

\+ Can export graphs

\+/- Statistical testing built in

\- Plot defaults are either uninformative or misleading

Pros and Cons
========================================================

### shinyCytokine

\+/- Cannot edit data from within app

\- Cannot save and load state (yet?)

***

### SPICE

\+/- Can add, remove, edit variables from within SPICE

\+ Can save and load workspaces

Similarities
========================================================

shinyCytokine's heatmap and SPICE's 'CoolPlot' are similar

* shinyCytokine is restricted to x-variable being cytokine
  combinations, and y-variable being samples
  
* 'CoolPlot' plots proportions as `overlay ~ categorical`

* shinyCytokine's heatmap can be facetted by a variable; 
  coolPlot cannot
  
* 'CoolPlot' can be resized, but the labels don't resize
  intelligently (become too big / too small)
  
Similarities
========================================================

shinyCytokine's 'cytokine profile' and SPICE's barchart
/ nPlot are similar

* shinyCytokine plots proportion by cytokine for a particular
  individual as points, connected by lines.
  
* SPICE produces barplots up to mean / median by default.
  Misleading. But can plot points instead.
  
* SPICE has more control over what defines a unique 'group'.

Conclusion
========================================================

shinyCytokine is a better tool for visualization of
cytokine data.

SPICE is a better workspace for the visualization,
analysis and modification of cytokine data.

shinyCytokine intends to integrate with other programmatic
workflows. SPICE acts as a stand-alone analysis application.

shinyCytokine can be easily deployed with a particular set of
data for view / usage by other researchers; all they need is a 
web browser.

Conclusion
=======================================================

Ultimately, SPICE is limited relative to using a standard
programming environment.

