	.. -*- mode: rst -*-
  
R data reduction and visualisation
==================================

R routines performing a principal components
analysis (PCA) on behavioral data, and regression
analysis with the outputed principal components.
Those routines also useful visualization functions
for correlation matrix plot, regression results plot,
principal components projections plot.

Prerequisites
=============

Data table
~~~~~~~~~~

You must have in the .xlsx format or .csv format, a **structured table**
containing for all your subjects their clinical description. The table 
can contanining categorical variables, or continous variable. It is advised
that you at least one column dedicated to the subject names, it 
will be useful for plotting purpose for example to pinpoint a subject
data point.

Software
~~~~~~~~

You will need to have on your machine the following software:

* R, at least the 3.6 version.

In R, the following libraries are mandatory:
* factoextra, 
* FactoMineR,
* readxl,
* psych,
* corrplot,
* easyGgplot2,
* ggrepel,
* GPArotation
* NbClust


Scripts
=======



The only script to execute is the `pipeline_launcher.sh` scripts. 
In this scripts, you will have to changes the data access path, 
to your image database, lobar atlas, and anatomical atlas. You'll
have to change the path to the `scripts` folder too.
The results will be created in a **Roblob** folder inside each
subject folder with all the sub-folder corresponding to each
step of the pipelines. I encourage you to check each step.

Contributors
============

`Dhaif BEKHA`_

.. _Dhaif BEKHA: dhaif@dhaifbekha.com
