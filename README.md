	.. -*- mode: rst -*-
  
R data reduction and visualisation
==================================

R routines performing a principal components
analysis (PCA) on behavioral data, and regression
analysis with the outputed principal components.
To help you with the components interpretation,
some rotation method are also an option, alongside
clutering tecnics.
Those routines contains also useful visualization functions
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

* R, at least the 3.6 version

R libraies
~~~~~~~~~~

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

* PCA_multiples_technics: this script perform a Principal Components Analysis,
  using the classical SVD decomposition, and, depending if the interpretation
  of the components are easy, you can rotate them using the psych packages.
  Rotating the principal components can sometimes help, but depending of 
  the rotation methods, resulting components can no longer be orthogonal.
* cluster_on_pca: This simple R function, takes a `PCA object` from the 
  FactoMineR package, and with a k-means clustering, segment into k 
  categories the projected data points. It can help to interpret the
  principal components meaning, along with your clinical data.
* plot_linear_model_pcs: Fit a linear model between the outputed
  Principal Components, and behavioral data. It's a very common
  technics, to perform a linear regression on the newly created variable
  as long as you can explain them ;)


Contributors
============

`Dhaif BEKHA`_

.. _Dhaif BEKHA: dhaif@dhaifbekha.com
