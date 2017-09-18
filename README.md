# Mechanical vibration 2016-2017 #

## file list ##

* dataMass load provided data of the 3-mass system
* error_statespace_full/proportional.m function which return the error between the data provided in dataMass and a system made with the argument of the function (masses and damping). full and proportional correspond to full damping or proportional damping
* pp_(various).m post processing scripts, they create figures and numbers in the report folder to write the report
* smartColorPlot.m bullshit: a script created to maximize the difference and the viewtifulness f the plot colors
* step_response_analysis.m perform the step response analysis 
* workspace_load_complsory.m useful quantities to load like the ratios required
* **main.m** it contains all the things to do and it call autonomously all the scripts
* Matrix_iteration_method.m it performs the matrix iteration method