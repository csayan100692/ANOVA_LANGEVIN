# PBTVMF
# PBTVMF is the R package for testing the homogeneity of mean directions of several Langevin populations.
# install this package in Rstudio by using the command install.packages("PBTVMF")
# Import this package by using library(PBTVMF)
# This package contains 3 functions "PBTVMF3", "PBTVMF4" and "PBTVMF5"
# PBTVMF3 computes the test statistic value and critical value of a given data matrix of three sample Langevin populations 
# PBTVMF4 computes the test statistic value and critical value of a given data matrix of four sample Langevin populations 
# PBTVMF5 computes the test statistic value and critical value of a given data matrix of five sample Langevin populations 
# we have to give two inputs for running the above three functions
# first input is the data matrix which contains a set of matrices where each matrix is a random sample of size n from a Langevin distribution of dimension 3
# second input is the significance level
# these functions give two outputs, test statistic and critical value
# if test statistic value is less than the critical value then we fail to reject the hypothesis of equal mean directions among the Langevin populations.
