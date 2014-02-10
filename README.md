ln
==

Linear-Nonlinear Systems Modeling in MATLAB

This MATLAB library provides a suite of linear-nonlinear (L-N) modeling functions. 
It is currently limited to continuous, SISO (single-input single-output) systems.

It allows you to:
*Estimate a linear kernel (filter) and a nonlinear function from a set of input-output trials.
*Use various preprocessing functions.
*Measure performance of the L-N model on estimate and probe/holdout datasets.
*Plot the model and simulated fits.
*Fit a parameterized function to the filter to make a lower-parameter model.

It does not currently provide functionality for:
*Fitting spike data
*Multiple-input systems. 

This code is fairly general, but it was originally developed for analyzing calcium imaging of 
nonspiking neurons in this paper:  http://dx.doi.org/10.1016/j.neuron.2013.11.020



Saul Kato
saulkato.com
