function results=lnest(iostruct,metadata,options)
%results=LNEST(iostruct,metadata,options)
%
%Estimate SISO (single-input single-output) L-N (linear-nonlinear) model from a set of input-output traces
%Also performs cross-validation.
%
%
%INPUT
%
%iostruct: a struct containing the input output data
%metadata: a struct containing metadata for preprocessing
%options: options for doing the LN estimation
%
%required iostruct fields:
% .inputtraces  matrix of time series COLUMN vectors
% .outputtraces matrix of time series COLUMN vectors
%
%  [inputtraces and outputtraces must be the same size.]
%
%optional iostruct fields:
% .times  A matrix of timestamp vectors, same size as .inputtraces
%         used for jitter correction
%
% .dt (default: .1) the [intended] time interval of the traces, in seconds
%         
%OUTPUT
%results: a struct containing the L-N model estimation results
%
%
%LNEST needs export_fig.m (third party) to save out figures
%
% Saul Kato  saul - at - kato.com
% git:saulkato/wb
%


