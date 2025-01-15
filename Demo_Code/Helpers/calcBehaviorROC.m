function out = calcBehaviorROC(data_1, data_2, varargin)
% This function takes the behavioral differentiation data between
% anticipatory and non-reward control zones and creates a
% receiver-operating charateristics model. Outputs the model speficiation,
% scores, and AUC results
% NJ created 04/01/23


%reshape to Nx1 structure and combine AZ and CZ data into one
data_all = [data_1; data_2]; %put AZ and CZ data together ready for ROC

%first half-ish is from familiar AZ, rest is from familiar CZ
resp = (1:length(data_all))' <= length(data_1);

%create model for data 
out.mdl = fitglm(data_all, resp,'Distribution','binomial','Link','logit');
out.scores = out.mdl.Fitted.Probability;
[out.X, out.Y, out.T, out.AUC] = perfcurve(resp,out.scores,'true','UseNearest','off');