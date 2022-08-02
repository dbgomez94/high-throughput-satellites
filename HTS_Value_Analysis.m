clear;clc;

% Plotting histogram
% F: number of elements in each bin
% X: center of each bins

%% Best Case
% load('HTS_Value_Analysis_Results_Best_V2.mat');
% [F_best,X_best] = hist(NPV,300);
% F_best = F_best/trapz(X_best,F_best);
% x_best = muhat-3*sigma:0.01:muhat+3*sigma;
% y_pdf_best = normpdf(x_best,muhat,sigma);
% y_cdf_best = normcdf(x_best,muhat,sigma);
% 
% fig = figure(2);
% Plot_PMF_CDF(X_best, F_best, x_best, y_pdf_best, y_cdf_best,fig)
% hold on

%% Nominal Case
% load('HTS_Value_Analysis_Results_Nominal_V2.mat');
% [F_nom,X_nom] = hist(NPV,300);
% F_nom = F_nom/trapz(X_nom,F_nom);
% x_nom = muhat-3*sigma:0.01:muhat+3*sigma;
% y_pdf_nom = normpdf(x_nom,muhat,sigma);
% y_cdf_nom = normcdf(x_nom,muhat,sigma);
% 
% Plot_PMF_CDF(X_nom, F_nom, x_nom, y_pdf_nom, y_cdf_nom,fig)

%% Worst Case
load('HTS_Value_Analysis_Results_Worst_V2.mat');
[F_worst,X_worst] = hist(NPV,300);
F_worst = F_worst/trapz(X_worst,F_worst);
x_worst = muhat-3*sigma:0.01:muhat+3*sigma;
y_pdf_worst = normpdf(x_worst,muhat,sigma);
y_cdf_worst = normcdf(x_worst,muhat,sigma);

Plot_PMF_CDF(X_worst, F_worst, x_worst, y_pdf_worst, y_cdf_worst,fig)