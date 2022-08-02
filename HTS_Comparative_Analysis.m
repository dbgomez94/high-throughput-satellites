clc;close all;

% Plotting histogram
% F: number of elements in each bin
% X: center of each bins

%% HTS Results
% % Best Case
load('HTS_Value_Analysis_Results_Best_V2.mat');
disp(muhat)
disp(sigma)
fig = figure(1);
Plot_PMF_CDF(X_best, F_best, x_best, y_pdf_best, y_cdf_best,fig)
hold on

% Nominal Case
load('HTS_Value_Analysis_Results_Nominal_V2.mat');
disp(muhat)
disp(sigma)
Plot_PMF_CDF(X_nom, F_nom, x_nom, y_pdf_nom, y_cdf_nom,fig)

% Worst Case
load('HTS_Value_Analysis_Results_Worst_V2.mat');
disp(muhat)
disp(sigma)
Plot_PMF_CDF(X_worst, F_worst, x_worst, y_pdf_worst, y_cdf_worst,fig)

%% HTS and WB Comparative Analysis
% WB Results

load('WB_Value_Analysis_Results_V3.mat')
[muhat,sigma] = normfit(NPV);
[F,X] = hist(NPV,300);
F = F/trapz(X,F);
x = muhat-3*sigma:0.01:muhat+3*sigma;
y_pdf = normpdf(x,muhat,sigma);
y_cdf = normcdf(x,muhat,sigma);

% WB vs Nominal Case
fig_com = figure();
% Wide Beam
Plot_PMF_CDF(X, F, x, y_pdf, y_cdf,fig_com)
hold on
% Nominal Case

load('HTS_Value_Analysis_Results_Nominal_V2.mat');
Plot_PMF_CDF(X_nom, F_nom, x_nom, y_pdf_nom, y_cdf_nom,fig_com)
