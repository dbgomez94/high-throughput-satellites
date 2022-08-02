%% Plotting Joint Distribution
close all
load('WB_Value_Analysis_Results_V3.mat');

figure(1)
histogram2(Roic,NPV,x_axis,y_axis,'Normalization','pdf','FaceColor','flat');
h = colorbar;
colormap bone
set(get(h,'title'),'string','PMF');
view(2)
xlabel('Discounted Return on Invested Capital (DROIC)','fontname','times new roman')
ylabel('Net Present Value ($ million)','fontname','times new roman')
zlabel('Probability Mass Function (PMF)','fontname','times new roman')

figure(2)
histogram2(Roic,NPV,x_axis,y_axis,'Normalization','pdf','facecolor','flat');
h = colorbar;
colormap bone
set(get(h,'title'),'string','PMF');
view(3)
xlabel('Discounted Return on Invested Capital (DROIC)','fontname','times new roman')
ylabel('Net Present Value ($ million)','fontname','times new roman')
zlabel('Probability Mass Function (PMF)','fontname','times new roman')


%% Plotting PMF & CDF vs NPV

% Plotting histogram
% F: number of elements in each bin
% X: center of each bins

[muhat,sigma] = normfit(NPV)
[F,X] = hist(NPV,300);
F = F/trapz(X,F);
x = muhat-3*sigma:0.01:muhat+3*sigma;
y_pdf = normpdf(x,muhat,sigma);
y_cdf = normcdf(x,muhat,sigma);

f_wb = figure();
Plot_PMF_CDF(X, F, x, y_pdf, y_cdf,f_wb)
