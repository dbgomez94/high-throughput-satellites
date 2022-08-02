%% Just HTS
load('HTS_Value_Analysis_Results_Nominal_V2.mat');
[mu_NPV,sig_NPV] = normfit(NPV);
[mu_RIOC,sig_RIOC] = normfit(Roic);
Roic=R./(C_aioc+C_ioc);
Roic_max=max(Roic);
Roic_min=min(Roic);
seg_x=(Roic_max-Roic_min)/100;

NPV_max=max(NPV);
NPV_min=min(NPV);
seg_y=(NPV_max-NPV_min)/100;

x_axis=Roic_min:seg_x:Roic_max;
y_axis=NPV_min:seg_y:NPV_max;


figure(1)
histogram2(Roic,NPV,x_axis,y_axis,'Normalization','pdf','FaceColor','flat');
h = colorbar;
colormap bone
set(get(h,'title'),'string','PMF');
view(2)
xlabel('Discounted Return on Invested Capital (ROIC)','fontname','times new roman')
ylabel('Net Present Value ($ million)','fontname','times new roman')
zlabel('Probability Mass Function (PMF)','fontname','times new roman')

%% Comparison of HTS and WB
figure(2)
yyaxis left
histogram2(Roic,NPV,x_axis,y_axis,'Normalization','pdf','FaceColor','flat');
h = colorbar;
colormap jet
set(get(h,'title'),'string','PMF');
view(2)
xlabel('Discounted Return on Invested Capital','fontname','times new roman')
ylabel('Net Present Value ($ million)','fontname','times new roman')
zlabel('Probability Mass Function (PMF)','fontname','times new roman')
hold on

clear;
load('WB_Value_Analysis_Results.mat');
yyaxis right
histogram2(Roic,NPV,x_axis,y_axis,'Normalization','pdf','FaceColor','flat');
h = colorbar;
colormap jet
set(get(h,'title'),'string','PMF');
view(2)
xlabel('Discounted Return on Invested Capital','fontname','times new roman')
ylabel('Net Present Value ($ million)','fontname','times new roman')
zlabel('Probability Mass Function (PMF)','fontname','times new roman')