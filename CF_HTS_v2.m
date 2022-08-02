clear;
clc;
%%% Pick Scenorio 1
tic;
prompt1='number of runs=';
n=input(prompt1);
prompt2='Throughput (Gbps)=';
TP=input(prompt2);
prompt3='Advitisement Highest Speed (Mbps)=';
ADV=input(prompt3);
R1=0.543;
R2=0.373;
user_act=0.0591;
N_max=floor(TP*1024/(ADV*R1*R2)/user_act);

t1=['The Maximum number of subscribers is '  num2str(N_max) '!'];
disp(t1);

prompt4='Discount rate (in decimal)=';
DR=input(prompt4);
prompt5='Terminal Cost (in US $)=';
R_t=input(prompt5);

prompt6='Scenorio Selection (Type in 1 for Best Scenorio; 2 for Norminal; 3 for worst) =';
Scenorio=input(prompt6);

if Scenorio==3;
prompt7='If you select worst Scenorio, please specific your Decay Rate per year (in decimal)=';
Decay=input(prompt7)/4;
end

prompt8='what is estimated mean of your ARPU? (in US $)';
ARPU_mu=input(prompt8);

prompt9='what is estimated STD of your ARPU? (in US $)';
ARPU_sigma=input(prompt9);

L=[];
for t=1:60
if Scenorio==1
    if t<6
    L_single=0.017*t;
    L=[L L_single];
    elseif t<41 && t>5
    L_single=0.3959*log(t/4)-0.0028;
    L=[L L_single];
    else 
    L_single=0.909;
    L=[L L_single];
    end
elseif Scenorio==2
    if t<6
    L_single=0.017*t;
    L=[L L_single];
    elseif t<20 && t>5
    L_single=0.3959*log(t/4)-0.0028;
    L=[L L_single];
    else 
    L_single=0.614;
    L=[L L_single];
    end
else
    if t<6
    L_single=0.017*t;
    L=[L L_single]; 
    elseif t<20 && t>5
    L_single=0.3959*log(t/4)-0.0028;
    L=[L L_single];
    else 
    L_single=0.614-Decay*(t-19);
    L=[L L_single];
    end
end
end

prompt10='Initial annual operational cost (in million US$)=';
C_ops=input(prompt10);

prompt11='Annual operational growth rate (in decimal)=';
gamma_ops=input(prompt11);

prompt12='subscriber acquisition cost ARPU factor (integer)=';
CaAR=input(prompt12);

prompt13='IPS Initial annual cost (in million US$)=';
a=input(prompt13);

prompt14='IPS annual growth rate (in decimal)=';
gamma_ips=input(prompt14);

prompt15='Percentage factor of ARPU per quarterly (in decimal)=';
b=input(prompt15);

prompt16='Insurance rate (in decimal)=';
IR=input(prompt16);

prompt20='lowerbound of launch (in million US$)=';
c_ln_l=input(prompt20);

prompt21='upperbound of launch (in million US$)=';
c_ln_u=input(prompt21);

C_ops_f=[];
for i1=1:60
    C_ops_f1=(1+gamma_ops)^(i1/4);
    C_ops_f=[C_ops_f C_ops_f1];
end

C_ips_f=[];
for i2=1:60
    C_ips_f1=(1+gamma_ips)^(i2/4);
    C_ips_f=[C_ips_f C_ips_f1];
end

if Scenorio~3
N=N_max*L;
dN=N-[0 N(1:59)];
dDR=[];
for i3=1:60
    dDR1=(1+DR)^(i3/4);
    dDR=[dDR dDR1];
end

R=[];
C_aioc=[];

for simu=1:n   
    
ARPU=normrnd(ARPU_mu,ARPU_sigma); 

Re=(ARPU*3*N+R_t*dN)./dDR/10^6;
R_1=sum(Re);
R=[R R_1];

C_q=(C_ops*C_ops_f/4+(ARPU*CaAR*dN)/10^6+a*C_ips_f/4+b*ARPU*N/10^6)./dDR;
C_aioc_1=sum(C_q);
C_aioc=[C_aioc C_aioc_1];
end

else 
N=N_max*L;
dN=N-[0 N(1:59)];
dDR=[];
for i3=1:60
    dDR1=(1+DR)^(i3/4);
    dDR=[dDR dDR1];
end

R=[];
C_aioc=[];

for simu=1:n    
ARPU=normrnd(ARPU_mu,ARPU_sigma);    
Re=(ARPU*3*N+[R_t*dN(1:19) zeros(1,41)])./dDR/10^6;
R_1=sum(Re);
R=[R R_1];

C_q=(C_ops*C_ops_f/4+[ARPU*CaAR*dN(1:19) zeros(1,41)]/10^6+a*C_ips_f/4+b*ARPU*N/10^6)./dDR;
C_aioc_1=sum(C_q);
C_aioc=[C_aioc C_aioc_1];
end
        
end

C_ioc=[];
for i4=1:n
    c_ln=c_ln_l+(c_ln_u-c_ln_l)*rand;
    C_acq=167.28*TP^0.114;
    C_ioc_1=(1+IR)*(C_acq+c_ln);
    C_ioc=[C_ioc C_ioc_1];
end

NPV1=R-C_aioc-C_ioc;

P_neg_npv=sum(NPV1<0)/n;

t2=['The probability of NPV less than zero is '  num2str(P_neg_npv) '!'];
disp(t2);

%%% Pick Scenorio 2

prompt1='number of runs=';
n=input(prompt1);
prompt2='Throughput (Gbps)=';
TP=input(prompt2);
prompt3='Advitisement Highest Speed (Mbps)=';
ADV=input(prompt3);
R1=0.543;
R2=0.373;
user_act=0.0591;
N_max=floor(TP*1024/(ADV*R1*R2)/user_act);

t1=['The Maximum number of subscribers is '  num2str(N_max) '!'];
disp(t1);

prompt4='Discount rate (in decimal)=';
DR=input(prompt4);
prompt5='Terminal Cost (in US $)=';
R_t=input(prompt5);

prompt6='Scenorio Selection (Type in 1 for Best Scenorio; 2 for Norminal; 3 for worst) =';
Scenorio=input(prompt6);

if Scenorio==3;
prompt7='If you select worst Scenorio, please specific your Decay Rate per year (in decimal)=';
Decay=input(prompt7)/4;
end

prompt8='what is estimated mean of your ARPU? (in US $)';
ARPU_mu=input(prompt8);

prompt9='what is estimated STD of your ARPU? (in US $)';
ARPU_sigma=input(prompt9);

L=[];
for t=1:60
if Scenorio==1
    if t<6
    L_single=0.017*t;
    L=[L L_single];
    elseif t<41 && t>5
    L_single=0.3959*log(t/4)-0.0028;
    L=[L L_single];
    else 
    L_single=0.909;
    L=[L L_single];
    end
elseif Scenorio==2
    if t<6
    L_single=0.017*t;
    L=[L L_single];
    elseif t<20 && t>5
    L_single=0.3959*log(t/4)-0.0028;
    L=[L L_single];
    else 
    L_single=0.614;
    L=[L L_single];
    end
else
    if t<6
    L_single=0.017*t;
    L=[L L_single]; 
    elseif t<20 && t>5
    L_single=0.3959*log(t/4)-0.0028;
    L=[L L_single];
    else 
    L_single=0.614-Decay*(t-19);
    L=[L L_single];
    end
end
end

prompt10='Initial annual operational cost (in million US$)=';
C_ops=input(prompt10);

prompt11='Annual operational growth rate (in decimal)=';
gamma_ops=input(prompt11);

prompt12='subscriber acquisition cost ARPU factor (integer)=';
CaAR=input(prompt12);

prompt13='IPS Initial annual cost (in million US$)=';
a=input(prompt13);

prompt14='IPS annual growth rate (in decimal)=';
gamma_ips=input(prompt14);

prompt15='Percentage factor of ARPU per quarterly (in decimal)=';
b=input(prompt15);

prompt16='Insurance rate (in decimal)=';
IR=input(prompt16);

prompt20='lowerbound of launch (in million US$)=';
c_ln_l=input(prompt20);

prompt21='upperbound of launch (in million US$)=';
c_ln_u=input(prompt21);

C_ops_f=[];
for i1=1:60
    C_ops_f1=(1+gamma_ops)^(i1/4);
    C_ops_f=[C_ops_f C_ops_f1];
end

C_ips_f=[];
for i2=1:60
    C_ips_f1=(1+gamma_ips)^(i2/4);
    C_ips_f=[C_ips_f C_ips_f1];
end

if Scenorio~3
N=N_max*L;
dN=N-[0 N(1:59)];
dDR=[];
for i3=1:60
    dDR1=(1+DR)^(i3/4);
    dDR=[dDR dDR1];
end

R=[];
C_aioc=[];

for simu=1:n   
    
ARPU=normrnd(ARPU_mu,ARPU_sigma); 

Re=(ARPU*3*N+R_t*dN)./dDR/10^6;
R_1=sum(Re);
R=[R R_1];

C_q=(C_ops*C_ops_f/4+(ARPU*CaAR*dN)/10^6+a*C_ips_f/4+b*ARPU*N/10^6)./dDR;
C_aioc_1=sum(C_q);
C_aioc=[C_aioc C_aioc_1];
end

else 
N=N_max*L;
dN=N-[0 N(1:59)];
dDR=[];
for i3=1:60
    dDR1=(1+DR)^(i3/4);
    dDR=[dDR dDR1];
end

R=[];
C_aioc=[];

for simu=1:n    
ARPU=normrnd(ARPU_mu,ARPU_sigma);    
Re=(ARPU*3*N+[R_t*dN(1:19) zeros(1,41)])./dDR/10^6;
R_1=sum(Re);
R=[R R_1];

C_q=(C_ops*C_ops_f/4+[ARPU*CaAR*dN(1:19) zeros(1,41)]/10^6+a*C_ips_f/4+b*ARPU*N/10^6)./dDR;
C_aioc_1=sum(C_q);
C_aioc=[C_aioc C_aioc_1];
end
        
end

C_ioc=[];
for i4=1:n
    c_ln=c_ln_l+(c_ln_u-c_ln_l)*rand;
    C_acq=167.28*TP^0.114;
    C_ioc_1=(1+IR)*(C_acq+c_ln);
    C_ioc=[C_ioc C_ioc_1];
end

NPV2=R-C_aioc-C_ioc;

P_neg_npv=sum(NPV2<0)/n;

t2=['The probability of NPV less than zero is '  num2str(P_neg_npv) '!'];
disp(t2);

%%% Pick Scenorio 3

prompt1='number of runs=';
n=input(prompt1);
prompt2='Throughput (Gbps)=';
TP=input(prompt2);
prompt3='Advitisement Highest Speed (Mbps)=';
ADV=input(prompt3);
R1=0.543;
R2=0.373;
user_act=0.0591;
N_max=floor(TP*1024/(ADV*R1*R2)/user_act);

t1=['The Maximum number of subscribers is '  num2str(N_max) '!'];
disp(t1);

prompt4='Discount rate (in decimal)=';
DR=input(prompt4);
prompt5='Terminal Cost (in US $)=';
R_t=input(prompt5);

prompt6='Scenorio Selection (Type in 1 for Best Scenorio; 2 for Norminal; 3 for worst) =';
Scenorio=input(prompt6);

if Scenorio==3;
prompt7='If you select worst Scenorio, please specific your Decay Rate per year (in decimal)=';
Decay=input(prompt7)/4;
end

prompt8='what is estimated mean of your ARPU? (in US $)';
ARPU_mu=input(prompt8);

prompt9='what is estimated STD of your ARPU? (in US $)';
ARPU_sigma=input(prompt9);

L=[];
for t=1:60
if Scenorio==1
    if t<6
    L_single=0.017*t;
    L=[L L_single];
    elseif t<41 && t>5
    L_single=0.3959*log(t/4)-0.0028;
    L=[L L_single];
    else 
    L_single=0.909;
    L=[L L_single];
    end
elseif Scenorio==2
    if t<6
    L_single=0.017*t;
    L=[L L_single];
    elseif t<20 && t>5
    L_single=0.3959*log(t/4)-0.0028;
    L=[L L_single];
    else 
    L_single=0.614;
    L=[L L_single];
    end
else
    if t<6
    L_single=0.017*t;
    L=[L L_single]; 
    elseif t<20 && t>5
    L_single=0.3959*log(t/4)-0.0028;
    L=[L L_single];
    else 
    L_single=0.614-Decay*(t-19);
    L=[L L_single];
    end
end
end

prompt10='Initial annual operational cost (in million US$)=';
C_ops=input(prompt10);

prompt11='Annual operational growth rate (in decimal)=';
gamma_ops=input(prompt11);

prompt12='subscriber acquisition cost ARPU factor (integer)=';
CaAR=input(prompt12);

prompt13='IPS Initial annual cost (in million US$)=';
a=input(prompt13);

prompt14='IPS annual growth rate (in decimal)=';
gamma_ips=input(prompt14);

prompt15='Percentage factor of ARPU per quarterly (in decimal)=';
b=input(prompt15);

prompt16='Insurance rate (in decimal)=';
IR=input(prompt16);

prompt20='lowerbound of launch (in million US$)=';
c_ln_l=input(prompt20);

prompt21='upperbound of launch (in million US$)=';
c_ln_u=input(prompt21);

C_ops_f=[];
for i1=1:60
    C_ops_f1=(1+gamma_ops)^(i1/4);
    C_ops_f=[C_ops_f C_ops_f1];
end

C_ips_f=[];
for i2=1:60
    C_ips_f1=(1+gamma_ips)^(i2/4);
    C_ips_f=[C_ips_f C_ips_f1];
end

if Scenorio~3
N=N_max*L;
dN=N-[0 N(1:59)];
dDR=[];
for i3=1:60
    dDR1=(1+DR)^(i3/4);
    dDR=[dDR dDR1];
end

R=[];
C_aioc=[];

for simu=1:n   
    
ARPU=normrnd(ARPU_mu,ARPU_sigma); 

Re=(ARPU*3*N+R_t*dN)./dDR/10^6;
R_1=sum(Re);
R=[R R_1];

C_q=(C_ops*C_ops_f/4+(ARPU*CaAR*dN)/10^6+a*C_ips_f/4+b*ARPU*N/10^6)./dDR;
C_aioc_1=sum(C_q);
C_aioc=[C_aioc C_aioc_1];
end

else 
N=N_max*L;
dN=N-[0 N(1:59)];
dDR=[];
for i3=1:60
    dDR1=(1+DR)^(i3/4);
    dDR=[dDR dDR1];
end

R=[];
C_aioc=[];

for simu=1:n    
ARPU=normrnd(ARPU_mu,ARPU_sigma);    
Re=(ARPU*3*N+[R_t*dN(1:19) zeros(1,41)])./dDR/10^6;
R_1=sum(Re);
R=[R R_1];

C_q=(C_ops*C_ops_f/4+[ARPU*CaAR*dN(1:19) zeros(1,41)]/10^6+a*C_ips_f/4+b*ARPU*N/10^6)./dDR;
C_aioc_1=sum(C_q);
C_aioc=[C_aioc C_aioc_1];
end
        
end

C_ioc=[];
for i4=1:n
    c_ln=c_ln_l+(c_ln_u-c_ln_l)*rand;
    C_acq=167.28*TP^0.114;
    C_ioc_1=(1+IR)*(C_acq+c_ln);
    C_ioc=[C_ioc C_ioc_1];
end

NPV3=R-C_aioc-C_ioc;

P_neg_npv=sum(NPV3<0)/n;

t2=['The probability of NPV less than zero is '  num2str(P_neg_npv) '!'];
disp(t2);

figure

[muhat1,sigma1] = normfit(NPV1);
% construct histogram, F number of elements in each bin, X center of each
% bin.
% [F1,X1] = hist(NPV1,300);
% % get ready to plot as probability histogram
% F1 = F1/trapz(X1,F1);
% bar(X1,F1); hold on;
% use muhat and sigma to construct pdf
x1 = muhat1-3*sigma1:0.01:muhat1+3*sigma1;
% plot PDF over histogram
y1 = normpdf(x1,muhat1,sigma1);
plot(x1,y1,'r','linewidth',2);
xlabel('Net Present Value (in million US$)');
ylabel('Probability Density Function')

hold on

[muhat2,sigma2] = normfit(NPV2);
% construct histogram, F number of elements in each bin, X center of each
% bin.
% [F2,X2] = hist(NPV2,300);
% % get ready to plot as probability histogram
% F2 = F2/trapz(X2,F2);
% bar(X2,F2); hold on;
% use muhat and sigma to construct pdf
x2 = muhat2-3*sigma2:0.01:muhat2+3*sigma2;
% plot PDF over histogram
y3 = normpdf(x2,muhat2,sigma2);
plot(x2,y3,'r','linewidth',2);


hold on

[muhat3,sigma3] = normfit(NPV3);
% construct histogram, F number of elements in each bin, X center of each
% bin.
% [F3,X3] = hist(NPV3,300);
% % get ready to plot as probability histogram
% F3 = F3/trapz(X3,F3);
% bar(X3,F3); hold on;
% use muhat and sigma to construct pdf
x3 = muhat3-3*sigma3:0.01:muhat3+3*sigma3;
% plot PDF over histogram
y5 = normpdf(x3,muhat3,sigma3);
plot(x3,y5,'r','linewidth',2);



figure
y2 = normcdf(x1,muhat1,sigma1);
plot(x1,y2,'r','linewidth',2);
xlabel('Net Present Value (in million US$)');
ylabel('Culmulative Distribution Function');

hold on

y4 = normcdf(x2,muhat2,sigma2);
plot(x2,y4,'r','linewidth',2);

hold on

y6 = normcdf(x3,muhat3,sigma3);
plot(x3,y6,'r','linewidth',2);



% Roic=R./(C_aioc+C_ioc);
% 
% figure
% [muhat1,sigma1] = normfit(Roic);
% % construct histogram, F number of elements in each bin, X center of each
% % bin.
% [F1,X1] = hist(Roic,300);
% % get ready to plot as probability histogram
% F1 = F1/trapz(X1,F1);
% bar(X1,F1); hold on;
% % use muhat and sigma to construct pdf
% x1 = muhat1-3*sigma1:0.01:muhat1+3*sigma1;
% % plot PDF over histogram
% y1 = normpdf(x1,muhat1,sigma1);
% plot(x1,y1,'r','linewidth',2);
% xlabel('Ratio between Revenue and Life Cycle Cost');
% ylabel('Probability')


toc;