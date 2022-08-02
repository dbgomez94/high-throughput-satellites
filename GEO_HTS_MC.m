clear;
clc;
tic;
% prompt1='number of runs=';
% n=input(prompt1);
n = 100000;

% prompt2='Throughput (Gbps)=';
% TP=input(prompt2);
TP = 100;

% prompt3='Advitisement Highest Speed (Mbps)=';
% ADV=input(prompt3);
ADV = 12;

R1=0.543;
R2=0.373;
user_act=0.0591;
N_max=floor(TP*1024/(ADV*R1*R2)/user_act);

t1=['The Maximum number of subscribers is '  num2str(N_max) '!'];
disp(t1);

% prompt4='Discount rate (in decimal)=';
% DR=input(prompt4);
DR = 0.1;

% prompt5='Terminal Cost (in US $)=';
% R_t=input(prompt5);
R_t = 300;

%************************************************************************
% prompt6='Scenorio Selection (Type in 1 for Best Scenorio; 2 for Norminal; 3 for worst) =';
% Scenorio=input(prompt6);
Scenorio = 1;
%************************************************************************

if Scenorio==3;
% prompt7='If you select worst Scenorio, please specific your Decay Rate per year (in decimal)=';
% Decay=input(prompt7)/4;
Decay = 0.05/4;
end

% prompt8='what is estimated mean of your ARPU? (in US $)';
% ARPU_mu=input(prompt8);
ARPU_mu = 40;

% prompt9='what is estimated STD of your ARPU? (in US $)';
% ARPU_sigma=input(prompt9);
ARPU_sigma = 2.67;

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

% prompt10='Initial annual operational cost (in million US$)=';
% C_ops=input(prompt10);
C_ops = 10;

% prompt11='Annual operational growth rate (in decimal)=';
% gamma_ops=input(prompt11);
gamma_ops = 0;

% prompt12='subscriber acquisition cost ARPU factor (integer)=';
% CaAR=input(prompt12);
CaAR = 8;

% prompt13='IPS Initial annual cost (in million US$)=';
% a=input(prompt13);
a = 10;

% prompt14='IPS annual growth rate (in decimal)=';
% gamma_ips=input(prompt14);
gamma_ips = 0;

% prompt15='Percentage factor of ARPU per quarterly (in decimal)=';
% b=input(prompt15);
b = 0.05;

% prompt16='Insurance rate (in decimal)=';
% IR=input(prompt16);
IR = 0.12;

% prompt20='lowerbound of launch (in million US$)=';
% c_ln_l=input(prompt20);
c_ln_l = 60;

% prompt21='upperbound of launch (in million US$)=';
% c_ln_u=input(prompt21);
c_ln_u = 100;

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

NPV=R-C_aioc-C_ioc;

P_neg_npv=sum(NPV<0)/n;

t2=['The probability of NPV less than zero is '  num2str(P_neg_npv) '!'];
disp(t2);

figure

[muhat,sigma] = normfit(NPV);
% construct histogram, F number of elements in each bin, X center of each
% bin.
[F,X] = hist(NPV,300);
% get ready to plot as probability histogram
F = F/trapz(X,F);
bar(X,F); hold on;
% use muhat and sigma to construct pdf
x = muhat-3*sigma:0.01:muhat+3*sigma;
% plot PDF over histogram
y = normpdf(x,muhat,sigma);
plot(x,y,'r','linewidth',2);
xlabel('Net Present Value (in million US$)');
ylabel('Probability Density Function')

figure
y2 = normcdf(x,muhat,sigma);
plot(x,y2,'r','linewidth',2);
xlabel('Net Present Value (in million US$)');
ylabel('Culmulative Distribution Function');

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




