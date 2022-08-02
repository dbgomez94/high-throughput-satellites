clear
clc
tic

prompt1='number of runs (how many times of simulation)=';
n=input(prompt1);

prompt2='number of transponders=';
N=input(prompt2);

prompt3='Design lifetime (in yrs) =';
td=input(prompt3);

%first MC option
prompt4='mean of marginal cost of durability (in decimal)=';
Cm_mu=input(prompt4);

prompt5='STD of marginal cost of durability (in decimal)=';
Cm_sigma=input(prompt5);

%second MC option, variation in launch
prompt6='lowerbound of launch (in million US$)=';
c_ln_l=input(prompt6);

prompt7='upperbound of launch (in million US$)=';
c_ln_u=input(prompt7);

%third MC option, catch the insurance variation
prompt8='mean of Insurance rate (in decimal)=';
IR_mu=input(prompt8);

prompt9='STD of Insurance rate (in decimal)=';
IR_sigma=input(prompt9);

Cioc=[];

for simu=1:n
    
if td>=15;
    Cm=normrnd(Cm_mu, Cm_sigma);
    c_ln=c_ln_l+(c_ln_u-c_ln_l)*rand;
    IR=normrnd(IR_mu, IR_sigma);
    Cioc_1=((63.1*log(N)-166.3)/((1+Cm)^(15-td))+c_ln)*(1+IR);
    Cioc=[Cioc Cioc_1];
    
else
    
    Cm=normrnd(Cm_mu, Cm_sigma);
    c_ln=c_ln_l+(c_ln_u-c_ln_l)*rand;
    IR=normrnd(IR_mu, IR_sigma);
    Cioc_1=((63.1*log(N)-166.3)*((1+Cm)^(15-td))+c_ln)*(1+IR);
    Cioc=[Cioc Cioc_1];
    
end
end

t1=['The mean of cost of ioc is '  num2str(mean(Cioc)) ' million US dollars!'];
disp(t1);

prompt10='time delay between acquisition and operation cost begin (in yrs)=';
dT=input(prompt10);

prompt11='initial operation cost (in million US$)=';
Ciops=input(prompt11);

%fourth MC option, can simulate the various different economical
%perspectives
prompt12='average of annual cost growth rate of the operation cost (in decimal)=';
rips_mu=input(prompt12);

prompt13='STD of annual cost growth rate of the operation cost (in decimal)=';
rips_sigma=input(prompt13);

prompt14='discount rate (in decimal)=';
DR=input(prompt14);

prompt15='estimated service year after launch (in yrs)=';
to=input(prompt15);

Cops=[];
for simu=1:n
    Copss=[];
for yr=1:to

rips=normrnd(rips_mu, rips_sigma);
Cops_1=Ciops*((1+rips)^(yr-1))/((1+DR)^(yr+dT));
Copss=[Copss Cops_1];
Cops_1=sum(Copss);

end

Cops=[Cops Cops_1];

end

LCC=Cioc+Cops;

t2=['The mean of cost of operation and Life Cycle Cost (LCC) are '  num2str(mean(Cops)) ' and ' num2str(mean(LCC)) 'million US dollars!'];
disp(t2);

% figure
% [muhat,sigma] = normfit(LCC);
% % construct histogram, F number of elements in each bin, X center of each
% % bin.
% [F,X] = hist(LCC,300);
% % get ready to plot as probability histogram
% F = F/trapz(X,F);
% bar(X,F); hold on;
% % use muhat and sigma to construct pdf
% x = muhat-3*sigma:0.01:muhat+3*sigma;
% % plot PDF over histogram
% y = normpdf(x,muhat,sigma);
% plot(x,y,'r','linewidth',2);
% 
% figure
% y2 = normcdf(x,muhat,sigma);
% plot(x,y2,'r','linewidth',2);

prompt16='steady-state load factor (0 to 1)=';
L0=input(prompt16);

%fifth MC option, give the flexibilities of control the obsolescence
prompt17='estimated mean year of obsolescence occurs (in yrs)=';
Tobs_mu=input(prompt17);

prompt18='estimated STD year of obsolescence occurs (in yrs)=';
Tobs_sigma=input(prompt18);

prompt19='estimated mean of intensity of obsolescence (in yrs)=';
theta_mu=input(prompt19);

prompt20='estimated STD of intensity of obsolescence (in yrs)=';
theta_sigma=input(prompt20);

prompt21='estimated loading factor=';
Lf=input(prompt21);

ts=to+dT;
L=[];
for simu=1:n
    Ls=[];
    for yr=1:ts
        Tobs=normrnd(Tobs_mu, Tobs_sigma);
        theta=normrnd(theta_mu, theta_sigma);
        if yr<=dT;
            L_1=0;
        elseif yr>dT&yr<Tobs;
            L_1=L0*(1-exp(-(yr-dT))/Lf);
        else
            L_1=L0*(1-exp(-(Tobs-dT))/Lf)*exp(-((yr-Tobs)/theta)^2);
        end
        Ls=[Ls L_1];
    end
    L=[L; Ls];
end

Ly=mean(L);
figure
tsn=0:17;
plot(tsn, [0 Ly]);

%sixth MC option monitor the different combination of service, and
%the price flunctuation for each of the revenue streams.
prompt22='mean of fraction of video (in decimal)=';
dv_mu=input(prompt22);

prompt23='STD of fraction video (in decimal)=';
dv_sigma=input(prompt23);

prompt24='mean of fraction of audio (in decimal)=';
da_mu=input(prompt24);

prompt25='STD of fraction audio (in decimal)=';
da_sigma=input(prompt25);

prompt26='mean of price of video per transponder per year (in million US$)=';
Pv_mu=input(prompt26);

prompt27='STD of price of video per transponder per year (in million US$)=';
Pv_sigma=input(prompt27);

prompt28='mean of price of audio per transponder per year (in million US$)=';
Pa_mu=input(prompt28);

prompt29='STD of price of audio per transponder per year (in million US$)=';
Pa_sigma=input(prompt29);

prompt30='mean of price of digital per transponder per year (in million US$)=';
Pd_mu=input(prompt30);

prompt31='STD of price of digital per transponder per year (in million US$)=';
Pd_sigma=input(prompt31);

PVR=[];

for simu=1:n
    dv=normrnd(dv_mu, dv_sigma);
    Pv=normrnd(Pv_mu, Pv_sigma);
    da=normrnd(da_mu, da_sigma);
    Pa=normrnd(Pa_mu, Pa_sigma);
    dd=1-dv-da;
    Pd=normrnd(Pd_mu, Pd_sigma);
    P=dv*Pv+da*Pa+dd*Pd;
    Pick=ceil(rand*n);
    PVR_1=[];
    for yr=1:ts
    PVRs=N*P*L(Pick,yr)/((1+DR)^yr);
    PVR_1=[PVR_1 PVRs];
    PVR_2=sum(PVR_1);
    end
    PVR=[PVR PVR_2];
end
% 
% figure
% subplot(2,1,1)
% 
% [muhat,sigma] = normfit(PVR);
% % construct histogram, F number of elements in each bin, X center of each
% % bin.
% [F,X] = hist(PVR,300);
% % get ready to plot as probability histogram
% F = F/trapz(X,F);
% bar(X,F); hold on;
% % use muhat and sigma to construct pdf
% x = muhat-3*sigma:0.01:muhat+3*sigma;
% % plot PDF over histogram
% y = normpdf(x,muhat,sigma);
% plot(x,y,'r','linewidth',2);
% 
% subplot(2,1,2)
% y2 = normcdf(x,muhat,sigma);
% plot(x,y2,'r','linewidth',2);

NPV=PVR-LCC;

t3=['The mean of Net Present Value is '  num2str(mean(NPV)) ' million US dollars!'];
disp(t3);

figure
subplot(2,1,1)

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

subplot(2,1,2)
y2 = normcdf(x,muhat,sigma);
plot(x,y2,'r','linewidth',2):


toc