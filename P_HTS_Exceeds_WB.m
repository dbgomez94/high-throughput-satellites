x = 1:15;
x_99 = x.*99;
P = [98 95 90 83 72 59 45 32 21 12 6 3 1 0.5 0.2]./100;
P_interp = interp1(x,P,1:0.001:15);
figure(1)
plot(1:0.001:15,P_interp)
xlabel('E(NPV_HTS)/E(NPV_WB)');
ylabel('Cumulative Distribution Curve (CDF)');
grid on

