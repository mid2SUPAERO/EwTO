clear all
close all
 
%0.3 poisson coefficient mean value for metals in top88
%Al alloy Stainless Steel Ti alloy inconel 
Emat=[70.8, 197, 115, 205].*10^9;
rhomat=[2795, 7915, 4425, 7900];
co2mat=[13, 6.15, 40.4, 15.5];
L=2; %m
h=0.5; %m
delta_max=0.005;
F=20000;
lveh=100000000;  %km
 
life=25;
FRC=103; %tonco2/100kg/year
lifekero=FRC*life*1000/100;  %kgco2/kg
%from ADEME : jet A in France or europe :  3.83kgeCO2/kg.
emitkero=3.83; %kgco2 / kg kerosen
lifeco2=lifekero*emitkero;
co2veh=lifeco2/lveh; %kgCO2/km

%raw data from top88 Pareto
load('complHRr3.csv')
cP1=complHRr3(:,2);
%filtering
   win=1000;
   xgauss=0:1:win-1;
   sig=win/8;
   ygauss=1/(sig*sqrt(2*pi))*exp(-0.5*(xgauss-(win-1)/2).^2./sig^2);
   cFiltG=conv(cP1,ygauss);
   cFiltGT=cFiltG(win:end-win);
%
%cpareto=complHRr3(:,2); %raw pareto %Pando % multistart
cpareto=cFiltGT;
vpareto= 0.01:0.0001:1;
vpareto=vpareto(win/2:end-win/2-1);% 
figure(1)
plot(vpareto, cpareto);
 
figure(2)
plot(vpareto,vpareto'.*cpareto);
 
[optimalvfv,optimalvf]=min(vpareto'.*cpareto);
optVf=vpareto(optimalvf);
 
%Al alloy Stainless Steel Ti alloy inconel 
for material=1:4 % 2 or 3 or 4
   
thick(material)=cpareto(optimalvf)*F/delta_max/Emat(material);
mass(material)=L*h*thick(material)*optVf*rhomat(material);

%Idx_veh=(Emat(material)/rhomat(material))*(co2mat(material)+lveh*co2veh)
%Idx_bridge=(Emat(material)/rhomat(material))*(co2mat(material))

Impact_CO2veh(material)=(co2mat(material)+lveh*co2veh)*mass(material);
 
Impact_CO2bridge(material)=(co2mat(material))*mass(material);
end
 
thick
mass
Impact_CO2veh
Impact_CO2bridge
