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
load('complHRr3Pstrain.csv')
cP1strain=complHRr3Pstrain(:,2);
%filtering
   win=1000;
   xgauss=0:1:win-1;
   sig=win/8;
   ygauss=1/(sig*sqrt(2*pi))*exp(-0.5*(xgauss-(win-1)/2).^2./sig^2);
   cFiltG=conv(cP1,ygauss);
   cFiltGT=cFiltG(win:end-win);
   cFiltGstrain=conv(cP1strain,ygauss);
   cFiltGTstrain=cFiltGstrain(win:end-win);
%
%cpareto=complHRr3(:,2); %raw pareto %Pando % multistart
cpareto=cFiltGT;
cparetostrain=cFiltGTstrain;
vpareto= 0.01:0.0001:1;
vpareto=vpareto(win/2:end-win/2-1);% 
figure(1)
plot(vpareto, cpareto);
hold on
plot(vpareto, cparetostrain);
hold off
 
figure(2)
plot(vpareto,vpareto'.*cpareto);
hold on
plot(vpareto, vpareto'.*cparetostrain);
hold off
legend('plane stress','plane strain')
xlabel('volume fraction (V_f)')
ylabel('V_f*f(V_f)')
 
[optimalvfv,optimalvf]=min(vpareto'.*cpareto);
optVf=vpareto(optimalvf);

[optimalvfvstrain,optimalvfstrain]=min(vpareto'.*cparetostrain);
optVfstrain=vpareto(optimalvfstrain);

%Al alloy Stainless Steel Ti alloy inconel 
%1         2               3         4

thickinit=cpareto(optimalvf)*F/delta_max./Emat;

thick=thickinit;
change=1;
while change>0.001
    r=thick/L  % t over L ratio
    compliance=cpareto(optimalvf)*(r.^(-0.5))/(r.^(0.5)+r.^(-0.5))+cparetostrain(optimalvf)*(r.^(0.5))/(r.^(0.5)+r.^(-0.5));
    thick_new=compliance*F/delta_max./Emat;
    change=abs(thick_new-thick)./thick;
    thick=thick_new;
end

thick_error=(thick-thickinit)./thickinit;

mass=L*h*thick*optVf.*rhomat;

%Idx_veh=(Emat(material)/rhomat(material))*(co2mat(material)+lveh*co2veh)
%Idx_bridge=(Emat(material)/rhomat(material))*(co2mat(material))
 
Impact_CO2veh=(co2mat+lveh*co2veh).*mass;
Impact_CO2bridge=(co2mat).*mass;

thick
mass
Impact_CO2veh
Impact_CO2bridge
