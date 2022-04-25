clear

%% Input data in 10^12 dyn/cm^2 (1964Fisher)

Temp=[4;23;73;123;173;223;273;298;323;373;423;473;523;573;623;673;723;773;823;873;923;973;1023;1073];
c_11=[1.761;1.759;1.749;1.726;1.699;1.668;1.639;1.624;1.609;1.579;1.551;1.522;1.495;1.468;1.442;1.416;1.392;1.368;1.345;1.322;1.299;1.276;1.253;1.231];
c_33=[1.905;1.905;1.894;1.876;1.857;1.837;1.816;1.807;1.795;1.774;1.753;1.734;1.715;1.696;1.678;1.661;1.644;1.627;1.610;1.593;1.576;1.560;1.545;1.529];
c_44=[0.508;0.508;0.505;0.499;0.490;0.481;0.472;0.467;0.462;0.453;0.444;0.434;0.424;0.414;0.404;0.392;0.381;0.370;0.359;0.348;0.337;0.326;0.316;0.307];
c_66=[0.446;0.446;0.439;0.425;0.405;0.384;0.363;0.352;0.342;0.323;0.304;0.285;0.267;0.250;0.234;0.219;0.205;0.191;0.178;0.166;0.154;0.142;0.130;0.118];
c_13=[0.683;0.682;0.680;0.681;0.684;0.687;0.689;0.690;0.691;0.694;0.695;0.695;0.692;0.692;0.691;0.690;0.692;0.688;0.688;0.688;0.688];
c_12=[0.869;0.867;0.871;0.877;0.889;0.901;0.913;0.920;0.925;0.934;0.943;0.952;0.961;0.967;0.973;0.978;0.983;0.985;0.988;0.991;0.992;0.993;0.994;0.996];

%% Convert from 10^12 dyn/cm^2 to GPa

c_11=1E-9*1E-5*1E12*1E4*c_11; % GPa
c_33=1E-9*1E-5*1E12*1E4*c_33; % GPa
c_44=1E-9*1E-5*1E12*1E4*c_44; % GPa
c_66=1E-9*1E-5*1E12*1E4*c_66; % GPa
c_13=1E-9*1E-5*1E12*1E4*c_13; % GPa
c_12=1E-9*1E-5*1E12*1E4*c_12; % GPa

%% Compile C matrices as a f(T)  (2011Tromans)

maxdata=21; % length of c_13 matrix

for T=1:maxdata 
    C_T{T,1}=Temp(T);
    C_T{T,2}=[c_11(T),c_12(T),c_13(T),    0,      0,      0; ...
              c_12(T),c_11(T),c_13(T),    0,      0,      0; ...
              c_13(T),c_13(T),c_33(T),    0,      0,      0; ...
              0,      0,      0,          c_44(T),0,      0; ...
              0,      0,      0,          0,      c_44(T),0; ...
              0,      0,      0,          0,      0,      c_66(T)];
end

%% Calculate prismatic energy factor (1976Savin)

K3_T=C_T;

for T=1:maxdata
    K3_T{T,2}=(C_T{T,2}(1,1)^2-C_T{T,2}(1,2)^2)/(2*C_T{T,2}(1,1));
end

K3=cell2mat(K3_T);

%% Calculate basal energy factor (1976Savin)

lambdasq_T=K3_T;

for T=1:maxdata
    lambdasq_T{T,2}=(C_T{T,2}(1,1)/C_T{T,2}(3,3))^0.5;
end

K1_T=K3_T;

for T=1:maxdata
    K1_T{T,2}=(lambdasq_T{T,2}*C_T{T,2}(3,3)+C_T{T,2}(1,3))*...
                ((C_T{T,2}(4,4)*(lambdasq_T{T,2}*C_T{T,2}(3,3)-C_T{T,2}(1,3)))/...
                (C_T{T,2}(3,3)*(lambdasq_T{T,2}*C_T{T,2}(3,3)+C_T{T,2}(1,3)+2*C_T{T,2}(4,4))))^0.5;
end

K1=cell2mat(K1_T);

%% Evaluate K_eff at a given temperature

chosenTemp=823;                 % K
[chosenTempRow] = find(K3(:,1)==chosenTemp,1);

% disp(strcat('K3=',num2str(K3(chosenTempRow,2))))
% disp(strcat('K1=',num2str(K1(chosenTempRow,2))))

maxFeret=19.2;                  % nm, Average maximum Feret diameter from the as-irradiated sample.
minFeret=11.2;                  % nm, Average minimum Feret diameter from the as-irradiated sample.
FeretRatio=maxFeret/minFeret;   % ratio of K3=c to K1=a

K_eff=  FeretRatio*K3(chosenTempRow,2)/(FeretRatio+1) + ...
        1*K1(chosenTempRow,2)/(FeretRatio+1);
    
    %% Calculating the energy per <a> loop in Ti.

R=19E-9/2;                      % m, average radius of dislocation loop.
b=0.295E-9;                     % m, Burgers vector.
K_eff=39.1598E9;                % N/m^2 - from ti.m
alpha=1.25;                     % between 0.5-2, dislocation core parameter, see Hirth + Lothe p.232.
expy=1.17;                      % between 1.13-1.21 in metals. See Hirth + Lothe p.161.

circumference=2*pi*R;           % m, circumference of dislocation loop.

eperlength_J_m=((K_eff*b^2)/(4*pi))*(log((4*R)/((b*expy)/(2*alpha)))-1);        % J/m, energy per length.
eperlength_eV_nm=1E-9*eperlength_J_m/(1.6E-19);                                 % eV/nm, energy per length.      

eperloop_J=circumference*eperlength_J_m;  % J / loop, total energy per dislocation loop.
eperloop_eV=eperloop_J/1.6E-19;       % eV / loop, total energy per dislocation loop.

massdensity=4.5E6;                                  % g/m^3
atomicmass=47.867;                                  % g/mole
mole=6.02E23;                                       % atoms/mole
numberdensity=massdensity*mole/atomicmass;          % atoms/m^3


%% Calculate energy from TEM number density of loops.

A1_loopdensity=4E21;                                     % /m^3 ± 19% due to uncertainty in sample thickness.
A1_loopenergy=A1_loopdensity*eperloop_J/massdensity      % J/g ± 19% due to uncertainty in sample thickness.

