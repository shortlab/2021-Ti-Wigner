clear

load('A1_loop_dia.mat')

%% Calculating the energy per <a> loop in Ti.

R=1E-9*mean(A1)/2;              % m, average radius of dislocation loop.
G=45E9;                         % N/m^2, shear modulus.
b=0.295E-9;                     % m, Burgers vector.
nu=0.37;                        % Poisson's ratio.
alpha=1.25;                     % between 0.5-2, dislocation core parameter, H+L p.232.
expy=1.17;                      % between 1.13-1.21 in metals. See Hirth + Lothe p.161.
gammab=0.145E-3;                % J/m^2, prismatic plane generalized SFE.

circumference=2*pi*R;           % m, circumference of dislocation loop.
area=pi*R^2;                    % m^2, area of dislocation loop.

eperlength=((G*b^2)/(4*pi*(1-nu)))*(log((4*R)/((b*expy)/(2*alpha)))-1);     % J/m, energy per length.
eperlength_eV_nm=1E-9*eperlength/(1.6E-19);                                 % eV/nm, energy per length.

e_line=circumference*eperlength;    % J, loop line energy.
e_area=area*gammab;                 % J, loop stacking fault energy.
e_area_eV=e_area/1.6E-19;           % eV, loop stacking fault energy.         

e_tot_J=e_line+e_area;          % J / loop, total energy per dislocation loop.
e_tot_eV=e_tot_J/1.6E-19;       % eV / loop, total energy per dislocation loop.


%% Calculate the number density of dislocation loops from DSC.

dscpeak1=0.36;                  % J/g, stored energy in ROI 1 from DSC.

massdensity=4.5E6;                                  % g/m^3
atomicmass=47.867;                                  % g/mole
mole=6.02E23;                                       % atoms/mole
numberdensity=massdensity*mole/atomicmass;          % atoms/m^3

loopdensity=dscpeak1*massdensity/e_tot_J;           % = 1.19E22 /m^3


%% Calculate loop energy & density from INL TEM number density of loops.

inldensity=3.5E21;                                  % /m^3 ±
inlenergy=inldensity*e_tot_J/massdensity            % = 0.106 J/g      *****


%% Calculate upper limit of loop energy

% Different values
max_R=1E-9*(mean(A1)+std(A1)/sqrt(size(A1,1)))/2;        % m, average radius of dislocation loop.
max_alpha=2;                            % between 0.5-2, dislocation core parameter, H+L p.232.
max_expy=1.13;                          % between 1.13-1.21 in metals. See Hirth + Lothe p.161.
max_gammab=0.17E-3;                     % J/m^2, prismatic plane generalized SFE.
max_inldensity=inldensity/0.9;          % Thickness/mean free path accuracy ± 10%

max_circumference=2*pi*max_R;           % m, circumference of dislocation loop.
max_area=pi*max_R^2;                    % m^2, area of dislocation loop.

max_eperlength=((G*b^2)/(4*pi*(1-nu)))*(log((4*max_R)/((b*max_expy)/(2*max_alpha)))-1);     % J/m, energy per length.
max_eperlength_eV_nm=1E-9*max_eperlength/(1.6E-19);                                 % eV/nm, energy per length.

max_e_line=max_circumference*max_eperlength;    % J, loop line energy.
max_e_area=max_area*max_gammab;                 % J, loop stacking fault energy.
max_e_area_eV=max_e_area/1.6E-19;               % eV, loop stacking fault energy.         
max_e_tot_J=max_e_line+max_e_area;              % J / loop, total energy per dislocation loop.
max_e_tot_eV=max_e_tot_J/1.6E-19;               % eV / loop, total energy per dislocation loop.

max_inlenergy=max_inldensity*max_e_tot_J/massdensity    % = 0.135 J/g  *****
inlenergy_pluserror=max_inlenergy-inlenergy             % = 0.029 J/g  *****


%% Calculate lower limit of loop energy

% Different values
min_R=1E-9*(mean(A1)-std(A1)/sqrt(size(A1,1)))/2;        % m, average radius of dislocation loop.
min_alpha=0.5;                          % between 0.5-2, dislocation core parameter, H+L p.232.
min_expy=1.21;                          % between 1.13-1.21 in metals. See Hirth + Lothe p.161.
min_gammab=0.12E-3;                     % J/m^2, prismatic plane generalized SFE.
min_inldensity=inldensity/1.1;          % Thickness/mean free path accuracy ± 10%

min_circumference=2*pi*min_R;           % m, circumference of dislocation loop.
min_area=pi*min_R^2;                    % m^2, area of dislocation loop.

min_eperlength=((G*b^2)/(4*pi*(1-nu)))*(log((4*min_R)/((b*min_expy)/(2*min_alpha)))-1);     % J/m, energy per length.
min_eperlength_eV_nm=1E-9*min_eperlength/(1.6E-19);                                 % eV/nm, energy per length.

min_e_line=min_circumference*min_eperlength;    % J, loop line energy.
min_e_area=min_area*min_gammab;                 % J, loop stacking fault energy.
min_e_area_eV=min_e_area/1.6E-19;               % eV, loop stacking fault energy.         
min_e_tot_J=min_e_line+min_e_area;              % J / loop, total energy per dislocation loop.
min_e_tot_eV=min_e_tot_J/1.6E-19;               % eV / loop, total energy per dislocation loop.

min_inlenergy=min_inldensity*min_e_tot_J/massdensity    % = 0.074 J/g  *****
inlenergy_minuserror=inlenergy-min_inlenergy            % = 0.032 J/g  *****
