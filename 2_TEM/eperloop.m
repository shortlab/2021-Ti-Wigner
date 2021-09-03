clear

%% Calculating the energy per <a> loop in Ti.

R=25E-9/2;                      % m, radius of dislocation loop.
G=45E9;                         % N/m^2, shear modulus.
b=0.295E-9;                     % m, Burgers vector.
nu=0.37;                        % Poisson's ratio.
alpha=1.25;                     % between 0.5-2, dislocation core parameter, H+L p.232.
expy=1.17;                      % between 1.13-1.21 in metals? See Hirth + Lothe p.161.
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

loopdensity=dscpeak1*massdensity/e_tot_J;           % = 85E20 /m^3


%% Calculate line energy & density from 1983Griffiths' number density of loops.

griffdensity=6E20;                                  % /m^3 ± 0.5E19
griffenergy=griffdensity*e_tot_J/massdensity;       % = 0.0254 J/g      *****
griffdisndensity=griffdensity*circumference;        % = 4.7E13 /m^2


%% Calculate (isolated) vacancy concentration corresponding to unaccounted for energy.

vacenergy=1.43*1.6E-19;                             % J, from 1992Ackland.
edifference=dscpeak1-griffenergy;                   % <insert> J/g unaccounted for energy.
vacdensity=edifference*massdensity/vacenergy;       % = <insert> /m^3
vacpercent=100*vacdensity/numberdensity;            % = <insert> at.%


%% Calculate dislocation line energy & density from simulations.

natoms=492800;                      % #, Number of atoms in simulation cell.
cellmass=natoms*atomicmass/mole;    % g, Mass of simulation cell.
totallinelength=1200E-10;           % m, total dislocation line length.
cellvolume=8956310;                 % A^3, Volume of simulation cell.

totallineenergy=totallinelength*eperlength/cellmass;    % = 7.4 J/g
totallinedensity=totallinelength/(cellvolume*1E-10^3);  % = 1.3E16 /m^2

