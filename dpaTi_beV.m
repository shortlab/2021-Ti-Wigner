%% Determining the dpa for Ti in the MITR using sigma_D from JANIS

clear

% Atomic mass of each isotope = COLUMN 1

data={46;47;48;49;50};


% Input: MITR energies (MITR_E in eV) and flux (MITR_Phi in n/cm^2/s) 
% taken from "WATF2 Detailed Spectrum.xslx"

load('MITR.mat')

% Flux in "WATF2 Detailed Spectrum.xslx" is for 6MW not 5.5MW

MITR_Phi=(5.5/6)*MITR_Phi;


% Input: JANIS cross sections for isotopes of Ti (eV and b.eV) = COLUMN 2
% ENDF/B-VIII.0 MT=444 (D) interpolated with (lots of) values per decade.

% Within this cell array column, energy = double column 1 
%                                cross-section = double column 2

for i=1:size(data,1) 
    data(i,2)={readmatrix(strcat('Ti',num2str(data{i,1}),'-beV.csv'),'FileType','text')};
end


% Isotopic abundances for Ti (atomic fraction) = COLUMN 3
% https://ciaaw.org/titanium.htm

data(:,3)={0.0825;0.0744;0.7372;0.0541;0.0518};


% Energy transferred gamma: 4mM/(m+M)^2 = COLUMN 4

for i=1:size(data,1)                        
    data(i,4)={4*data{i,1}/(data{i,1}+1)^2};
end


% Input parameters for DPA equation

time=3600*9619/5.5; % (s) numbers from Dave Carpenter's email: 9619 MWh at 5.5 MW.

E_d=30; % Displacement energy (eV) 
        % from G Was Fundamentals of Radiation Materials Science (Table 2.2, p88)


for i=1:size(data,1)
    
    % Define integral limits minE to maxE
    
    minE=E_d/data{i,4}; % = E_d/gamma (eV)
    maxE=19E6;      % upper limit of JANIS data (eV)
    
    
    % Turn the energy limits into their corresponding row index
    
    [minindex]=find(MITR_E>minE,1);
    [maxindex]=find(MITR_E>maxE,1);
    
    
    % Evaluate the JANIS cross-section at MITR flux energies
    
    JANIS_XS=interp1(data{i,2}(:,1),data{i,2}(:,2),MITR_E(minindex:maxindex));
    
    
    % DPA calculation: = COLUMN 5
    % from G Was Fundamentals of Radiation Materials Science (Eq. 2.124, p117)

    % DPA = kappa*(b-to-cm^2)*time*sum(DamageCross-section*Flux*)/2*E_d;

    data(i,5)={0.8*1E-24*time*sum(JANIS_XS.*MITR_Phi(minindex:maxindex))/(2*E_d)};
    
end 

%% Average the DPA over the different Ti isotopes and display!

DPA=cell2mat(data(:,3))'*cell2mat(data(:,5))


%% Legend for 'data' cell array:

% 1: isotopic mass
% 2: 1 - energy
%    2 - cross-section
% 3: isotopic abundance
% 4: energy transferred (gamma)
% 5: isotopic DPA
