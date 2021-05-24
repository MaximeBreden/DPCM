%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adimensionnement des parametres physico-chimiques de DPCM                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

physicochimie

gamma = Faraday/(R*T);

%% chimie

% adimensionnement des parametres globaux
Eredox_ad=Eredox*gamma;
Vad=Va*gamma;

% dimensionless constants for ferric release
mP_0= m0_Fe*L0/D1*aFe3;
kP_0= k0_Fe*L0/D1;
	
% dimensionless constants for proton reduction and ferrous release
mN_0= m0_e*L0*ohmega_ox/D2*10^(-npH_e*pH)*exp(-Eredox_ad);
kN_0= k0_e*L0/D2*10^(-npH_e*pH);
nN_0= mr_0*L0*ohmega_ox/D2*aFe2;
pN_0= kr_0*L0/D2*aFe3;
	
% dimensionless outer potential of zero charge
delta_Psi0_pzc = delta_Psi0_pzc_d*gamma;
	
% dimensionless constants for oxygen exchange
mC_0= m0_ox*ohmega_ox*L0/D3*10^(-npH_ox*pH);
kC_0= k0_ox*ohmega_ox*L0/D3;
	
% dimensionless constant for oxide host lattice dissolution
kd_0_ad= ohmega_ox*L0*kd_0*10^(-npH_d*pH)/D1;

% dimensionless constants for Poisson equation
alpha_0= Ki*Ki_0/(Gama_0*L0);

% rapport de Pilling-Bedworth
PB_ratio = ohmega_ox/ohmega_fer;

%% physique

% dimensionless constants for iron oxidation
mP_1= m1_Fe*L0/D1;
kP_1= k1_Fe*L0/D1;
	
% dimensionless constants for electronic exchange
k1_e= sqrt(k_bol*T/2./pi/(1.*m_elec));
m1_e= sqrt(k_bol*T/2./pi/(1.*m_elec));
mN_1= m1_e*L0/D2;
kN_1= k1_e*L0/D2;
N_metal= k_bol*T*n_DOS*ohmega_ox;
	
% dimensionless constants for oxide host lattice growth
mC_1= 4.*L0*m1_ox*ohmega_ox/D3;
kC_1= 4.*L0*k1_ox*ohmega_ox/D3;
	
% dimensionless inner potential of zero charge
delta_Psi1_pzc = delta_Psi1_pzc_d*gamma;
	
% dimmensionless constants for Poisson equation
Lambda_2= Ki*Ki_0*R*T*ohmega_ox/(Faraday^2*L0^2);
alpha_1= Ki*Ki_0/(Gama_1*L0);






