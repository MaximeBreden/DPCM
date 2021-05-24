%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametres physico-chimiques de DPCM                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%% Constantes physiques universelles

m_elec = 9.11e-31 ; 	%masse de l'electron en kg
k_bol = 1.38e-23 ;  	%constante de Boltzmann en J/K
Ki_0 = 8.854e-12 ;  	%constante dielectrique du vide en F/m^2
R = 8.314 ;      	    %constante des gaz parfaits en J/K/mol
Faraday = 9.6485e4;     %constante de Faraday en C/mol 	


%% B1. Caract. physico-chimiques de la solution  

% Temperature                                T [K]                  
T = 0.298000D+03;
% Potentiel redox de la solution             Eredox [V/ENH]         
Eredox = -0.300000D+00;
% Activite du cation ferreux                 aFe2 [mol/m^3]         
aFe2 = 0.000000D+00;
% Activite du cation ferrique                aFe3 [mol/m^3]        
aFe3 = 0.000000D+00;


%% B3. Caracteristiques du metal                                                   
% Volume molaire du fer                      ohmega_fer [m^3/mol]   
ohmega_fer = 0.710500D-05;
% Densite d'etat electronique du fer         n_DOS [mol/(m^3.J)]    
n_DOS = 0.135000D+25;


%% B4. Caracteristique de la couche d oxyde  

% Longueur carateristique                    L0 [m]                 
L0 = 0.100000D-08;
% Volume molaire de l'oxyde                  ohmega_ox [m^3/mol]    
ohmega_ox = 0.447400D-04;
% Constante dielectrique rel. de l'oxyde     
Ki = 0.100000D+02;

%% B5. Parametres adimensionnes des equations  

% Charge du reseau hote                      Rho  
Rho= -0.500000D+01 ; 
% Concentr. max. en cations fer              Pm              
Pm = 0.200500D+01 ; 


%% B6. Autres parametres des equations 

% Coeff. de diffusion (cations fer)          D1 [m^2/s]    
D1 = 0.100000D-22;
% Coeff. de diffusion (electrons)            D2 [m^2/s]   
D2 = 0.100000D-05;
% Coeff. de diffusion (lacunes d'oxygene)    D3 [m^2/s]        
D3 = 0.100000D-19;
% Capacite diff. de l'interface ext.         Gama_0 [F/m^2]        
Gama_0= 0.500000D+00;
% Capacite diff. de l'interface int.         Gama_1 [F/m^2]      
Gama_1 = 0.100000D+01; 
% Pot. de charge nulle en x0                 delta_Psi0_pzc_d [V]     
delta_Psi0_pzc_d = 0.190302D+00 -log(10)*R*T/Faraday*pH; 
% Pot. de charge nulle en x1                 delta_Psi1_pzc_d [V]   
delta_Psi1_pzc_d = -0.105302D+00;


%% B7. Fonctions cinetiques (cations fer)  

% Relachement ferrique                       k0_Fe [m/s]      
k0_Fe = 1.000000D-04;
% Coefficient de Butler-Volmer               aP_0         
aP_0 = 0.500000D+00;
% Insertion ferrique                         m0_Fe [m/s]  
m0_Fe = 1.32927D+00;
% Coefficient de Butler-Volmer               bP_0      
bP_0 = 0.500000D+00;
% Direction aval de l'oxydation du metal     m1_Fe [m/s] 
m1_Fe = 1.000000D-09;
% Coefficient de Butler-Volmer               aP_1       
aP_1= 0.500000D+00;
% Direction amont de l'oxydation du metal    k1_Fe [m/s]   
k1_Fe=1.000000D-03;
% Coefficient de Butler-Volmer               bP_1         
bP_1 = 0.500000D+00;


%% B8. Fonctions cinetique (electrons) 

% Reduction de l'eau                         k0_e [m/s]      
k0_e=1.3e-5;
% Coefficient de Butler-Volmer               aN_0            
aN_0 = 0.500000D+00;
% Dependance de k0_e au pH                   npH_e   
npH_e = 0.100000D+01;
% Oxydation de l'hydrogene                   m0_e [mol/(m^2.s)]  
m0_e = 0.000000D+00;
% Coefficient de Butler-Volmer               bN_0             
bN_0  = 0.500000D+00;
% Reduction ferrique                         kr_0 [m/s]  
kr_0 = 2.54255D+06;
% Coefficient de Butler-Volmer               ar_0                
ar_0 = 0.500000D+00;
% Oxydation ferreuse                         mr_0 [mol/(m^2.s)]           
mr_0 = 1.000000D-04;
% Coefficient de Butler-Volmer               br_0                 
br_0 = 0.500000D+00;

%% B9. Fonction cinetiques (lacunes d oxygene)

% D-amont de l'echg. de lacunes d'oxygene    k0_ox [mol/(m^2.s)]    
k0_ox = 1.00000D-02;
% Coefficient de Butler-Volmer               a0_ox                  
aC_0 = 0.500000D+00;
% D-aval de l'echg. de lacunes d'oxygene     m0_ox [mol/(m^2.s)]    
m0_ox = 2.71838D+01;
% Coefficient de Butler-Volmer               b0_ox                   
bC_0 = 0.500000D+00;
% Dependance au pH                           npH_ox                 
npH_ox = 0.200000D+01;
% Croissance de la couche d'oxyde            k1_ox [mol/(m^2.s)]    
k1_ox = 6.0000D-06;
% Coefficient de Butler-Volmer               a1_ox                   
aC_1 = 0.500000D+00;
% de-Croissance de la couche d'oxyde         m1_ox [mol/(m^2.s)]   
m1_ox = 2.73855D-08;
% Coefficient de Butler-Volmer               b1_ox                   
bC_1 = 0.500000D+00;
% Dissolution de la couche d'oxyde           kd_0 [mol/(m^2.s)]     
kd_0 = 2.85400D-09;
% Coefficient de Butler-Volmer               ad_0                   
ad_0 = 0.257000D+00;
% Dependance au pH de kd_0                   npH_d                
npH_d   = 0.500000D+00;

