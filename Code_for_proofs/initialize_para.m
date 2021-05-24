function para=initialize_para(pH,Va)

% Initialization of the (many) parameters used in the model. This is based 
% on code from Calipso (adimensionnement.m and physicochimie.m).

adimensionnement %Uses pH and Va

para.V=Vad;
para.OmegaFe=ohmega_fer;
para.OmegaOx=ohmega_ox;
para.pH=pH;
para.kd0=kd_0_ad;
para.ad0=ad_0;
para.kappa=para.OmegaFe/(4*para.OmegaOx); 
para.PB=PB_ratio;

%% Parameters for the Poisson equation

para.alpha0=alpha_0;
para.alpha1=alpha_1;
para.lambda=sqrt(Lambda_2);
para.rho=Rho;
para.DPsi0=delta_Psi0_pzc; 
para.DPsi1=delta_Psi1_pzc;
para.zC=2;
para.zN=-1;
para.zP=3;

%% Parameters for the convection-diffusion equations

para.epsC=D1/D3;
para.epsN=D1/D2;
para.epsP=1;

%% Parameters for the bc of the oxigen vacancies

para.mC0=mC_0;
para.kC0=kC_0;

para.mC1=mC_1;
para.kC1=kC_1;

para.aC0=aC_0;
para.bC0=bC_0;

para.aC1=aC_1;
para.bC1=bC_1;


%% Parameters for the bc of the electrons

para.mN0=mN_0;
para.pN0=pN_0;
para.nN0=nN_0;
para.kN0=kN_0;

para.mN1=mN_1;
para.kN1=kN_1;

para.aN0=aN_0;
para.bN0=bN_0;

para.ar0=ar_0;
para.br0=br_0;

para.Nmetal=N_metal;

%% Parameters for the bc of the Fe

para.mP0=mP_0;
para.kP0=kP_0;

para.mP1=mP_1;
para.kP1=kP_1;

para.aP0=aP_0;
para.bP0=bP_0;

para.aP1=aP_1;
para.bP1=bP_1;

para.Pm=Pm;

%% Some scaling parameters, mostly used to draw pictures in physical units               
para.resc_l=L0/10^-9;%from nondimensional to nm
para.resc_delta=D1/L0*3600*24*365.2425/10^-6;%from nondimensional to \mu m/year
para.resc_tildeJ=D1*Faraday/(ohmega_ox*L0)*100;%from nondimensional to \mu A/cm^2
para.resc_V=R*T/Faraday;%from nondimensional to V


%% For initializing the solution
para.coupling=0;% For the numerical continuation. The source term is the 
                % poisson equation for the electric potential is given by
                % para.coupling*(2C-N+3P+rho). Of course we are only really
                % interested in the case where para.coupling=1, but it is
                % easier to find numerical solutions when para.coupling=0
                % (using the paper from Claire and Thomas), and starting
                % from there we can numerically continue the solution up
                % to para.coupling=1.
                
%% For later
para.case='Poten'; 

%% BC order 0
%%% C
para.betaC0=@(x) para.mC0*exp(-2*para.bC0*x)/4+para.kC0*exp(2*para.aC0*x)/4;
para.betaC1=@(x) para.mC1*exp(-3*para.bC1*x)/4+para.kC1*exp(3*para.aC1*x)/4;
para.gammaC0=@(x) para.mC0*exp(-2*para.bC0*x);
para.gammaC1=@(x) para.kC1*exp(3*para.aC1*x);
para.rC0=@(X) para.betaC0(X(2))*X(1)-para.gammaC0(X(2));
para.rC1=@(X) para.betaC1(X(3)-X(2))*X(1)-para.gammaC1(X(3)-X(2));

%%% N
para.betaN0=@(x) para.kN0*exp(-para.aN0*x)+para.pN0*exp(-para.br0*x);
para.betaN1=@(x) para.mN1;
para.gammaN0=@(x) para.mN0*exp(para.bN0*x)+para.nN0*exp(para.ar0*x);
para.gammaN1=@(x) para.kN1*para.Nmetal*log(1+exp(-x));
para.rN0=@(X) para.betaN0(X(2))*X(1)-para.gammaN0(X(2));
para.rN1=@(X) para.betaN1(X(3)-X(2))*X(1)-para.gammaN1(X(3)-X(2));

%%% P
para.betaP0=@(x) para.mP0*exp(-3*para.bP0*x) + para.kP0*exp(3*para.aP0*x);
para.betaP1=@(x) para.mP1*exp(-3*para.bP1*x) + para.kP1*exp(3*para.aP1*x);
para.gammaP0=@(x) para.mP0*para.Pm*exp(-3*para.bP0*x);
para.gammaP1=@(x) para.kP1*para.Pm*exp(3*para.aP1*x);
para.rP0=@(X) para.betaP0(X(2))*X(1)-para.gammaP0(X(2));
para.rP1=@(X) para.betaP1(X(3)-X(2))*X(1)-para.gammaP1(X(3)-X(2));

%% BC order 1
%%% C
para.derbetaC0=@(x) -2*para.bC0*para.mC0*exp(-2*para.bC0*x)/4 + 2*para.aC0*para.kC0*exp(2*para.aC0*x)/4;
para.derbetaC1=@(x) -3*para.bC1*para.mC1*exp(-3*para.bC1*x)/4 + 3*para.aC1*para.kC1*exp(3*para.aC1*x)/4;
para.dergammaC0=@(x) -2*para.bC0*para.mC0*exp(-2*para.bC0*x);
para.dergammaC1=@(x) 3*para.aC1*para.kC1*exp(3*para.aC1*x);
para.der1rC0=@(X) para.betaC0(X(2));
para.der2rC0=@(X) para.derbetaC0(X(2))*X(1)-para.dergammaC0(X(2));
para.der1rC1=@(X) para.betaC1(X(3)-X(2));
para.der2rC1=@(X) -para.derbetaC1(X(3)-X(2))*X(1)+para.dergammaC1(X(3)-X(2));
para.der3rC1=@(X) para.derbetaC1(X(3)-X(2))*X(1)-para.dergammaC1(X(3)-X(2));

%%% N
para.derbetaN0=@(x) -para.aN0*para.kN0*exp(-para.aN0*x)-para.br0*para.pN0*exp(-para.br0*x);
para.derbetaN1=@(x) 0;
para.dergammaN0=@(x) para.bN0*para.mN0*exp(para.bN0*x)+para.ar0*para.nN0*exp(para.ar0*x);
para.dergammaN1=@(x) -para.kN1*para.Nmetal*exp(-x)/(1+exp(-x));
para.der1rN0=@(X) para.betaN0(X(2));
para.der2rN0=@(X) para.derbetaN0(X(2))*X(1)-para.dergammaN0(X(2));
para.der1rN1=@(X) para.betaN1(X(3)-X(2));
para.der2rN1=@(X) -para.derbetaN1(X(3)-X(2))*X(1)+para.dergammaN1(X(3)-X(2));
para.der3rN1=@(X) para.derbetaN1(X(3)-X(2))*X(1)-para.dergammaN1(X(3)-X(2));

%%% P
para.derbetaP0=@(x) -3*para.bP0*para.mP0*exp(-3*para.bP0*x) + 3*para.aP0*para.kP0*exp(3*para.aP0*x);
para.derbetaP1=@(x) -3*para.bP1*para.mP1*exp(-3*para.bP1*x) + 3*para.aP1*para.kP1*exp(3*para.aP1*x);
para.dergammaP0=@(x) -3*para.bP0*para.mP0*para.Pm*exp(-3*para.bP0*x);
para.dergammaP1=@(x) 3*para.aP1*para.kP1*para.Pm*exp(3*para.aP1*x);
para.der1rP0=@(X) para.betaP0(X(2));
para.der2rP0=@(X) para.derbetaP0(X(2))*X(1)-para.dergammaP0(X(2));
para.der1rP1=@(X) para.betaP1(X(3)-X(2));
para.der2rP1=@(X) -para.derbetaP1(X(3)-X(2))*X(1)+para.dergammaP1(X(3)-X(2));
para.der3rP1=@(X) para.derbetaP1(X(3)-X(2))*X(1)-para.dergammaP1(X(3)-X(2));

%% BC order 2
%%% C
para.der2betaC0=@(x) (-2*para.bC0)^2*para.mC0*exp(-2*para.bC0*x)/4 + (2*para.aC0)^2*para.kC0*exp(2*para.aC0*x)/4;
para.der2betaC1=@(x) (-3*para.bC1)^2*para.mC1*exp(-3*para.bC1*x)/4 + (3*para.aC1)^2*para.kC1*exp(3*para.aC1*x)/4;
para.der2gammaC0=@(x) (-2*para.bC0)^2*para.mC0*exp(-2*para.bC0*x);
para.der2gammaC1=@(x) (3*para.aC1)^2*para.kC1*exp(3*para.aC1*x);
para.der11rC0=@(X) 0;
para.der12rC0=@(X) para.derbetaC0(X(2));
para.der21rC0=@(X) para.derbetaC0(X(2));
para.der22rC0=@(X) para.der2betaC0(X(2))*X(1)-para.der2gammaC0(X(2));
para.der11rC1=@(X) 0;
para.der12rC1=@(X) -para.derbetaC1(X(3)-X(2));
para.der13rC1=@(X) para.derbetaC1(X(3)-X(2));
para.der21rC1=@(X) -para.derbetaC1(X(3)-X(2));
para.der22rC1=@(X) para.der2betaC1(X(3)-X(2))*X(1)-para.der2gammaC1(X(3)-X(2));
para.der23rC1=@(X) -para.der2betaC1(X(3)-X(2))*X(1)+para.der2gammaC1(X(3)-X(2));
para.der31rC1=@(X) para.derbetaC1(X(3)-X(2));
para.der32rC1=@(X) -para.der2betaC1(X(3)-X(2))*X(1)+para.der2gammaC1(X(3)-X(2));
para.der33rC1=@(X) para.der2betaC1(X(3)-X(2))*X(1)-para.der2gammaC1(X(3)-X(2));

%%% N
para.der2betaN0=@(x) (-para.aN0)^2*para.kN0*exp(para.aN0*x)+(-para.br0)^2*para.pN0*exp(-para.br0*x);
para.der2betaN1=@(x) 0;
para.der2gammaN0=@(x) (para.bN0)^2*para.mN0*exp(para.bN0*x)+(para.ar0)^2*para.nN0*exp(para.ar0*x);
para.der2gammaN1=@(x) para.kN1*para.Nmetal*exp(-x)/(1+exp(-x))^2;
para.der11rN0=@(X) 0;
para.der12rN0=@(X) para.derbetaN0(X(2));
para.der21rN0=@(X) para.derbetaN0(X(2));
para.der22rN0=@(X) para.der2betaN0(X(2))*X(1)-para.der2gammaN0(X(2));
para.der11rN1=@(X) 0;
para.der12rN1=@(X) -para.derbetaN1(X(3)-X(2));
para.der13rN1=@(X) para.derbetaN1(X(3)-X(2));
para.der21rN1=@(X) -para.derbetaN1(X(3)-X(2));
para.der22rN1=@(X) para.der2betaN1(X(3)-X(2))*X(1)-para.der2gammaN1(X(3)-X(2));
para.der23rN1=@(X) -para.der2betaN1(X(3)-X(2))*X(1)+para.der2gammaN1(X(3)-X(2));
para.der31rN1=@(X) para.derbetaN1(X(3)-X(2));
para.der32rN1=@(X) -para.der2betaN1(X(3)-X(2))*X(1)+para.der2gammaN1(X(3)-X(2));
para.der33rN1=@(X) para.der2betaN1(X(3)-X(2))*X(1)-para.der2gammaN1(X(3)-X(2));

%%% P
para.der2betaP0=@(x) (-3*para.bP0)^2*para.mP0*exp(-3*para.bP0*x) + (3*para.aP0)^2*para.kP0*exp(3*para.aP0*x);
para.der2betaP1=@(x) (-3*para.bP1)^2*para.mP1*exp(-3*para.bP1*x) + (3*para.aP1)^2*para.kP1*exp(3*para.aP1*x);
para.der2gammaP0=@(x) (-3*para.bP0)^2*para.mP0*para.Pm*exp(-3*para.bP0*x);
para.der2gammaP1=@(x) (3*para.aP1)^2*para.kP1*para.Pm*exp(3*para.aP1*x);
para.der11rP0=@(X) 0;
para.der12rP0=@(X) para.derbetaP0(X(2));
para.der21rP0=@(X) para.derbetaP0(X(2));
para.der22rP0=@(X) para.der2betaP0(X(2))*X(1)-para.der2gammaP0(X(2));
para.der11rP1=@(X) 0;
para.der12rP1=@(X) -para.derbetaP1(X(3)-X(2));
para.der13rP1=@(X) para.derbetaP1(X(3)-X(2));
para.der21rP1=@(X) -para.derbetaP1(X(3)-X(2));
para.der22rP1=@(X) para.der2betaP1(X(3)-X(2))*X(1)-para.der2gammaP1(X(3)-X(2));
para.der23rP1=@(X) -para.der2betaP1(X(3)-X(2))*X(1)+para.der2gammaP1(X(3)-X(2));
para.der31rP1=@(X) para.derbetaP1(X(3)-X(2));
para.der32rP1=@(X) -para.der2betaP1(X(3)-X(2))*X(1)+para.der2gammaP1(X(3)-X(2));
para.der33rP1=@(X) para.der2betaP1(X(3)-X(2))*X(1)-para.der2gammaP1(X(3)-X(2));


%% chi (corresponds to the nonlinear term in F^{\delta})
para.chi=@(x) para.OmegaFe/para.OmegaOx*para.kd0*exp(para.rho*para.ad0*x);
para.derchi=@(x) para.rho*para.ad0*para.OmegaFe/para.OmegaOx*para.kd0*exp(para.rho*para.ad0*x);
para.der2chi=@(x) (para.rho*para.ad0)^2*para.OmegaFe/para.OmegaOx*para.kd0*exp(para.rho*para.ad0*x);


