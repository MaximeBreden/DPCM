function ipara=para2intval(para)

%Turns all the nondimensionalized parameters into intvals. (The
%nondimensionalization itself could be done with interval arithmetic, but
%we did no do so here).

ipara.alpha0=intval(para.alpha0);
ipara.alpha1=intval(para.alpha1);
ipara.lambda=intval(para.lambda);
ipara.rho=intval(para.rho);
ipara.DPsi0=intval(para.DPsi0);
ipara.DPsi1=intval(para.DPsi1);
ipara.V=intval(para.V);
ipara.zC=intval(para.zC);
ipara.epsC=intval(para.epsC);
ipara.zN=intval(para.zN);
ipara.epsN=intval(para.epsN);
ipara.zP=intval(para.zP);
ipara.epsP=intval(para.epsP);
ipara.OmegaFe=intval(para.OmegaFe);
ipara.OmegaOx=intval(para.OmegaOx);
ipara.pH=intval(para.pH);
ipara.kd0=intval(para.kd0);
ipara.ad0=intval(para.ad0);
ipara.kappa=intval(para.kappa);
ipara.mC0=intval(para.mC0);
ipara.mC1=intval(para.mC1);
ipara.aC0=intval(para.aC0);
ipara.aC1=intval(para.aC1);
ipara.bC0=intval(para.bC0);
ipara.bC1=intval(para.bC1);
ipara.kC0=intval(para.kC0);
ipara.kC1=intval(para.kC1);
ipara.ar0=intval(para.ar0);
ipara.br0=intval(para.br0);
ipara.mN0=intval(para.mN0);
ipara.mN1=intval(para.mN1);
ipara.aN0=intval(para.aN0);
ipara.bN0=intval(para.bN0);
ipara.kN0=intval(para.kN0);
ipara.kN1=intval(para.kN1);
ipara.mP0=intval(para.mP0);
ipara.mP1=intval(para.mP1);
ipara.aP0=intval(para.aP0);
ipara.aP1=intval(para.aP1);
ipara.bP0=intval(para.bP0);
ipara.bP1=intval(para.bP1);
ipara.kP0=intval(para.kP0);
ipara.kP1=intval(para.kP1);
ipara.nN0=intval(para.nN0);
ipara.pN0=intval(para.pN0);
ipara.Pm=intval(para.Pm);
ipara.Nmetal=intval(para.Nmetal);
ipara.coupling=intval(para.coupling);%Artifical parameter in front of the coupling terms in the Poisson equation
ipara.case=para.case;

%% BC order 0
%%% C
ipara.betaC0=@(x) 1/4*(ipara.mC0*exp(-2*ipara.bC0*x) + ipara.kC0*exp(2*ipara.aC0*x));
ipara.betaC1=@(x) 1/4*(ipara.mC1*exp(-3*ipara.bC1*x) + ipara.kC1*exp(3*ipara.aC1*x));
ipara.gammaC0=@(x) ipara.mC0*exp(-2*ipara.bC0*x);
ipara.gammaC1=@(x) ipara.kC1*exp(3*ipara.aC1*x);
ipara.rC0=@(X) ipara.betaC0(X(2))*X(1)-ipara.gammaC0(X(2));
ipara.rC1=@(X) ipara.betaC1(ipara.V-X(2))*X(1)-ipara.gammaC1(ipara.V-X(2));

%%% N
ipara.betaN0=@(x) ipara.kN0*exp(-ipara.aN0*x) + ipara.pN0*exp(-ipara.br0*x);
ipara.betaN1=@(x) ipara.mN1;
ipara.gammaN0=@(x) ipara.mN0*exp(ipara.bN0*x) + ipara.nN0*exp(ipara.ar0*x);
ipara.gammaN1=@(x) ipara.kN1*ipara.Nmetal*log(1+exp(-x));
ipara.rN0=@(X) ipara.betaN0(X(2))*X(1)-ipara.gammaN0(X(2));
ipara.rN1=@(X) ipara.betaN1(ipara.V-X(2))*X(1)-ipara.gammaN1(ipara.V-X(2));

%%% P
ipara.betaP0=@(x) ipara.mP0*exp(-3*ipara.bP0*x) + ipara.kP0*exp(3*ipara.aP0*x);
ipara.betaP1=@(x) ipara.mP1*exp(-3*ipara.bP1*x) + ipara.kP1*exp(3*ipara.aP1*x);
ipara.gammaP0=@(x) ipara.mP0*ipara.Pm*exp(-3*ipara.bP0*x);
ipara.gammaP1=@(x) ipara.kP1*ipara.Pm*exp(3*ipara.aP1*x);
ipara.rP0=@(X) ipara.betaP0(X(2))*X(1)-ipara.gammaP0(X(2));
ipara.rP1=@(X) ipara.betaP1(ipara.V-X(2))*X(1)-ipara.gammaP1(ipara.V-X(2));

%% BC order 1
%%% C
ipara.derbetaC0=@(x) 1/4*(-2*ipara.bC0*ipara.mC0*exp(-2*ipara.bC0*x) + 2*ipara.aC0*ipara.kC0*exp(2*ipara.aC0*x));
ipara.derbetaC1=@(x) 1/4*(-3*ipara.bC1*ipara.mC1*exp(-3*ipara.bC1*x) + 3*ipara.aC1*ipara.kC1*exp(3*ipara.aC1*x));
ipara.dergammaC0=@(x) -2*ipara.bC0*ipara.mC0*exp(-2*ipara.bC0*x);
ipara.dergammaC1=@(x) 3*ipara.aC1*ipara.kC1*exp(3*ipara.aC1*x);
ipara.der1rC0=@(X) ipara.betaC0(X(2));
ipara.der2rC0=@(X) ipara.derbetaC0(X(2))*X(1)-ipara.dergammaC0(X(2));
ipara.der1rC1=@(X) ipara.betaC1(X(3)-X(2));
ipara.der2rC1=@(X) -ipara.derbetaC1(X(3)-X(2))*X(1)+ipara.dergammaC1(X(3)-X(2));
ipara.der3rC1=@(X) ipara.derbetaC1(X(3)-X(2))*X(1)-ipara.dergammaC1(X(3)-X(2));

%%% N
ipara.derbetaN0=@(x) -ipara.aN0*ipara.kN0*exp(-ipara.aN0*x) - ipara.br0*ipara.pN0*exp(-ipara.br0*x);
ipara.derbetaN1=@(x) 0;
ipara.dergammaN0=@(x) ipara.bN0*ipara.mN0*exp(ipara.bN0*x) + ipara.ar0*ipara.nN0*exp(ipara.ar0*x);
ipara.dergammaN1=@(x) -ipara.kN1*ipara.Nmetal*exp(-x)/(1+exp(-x));
ipara.der1rN0=@(X) ipara.betaN0(X(2));
ipara.der2rN0=@(X) ipara.derbetaN0(X(2))*X(1)-ipara.dergammaN0(X(2));
ipara.der1rN1=@(X) ipara.betaN1(X(3)-X(2));
ipara.der2rN1=@(X) -ipara.derbetaN1(X(3)-X(2))*X(1)+ipara.dergammaN1(X(3)-X(2));
ipara.der3rN1=@(X) ipara.derbetaN1(X(3)-X(2))*X(1)-ipara.dergammaN1(X(3)-X(2));

%%% P
ipara.derbetaP0=@(x) -3*ipara.bP0*ipara.mP0*exp(-3*ipara.bP0*x) + 3*ipara.aP0*ipara.kP0*exp(3*ipara.aP0*x);
ipara.derbetaP1=@(x) -3*ipara.bP1*ipara.mP1*exp(-3*ipara.bP1*x) + 3*ipara.aP1*ipara.kP1*exp(3*ipara.aP1*x);
ipara.dergammaP0=@(x) -3*ipara.bP0*ipara.mP0*ipara.Pm*exp(-3*ipara.bP0*x);
ipara.dergammaP1=@(x) 3*ipara.aP1*ipara.kP1*ipara.Pm*exp(3*ipara.aP1*x);
ipara.der1rP0=@(X) ipara.betaP0(X(2));
ipara.der2rP0=@(X) ipara.derbetaP0(X(2))*X(1)-ipara.dergammaP0(X(2));
ipara.der1rP1=@(X) ipara.betaP1(X(3)-X(2));
ipara.der2rP1=@(X) -ipara.derbetaP1(X(3)-X(2))*X(1)+ipara.dergammaP1(X(3)-X(2));
ipara.der3rP1=@(X) ipara.derbetaP1(X(3)-X(2))*X(1)-ipara.dergammaP1(X(3)-X(2));

%% BC order 2
%%% C
ipara.der2betaC0=@(x) 1/4*((-2*ipara.bC0)^2*ipara.mC0*exp(-2*ipara.bC0*x) + (2*ipara.aC0)^2*ipara.kC0*exp(2*ipara.aC0*x));
ipara.der2betaC1=@(x) 1/4*((-3*ipara.bC1)^2*ipara.mC1*exp(-3*ipara.bC1*x) + (3*ipara.aC1)^2*ipara.kC1*exp(3*ipara.aC1*x));
ipara.der2gammaC0=@(x) (-2*ipara.bC0)^2*ipara.mC0*exp(-2*ipara.bC0*x);
ipara.der2gammaC1=@(x) (3*ipara.aC1)^2*ipara.kC1*exp(3*ipara.aC1*x);
ipara.der11rC0=@(X) 0;
ipara.der12rC0=@(X) ipara.derbetaC0(X(2));
ipara.der21rC0=@(X) ipara.derbetaC0(X(2));
ipara.der22rC0=@(X) ipara.der2betaC0(X(2))*X(1)-ipara.der2gammaC0(X(2));
ipara.der11rC1=@(X) 0;
ipara.der12rC1=@(X) -ipara.derbetaC1(X(3)-X(2));
ipara.der13rC1=@(X) ipara.derbetaC1(X(3)-X(2));
ipara.der21rC1=@(X) -ipara.derbetaC1(X(3)-X(2));
ipara.der22rC1=@(X) ipara.der2betaC1(X(3)-X(2))*X(1)-ipara.der2gammaC1(X(3)-X(2));
ipara.der23rC1=@(X) -ipara.der2betaC1(X(3)-X(2))*X(1)+ipara.der2gammaC1(X(3)-X(2));
ipara.der31rC1=@(X) ipara.derbetaC1(X(3)-X(2));
ipara.der32rC1=@(X) -ipara.der2betaC1(X(3)-X(2))*X(1)+ipara.der2gammaC1(X(3)-X(2));
ipara.der33rC1=@(X) ipara.der2betaC1(X(3)-X(2))*X(1)-ipara.der2gammaC1(X(3)-X(2));

%%% N
ipara.der2betaN0=@(x) (-ipara.aN0)^2*ipara.kN0*exp(ipara.aN0*x) + (-ipara.br0)^2*ipara.pN0*exp(-ipara.br0*x);
ipara.der2betaN1=@(x) 0;
ipara.der2gammaN0=@(x) (ipara.bN0)^2*ipara.mN0*exp(ipara.bN0*x) + (ipara.ar0)^2*ipara.nN0*exp(ipara.ar0*x);
ipara.der2gammaN1=@(x) ipara.kN1*ipara.Nmetal*exp(-x)/(1+exp(-x))^2;
ipara.der11rN0=@(X) 0;
ipara.der12rN0=@(X) ipara.derbetaN0(X(2));
ipara.der21rN0=@(X) ipara.derbetaN0(X(2));
ipara.der22rN0=@(X) ipara.der2betaN0(X(2))*X(1)-ipara.der2gammaN0(X(2));
ipara.der11rN1=@(X) 0;
ipara.der12rN1=@(X) -ipara.derbetaN1(X(3)-X(2));
ipara.der13rN1=@(X) ipara.derbetaN1(X(3)-X(2));
ipara.der21rN1=@(X) -ipara.derbetaN1(X(3)-X(2));
ipara.der22rN1=@(X) ipara.der2betaN1(X(3)-X(2))*X(1)-ipara.der2gammaN1(X(3)-X(2));
ipara.der23rN1=@(X) -ipara.der2betaN1(X(3)-X(2))*X(1)+ipara.der2gammaN1(X(3)-X(2));
ipara.der31rN1=@(X) ipara.derbetaN1(X(3)-X(2));
ipara.der32rN1=@(X) -ipara.der2betaN1(X(3)-X(2))*X(1)+ipara.der2gammaN1(X(3)-X(2));
ipara.der33rN1=@(X) ipara.der2betaN1(X(3)-X(2))*X(1)-ipara.der2gammaN1(X(3)-X(2));

%%% P
ipara.der2betaP0=@(x) (-3*ipara.bP0)^2*ipara.mP0*exp(-3*ipara.bP0*x) + (3*ipara.aP0)^2*ipara.kP0*exp(3*ipara.aP0*x);
ipara.der2betaP1=@(x) (-3*ipara.bP1)^2*ipara.mP1*exp(-3*ipara.bP1*x) + (3*ipara.aP1)^2*ipara.kP1*exp(3*ipara.aP1*x);
ipara.der2gammaP0=@(x) (-3*ipara.bP0)^2*ipara.mP0*ipara.Pm*exp(-3*ipara.bP0*x);
ipara.der2gammaP1=@(x) (3*ipara.aP1)^2*ipara.kP1*ipara.Pm*exp(3*ipara.aP1*x);
ipara.der11rP0=@(X) 0;
ipara.der12rP0=@(X) ipara.derbetaP0(X(2));
ipara.der21rP0=@(X) ipara.derbetaP0(X(2));
ipara.der22rP0=@(X) ipara.der2betaP0(X(2))*X(1)-ipara.der2gammaP0(X(2));
ipara.der11rP1=@(X) 0;
ipara.der12rP1=@(X) -ipara.derbetaP1(X(3)-X(2));
ipara.der13rP1=@(X) ipara.derbetaP1(X(3)-X(2));
ipara.der21rP1=@(X) -ipara.derbetaP1(X(3)-X(2));
ipara.der22rP1=@(X) ipara.der2betaP1(X(3)-X(2))*X(1)-ipara.der2gammaP1(X(3)-X(2));
ipara.der23rP1=@(X) -ipara.der2betaP1(X(3)-X(2))*X(1)+ipara.der2gammaP1(X(3)-X(2));
ipara.der31rP1=@(X) ipara.derbetaP1(X(3)-X(2));
ipara.der32rP1=@(X) -ipara.der2betaP1(X(3)-X(2))*X(1)+ipara.der2gammaP1(X(3)-X(2));
ipara.der33rP1=@(X) ipara.der2betaP1(X(3)-X(2))*X(1)-ipara.der2gammaP1(X(3)-X(2));


%% chi
ipara.chi=@(x) ipara.OmegaFe/ipara.OmegaOx*ipara.kd0*exp(ipara.rho*ipara.ad0*x);
ipara.derchi=@(x) ipara.rho*ipara.ad0*ipara.OmegaFe/ipara.OmegaOx*ipara.kd0*exp(ipara.rho*ipara.ad0*x);
ipara.der2chi=@(x) (ipara.rho*ipara.ad0)^2*ipara.OmegaFe/ipara.OmegaOx*ipara.kd0*exp(ipara.rho*ipara.ad0*x);