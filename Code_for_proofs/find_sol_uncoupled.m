function [l,JC,JN,JP,delta,C,N,P,Psi_0,Psi_1]=find_sol_uncoupled(para,lmin,lmax,tol,it_max)

% Code from "C. Chainais-Hillairet and T.O. Gallouet. Study of a 
% pseudo-stationary state for a corrosion model: existence and numerical 
% approximation. Nonlinear Anal-Real., 31, 38-56, 2016." which computes a 
% solution without coupling, by applying a dichotomy algorithm on l.

%% DONNEES PHYSIQUES ET ADIMENSIONNEES

delta_Psi0_pzc=para.DPsi0;
delta_Psi1_pzc=para.DPsi1;
Vad=para.V;
alpha_0=para.alpha0;
alpha_1=para.alpha1;
ad_0=para.ad0;
PB_ratio=para.PB;
kappa=para.kappa;
kd_0_ad=para.kd0;
epsilonC=para.epsC;
epsilonN=para.epsN;
epsilonP=para.epsP;

%% FONCTIONS DEFINISSANT LES CONDITIONS AUX LIMITES
betaC0=para.betaC0;
betaC1=para.betaC1;
gammaC0=para.gammaC0;
gammaC1=para.gammaC1;

betaN0=para.betaN0;
betaN1=para.betaN1;
gammaN0=para.gammaN0;
gammaN1=para.gammaN1;

betaP0=para.betaP0;
betaP1=para.betaP1;
gammaP0=para.gammaP0;
gammaP1=para.gammaP1;


%% calcul des coefficients K1tilde et K2tilde

Psitilde=(alpha_1*delta_Psi0_pzc+alpha_0*(Vad-delta_Psi1_pzc))/(alpha_0+alpha_1);

K1tilde=PB_ratio*(kappa/epsilonC*gammaC1(delta_Psi1_pzc)-betaC1(delta_Psi1_pzc)/epsilonC)*exp(5*ad_0*Psitilde);

K2tilde=PB_ratio*kappa/epsilonC*(betaC0(Psitilde)*gammaC1(Vad-Psitilde)-betaC1(Vad-Psitilde)*gammaC0(Psitilde))...
            /(betaC0(Psitilde)+betaC1(Vad-Psitilde))*exp(5*ad_0*Psitilde);


%% Verification de la condition min(K1tilde,K2tilde)<kd_0_ad<max(K1tilde,K2tilde)

test=(kd_0_ad>min(K1tilde,K2tilde))&(kd_0_ad<max(K1tilde,K2tilde));

if not(test)
    disp('There might not be a physical solution with these parameters')
end

%  vall=logspace(-8,8,33);
%  valv=0*vall;
%  N=length(vall);
%  for n=1:N
%      [vn,~,~]=F(vall(n));
%      valv(n)=vn+1;
%  end
%  figure
%  loglog(vall,valv)
 
[vmin,~,~]=F(lmin);
[vmax,~,~]=F(lmax);
if vmin*vmax>0
    error('No possible dichotomy with these choices of lmin and lmax')
end

it=0;
while it<it_max && (lmax-lmin)>tol
    l=(lmin+lmax)/2;
    [v,~,~]=F(l);
    if v*vmin<=0
        lmax=l;
        vmax=v;
    else
        lmin=l;
        vmin=v;
    end
    it=it+1;
end

if (lmax-lmin)>tol
    disp('Maximal number of iteration reached, the dochotomy may not have converged')
end

[~,JC,delta]=F(l);

somlalpha=l+alpha_0+alpha_1;
Psi_0=((l+alpha_1)*delta_Psi0_pzc+alpha_0*(Vad-delta_Psi1_pzc))/somlalpha;
Psi_1=(alpha_1*delta_Psi0_pzc+(l+alpha_0)*(Vad-delta_Psi1_pzc))/somlalpha;
Psi= @(x)(Psi_0+(Psi_1-Psi_0)*(x+1)/2);

%% CALCUL DE C
GammaC1=gammaC1(Vad-Psi_1);
BetaC1=betaC1(Vad-Psi_1);
GammaC0=gammaC0(Psi_0);
BetaC0=betaC0(Psi_0);
dC=2*(Vad-delta_Psi1_pzc-delta_Psi0_pzc)/delta/somlalpha;
num=GammaC1/BetaC1-GammaC0/BetaC0*exp(-(dC+epsilonC)*l*delta);
denom=delta/BetaC0*exp(-(dC+epsilonC)*l*delta)+delta/BetaC1+(1-exp(-(dC+epsilonC)*l*delta))/(dC+epsilonC);
JCsurldelta=-num/denom;

c_0=(GammaC0-delta*JCsurldelta)*exp(2*Psi_0)/BetaC0;
c=@(x)(c_0-JCsurldelta*(exp((dC+epsilonC)*l*delta*(x+1)/2)-1)/(dC+epsilonC)*exp(2*Psi_0));
C=@(x) (c(x).*exp(-2*Psi(x)-epsilonC*l*delta*(x+1)/2));

%% CALCUL DE JN
GammaN1=gammaN1(Vad-Psi_1);
BetaN1=betaN1(Vad-Psi_1);
GammaN0=gammaN0(Psi_0);
BetaN0=betaN0(Psi_0);
dN=-(Vad-delta_Psi1_pzc-delta_Psi0_pzc)/delta/somlalpha;

numN=GammaN1/BetaN1-GammaN0/BetaN0*exp(-(dN+epsilonN)*l*delta);
denomN=delta/BetaN0*exp(-(dN+epsilonN)*l*delta)+delta/BetaN1+(1-exp(-(dN+epsilonN)*l*delta))/(dN+epsilonN);
JNsurldelta=-numN/denomN;
JN=2*JNsurldelta*l*delta; %2 à cause de [0,1] VS [-1,1]

%% CALCUL DE N
n_0=(GammaN0-delta*JNsurldelta)*exp(-Psi_0)/BetaN0;
n=@(x)(n_0-JNsurldelta*(exp((dN+epsilonN)*l*delta*(x+1)/2)-1)/(dN+epsilonC)*exp(-Psi_0));
N=@(x) (n(x).*exp(Psi(x)-epsilonN*l*delta*(x+1)/2));

%% CALCUL DE JP
GammaP1=gammaP1(Vad-Psi_1);
BetaP1=betaP1(Vad-Psi_1);
GammaP0=gammaP0(Psi_0);
BetaP0=betaP0(Psi_0);
dP=3*(Vad-delta_Psi1_pzc-delta_Psi0_pzc)/delta/somlalpha;

numP=GammaP1/BetaP1-GammaP0/BetaP0*exp(-(dP+epsilonP)*l*delta);
denomP=delta/BetaP0*exp(-(dP+epsilonP)*l*delta)+delta/BetaP1+(1-exp(-(dP+epsilonP)*l*delta))/(dP+epsilonP);
JPsurldelta=-numP/denomP;
JP=2*JPsurldelta*l*delta; %2 à cause de [0,1] VS [-1,1]

%% CALCUL DE P
p_0=(GammaP0-delta*JPsurldelta)*exp(3*Psi_0)/BetaP0;
p=@(x)(p_0-JPsurldelta*(exp((dP+epsilonP)*l*delta*(x+1)/2)-1)/(dP+epsilonP)*exp(3*Psi_0));
P=@(x) (p(x).*exp(-3*Psi(x)-epsilonP*l*delta*(x+1)/2));


function [val,JC,delta]=F(l)
% val=0 ssi pt fixe 

%% CALCUL DE Psi_0 et Psi_1
somlalpha=l+alpha_0+alpha_1;
Psi_0=((l+alpha_1)*delta_Psi0_pzc+alpha_0*(Vad-delta_Psi1_pzc))/somlalpha;
Psi_1=(alpha_1*delta_Psi0_pzc+(l+alpha_0)*(Vad-delta_Psi1_pzc))/somlalpha;

%% CALCUL DE delta
delta=kd_0_ad/PB_ratio*exp(-5*ad_0*Psi_0);

%% CALCUL DE JC
GammaC1=gammaC1(Vad-Psi_1);
BetaC1=betaC1(Vad-Psi_1);
GammaC0=gammaC0(Psi_0);
BetaC0=betaC0(Psi_0);
dC=2*(Vad-delta_Psi1_pzc-delta_Psi0_pzc)/delta/somlalpha;

num=GammaC1/BetaC1-GammaC0/BetaC0*exp(-(dC+epsilonC)*l*delta);
denom=delta/BetaC0*exp(-(dC+epsilonC)*l*delta)+delta/BetaC1+(1-exp(-(dC+epsilonC)*l*delta))/(dC+epsilonC);
JCsurldelta=-num/denom;
JC=2*JCsurldelta*l*delta; %2 à cause de [0,1] VS [-1,1]

val=-kappa/epsilonC*JCsurldelta-1;
end

end