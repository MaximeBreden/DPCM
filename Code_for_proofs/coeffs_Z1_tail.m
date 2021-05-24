function coeffs=coeffs_Z1_tail(X,para,nu,A)

if strcmp(para.case,'Poten')
    K=(length(X)-5)/5;
    V=para.V;
else
    error('Invalid value of para.case')
end

Psi=X(1:K);
E=X(K+1:2*K);
C=X(2*K+1:3*K);
N=X(3*K+1:4*K);
P=X(4*K+1:5*K);
JC=X(5*K+1);
JN=X(5*K+2);
JP=X(5*K+3);
delta=X(5*K+4);
l=X(5*K+5);

BC0=[1 2*(-1).^(1:K-1)];
BC1=[1 2*ones(1,K-1)];
e=[1;zeros(K-1,1)];
if exist('intval','file') && isintval(X(1))
    e=intval(e);
end

ind_BC=[1,K+1,2*K+1,3*K+1,4*K+1,5*K+1,5*K+2,5*K+3,5*K+4,5*K+5];

%ADFDPsi
coefBC_DPsi=[-2,2,-2*2*l*para.der2rC0(BC0*[C Psi]),-2*2*l*para.der2rN0(BC0*[N Psi]),-2*2*l*para.der2rP0(BC0*[P Psi]),2*l*para.der2rC1([BC1*[C Psi] V]),2*l*para.der2rN1([BC1*[N Psi] V]),2*l*para.der2rP1([BC1*[P Psi] V]),-2*(-para.derchi(BC0*Psi)),0;
             2,2,2*2*l*para.der2rC0(BC0*[C Psi]),2*2*l*para.der2rN0(BC0*[N Psi]),2*2*l*para.der2rP0(BC0*[P Psi]),2*l*para.der2rC1([BC1*[C Psi] V]),2*l*para.der2rN1([BC1*[N Psi] V]),2*l*para.der2rP1([BC1*[P Psi] V]),2*(-para.derchi(BC0*Psi)),0];
ADFPsiDPsi=max_norm(A(1:K,ind_BC)*coefBC_DPsi',nu)/(2*nu^(2*K));
ADFEDPsi=max_norm(A(K+1:2*K,ind_BC)*coefBC_DPsi',nu)/(2*nu^(2*K));
ADFCDPsi=max_norm(A(2*K+1:3*K,ind_BC)*coefBC_DPsi',nu)/(2*nu^(2*K));
ADFNDPsi=max_norm(A(3*K+1:4*K,ind_BC)*coefBC_DPsi',nu)/(2*nu^(2*K));
ADFPDPsi=max_norm(A(4*K+1:5*K,ind_BC)*coefBC_DPsi',nu)/(2*nu^(2*K));
ADFJCDPsi=max(A(5*K+1,ind_BC)*coefBC_DPsi')/(2*nu^(2*K));
ADFJNDPsi=max(A(5*K+2,ind_BC)*coefBC_DPsi')/(2*nu^(2*K));
ADFJPDPsi=max(A(5*K+3,ind_BC)*coefBC_DPsi')/(2*nu^(2*K));
ADFdeltaDPsi=max(A(5*K+4,ind_BC)*coefBC_DPsi')/(2*nu^(2*K));
ADFlDPsi=max(A(5*K+5,ind_BC)*coefBC_DPsi')/(2*nu^(2*K));

%ADFDE
coefBC_DE=[-2*(-2*para.alpha0/l),2*2*para.alpha1/l,0,0,0,0,0,0,0,0;
           2*(-2*para.alpha0/l),2*2*para.alpha1/l,0,0,0,0,0,0,0,0];
ADFPsiDE=max_norm(A(1:K,ind_BC)*coefBC_DE',nu)/(2*nu^(2*K))+1/2*(nu^(-1)/(2*K-1)+nu/(2*K+1));
ADFEDE=max_norm(A(K+1:2*K,ind_BC)*coefBC_DE',nu)/(2*nu^(2*K));
ADFCDE=max_norm(A(2*K+1:3*K,ind_BC)*coefBC_DE',nu)/(2*nu^(2*K))+abs(para.zC)*shift_norm(C,nu)/(2*K);
ADFNDE=max_norm(A(3*K+1:4*K,ind_BC)*coefBC_DE',nu)/(2*nu^(2*K))+abs(para.zN)*shift_norm(N,nu)/(2*K);
ADFPDE=max_norm(A(4*K+1:5*K,ind_BC)*coefBC_DE',nu)/(2*nu^(2*K))+abs(para.zP)*shift_norm(P,nu)/(2*K);
ADFJCDE=max(A(5*K+1,ind_BC)*coefBC_DE')/(2*nu^(2*K));
ADFJNDE=max(A(5*K+2,ind_BC)*coefBC_DE')/(2*nu^(2*K));
ADFJPDE=max(A(5*K+3,ind_BC)*coefBC_DE')/(2*nu^(2*K));
ADFdeltaDE=max(A(5*K+4,ind_BC)*coefBC_DE')/(2*nu^(2*K));
ADFlDE=max(A(5*K+5,ind_BC)*coefBC_DE')/(2*nu^(2*K));

%ADFDC
coefBC_DC=[0,0,-2*2*l*para.der1rC0(BC0*[C Psi]),0,0,2*2*l*para.der1rC1([BC1*[C Psi] V]),0,0,0,0;
           0,0,2*2*l*para.der1rC0(BC0*[C Psi]),0,0,2*2*l*para.der1rC1([BC1*[C Psi] V]),0,0,0,0];
ADFPsiDC=max_norm(A(1:K,ind_BC)*coefBC_DC',nu)/(2*nu^(2*K));
ADFEDC=max_norm(A(K+1:2*K,ind_BC)*coefBC_DC',nu)/(2*nu^(2*K))+para.coupling*l^2/(4*para.lambda^2)*2*1/2*(nu^(-1)/(2*K-1)+nu/(2*K+1));
ADFCDC=max_norm(A(2*K+1:3*K,ind_BC)*coefBC_DC',nu)/(2*nu^(2*K))+shift_norm(-para.zC*E-para.epsC/2*delta*l*e,nu)/(2*K);
ADFNDC=max_norm(A(3*K+1:4*K,ind_BC)*coefBC_DC',nu)/(2*nu^(2*K));
ADFPDC=max_norm(A(4*K+1:5*K,ind_BC)*coefBC_DC',nu)/(2*nu^(2*K));
ADFJCDC=max(A(5*K+1,ind_BC)*coefBC_DC')/(2*nu^(2*K));
ADFJNDC=max(A(5*K+2,ind_BC)*coefBC_DC')/(2*nu^(2*K));
ADFJPDC=max(A(5*K+3,ind_BC)*coefBC_DC')/(2*nu^(2*K));
ADFdeltaDC=max(A(5*K+4,ind_BC)*coefBC_DC')/(2*nu^(2*K));
ADFlDC=max(A(5*K+5,ind_BC)*coefBC_DC')/(2*nu^(2*K));

%ADFDN
coefBC_DN=[0,0,0,-2*2*l*para.der1rN0(BC0*[N Psi]),0,0,2*2*l*para.der1rN1([BC1*[N Psi] V]),0,0,0;
           0,0,0,2*2*l*para.der1rN0(BC0*[N Psi]),0,0,2*2*l*para.der1rN1([BC1*[N Psi] V]),0,0,0];
ADFPsiDN=max_norm(A(1:K,ind_BC)*coefBC_DN',nu)/(2*nu^(2*K));
ADFEDN=max_norm(A(K+1:2*K,ind_BC)*coefBC_DN',nu)/(2*nu^(2*K))+para.coupling*l^2/(4*para.lambda^2)*1/2*(nu^(-1)/(2*K-1)+nu/(2*K+1));
ADFCDN=max_norm(A(2*K+1:3*K,ind_BC)*coefBC_DN',nu)/(2*nu^(2*K));
ADFNDN=max_norm(A(3*K+1:4*K,ind_BC)*coefBC_DN',nu)/(2*nu^(2*K))+shift_norm(-para.zN*E-para.epsN/2*delta*l*e,nu)/(2*K);
ADFPDN=max_norm(A(4*K+1:5*K,ind_BC)*coefBC_DN',nu)/(2*nu^(2*K));
ADFJCDN=max(A(5*K+1,ind_BC)*coefBC_DN')/(2*nu^(2*K));
ADFJNDN=max(A(5*K+2,ind_BC)*coefBC_DN')/(2*nu^(2*K));
ADFJPDN=max(A(5*K+3,ind_BC)*coefBC_DN')/(2*nu^(2*K));
ADFdeltaDN=max(A(5*K+4,ind_BC)*coefBC_DN')/(2*nu^(2*K));
ADFlDN=max(A(5*K+5,ind_BC)*coefBC_DN')/(2*nu^(2*K));

%ADFDP
coefBC_DP=[0,0,0,0,-2*2*l*para.der1rP0(BC0*[P Psi]),0,0,2*2*l*para.der1rP1([BC1*[P Psi] V]),0,0;
           0,0,0,0,2*2*l*para.der1rP0(BC0*[P Psi]),0,0,2*2*l*para.der1rP1([BC1*[P Psi] V]),0,0];
ADFPsiDP=max_norm(A(1:K,ind_BC)*coefBC_DP',nu)/(2*nu^(2*K));
ADFEDP=max_norm(A(K+1:2*K,ind_BC)*coefBC_DP',nu)/(2*nu^(2*K))+para.coupling*l^2/(4*para.lambda^2)*3*1/2*(nu^(-1)/(2*K-1)+nu/(2*K+1));
ADFCDP=max_norm(A(2*K+1:3*K,ind_BC)*coefBC_DP',nu)/(2*nu^(2*K));
ADFNDP=max_norm(A(3*K+1:4*K,ind_BC)*coefBC_DP',nu)/(2*nu^(2*K));
ADFPDP=max_norm(A(4*K+1:5*K,ind_BC)*coefBC_DP',nu)/(2*nu^(2*K))+shift_norm(-para.zP*E-para.epsP/2*delta*l*e,nu)/(2*K);
ADFJCDP=max(A(5*K+1,ind_BC)*coefBC_DP')/(2*nu^(2*K));
ADFJNDP=max(A(5*K+2,ind_BC)*coefBC_DP')/(2*nu^(2*K));
ADFJPDP=max(A(5*K+3,ind_BC)*coefBC_DP')/(2*nu^(2*K));
ADFdeltaDP=max(A(5*K+4,ind_BC)*coefBC_DP')/(2*nu^(2*K));
ADFlDP=max(A(5*K+5,ind_BC)*coefBC_DP')/(2*nu^(2*K));

%ADFDJC
ADFPsiDJC=0;
ADFEDJC=0;
ADFCDJC=0;
ADFNDJC=0;
ADFPDJC=0;
ADFJCDJC=0;
ADFJNDJC=0;
ADFJPDJC=0;
ADFdeltaDJC=0;
ADFlDJC=0;

%ADFDJN
ADFPsiDJN=0;
ADFEDJN=0;
ADFCDJN=0;
ADFNDJN=0;
ADFPDJN=0;
ADFJCDJN=0;
ADFJNDJN=0;
ADFJPDJN=0;
ADFdeltaDJN=0;
ADFlDJN=0;

%ADFDJP
ADFPsiDJP=0;
ADFEDJP=0;
ADFCDJP=0;
ADFNDJP=0;
ADFPDJP=0;
ADFJCDJP=0;
ADFJNDJP=0;
ADFJPDJP=0;
ADFdeltaDJP=0;
ADFlDJP=0;

%ADFDdelta
ADFPsiDdelta=0;
ADFEDdelta=0;
ADFCDdelta=0;
ADFNDdelta=0;
ADFPDdelta=0;
ADFJCDdelta=0;
ADFJNDdelta=0;
ADFJPDdelta=0;
ADFdeltaDdelta=0;
ADFlDdelta=0;

%ADFDl
ADFPsiDl=0;
ADFEDl=0;
ADFCDl=0;
ADFNDl=0;
ADFPDl=0;
ADFJCDl=0;
ADFJNDl=0;
ADFJPDl=0;
ADFdeltaDl=0;
ADFlDl=0;

coeffs=[ADFPsiDPsi,ADFPsiDE,ADFPsiDC,ADFPsiDN,ADFPsiDP,ADFPsiDJC,ADFPsiDJN,ADFPsiDJP,ADFPsiDdelta,ADFPsiDl;
        ADFEDPsi,ADFEDE,ADFEDC,ADFEDN,ADFEDP,ADFEDJC,ADFEDJN,ADFEDJP,ADFEDdelta,ADFEDl;
        ADFCDPsi,ADFCDE,ADFCDC,ADFCDN,ADFCDP,ADFCDJC,ADFCDJN,ADFCDJP,ADFCDdelta,ADFCDl;
        ADFNDPsi,ADFNDE,ADFNDC,ADFNDN,ADFNDP,ADFNDJC,ADFNDJN,ADFNDJP,ADFNDdelta,ADFNDl;
        ADFPDPsi,ADFPDE,ADFPDC,ADFPDN,ADFPDP,ADFPDJC,ADFPDJN,ADFPDJP,ADFPDdelta,ADFPDl;
        ADFJCDPsi,ADFJCDE,ADFJCDC,ADFJCDN,ADFJCDP,ADFJCDJC,ADFJCDJN,ADFJCDJP,ADFJCDdelta,ADFJCDl;
        ADFJNDPsi,ADFJNDE,ADFJNDC,ADFJNDN,ADFJNDP,ADFJNDJC,ADFJNDJN,ADFJNDJP,ADFJNDdelta,ADFJNDl;
        ADFJPDPsi,ADFJPDE,ADFJPDC,ADFJPDN,ADFJPDP,ADFJPDJC,ADFJPDJN,ADFJPDJP,ADFJPDdelta,ADFJPDl;
        ADFdeltaDPsi,ADFdeltaDE,ADFdeltaDC,ADFdeltaDN,ADFdeltaDP,ADFdeltaDJC,ADFdeltaDJN,ADFdeltaDJP,ADFdeltaDdelta,ADFdeltaDl;
        ADFlDPsi,ADFlDE,ADFlDC,ADFlDN,ADFlDP,ADFlDJC,ADFlDJN,ADFlDJP,ADFlDdelta,ADFlDl];


function val=max_norm(v,nu)
val=max([1 2*nu.^(1:K-1)]*abs(v));
end

function val=shift_norm(u,nu)
u_shift=[u(2:end);0;0]-[u(2);u];
val=[1 (nu.^(1:K)+nu.^(-(1:K)))]*abs(u_shift);
end

end
