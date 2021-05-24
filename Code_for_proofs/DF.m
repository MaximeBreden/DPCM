function DFX=DF(X,para)

K=(length(X)-5)/5;
V=para.V;

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
Id=spdiags(ones(K,1),1,K-1,K);
if exist('intval','file') && isintval(X(1))
    Int=diag(1./intval(2*(1:K-1)))*spdiags([-ones(K,1),ones(K,1)],[0,2],K-1,K);
    ZerosKm11=intval(zeros(K-1,1));
    ZerosKm1K=intval(zeros(K-1,K));
    Zeros1K=intval(zeros(1,K));
    ZerosKK=intval(zeros(K,K));
    ZerosK1=intval(zeros(K,1));
else
    Int=diag(1./(2*(1:K-1)))*spdiags([-ones(K,1),ones(K,1)],[0,2],K-1,K);
    ZerosKm11=zeros(K-1,1);
    ZerosKm1K=zeros(K-1,K);
    Zeros1K=zeros(1,K);
    ZerosKK=zeros(K,K);
    ZerosK1=zeros(K,1);
end
e=[1;ZerosKm11];
    
Me=convomat(e);
MC=convomat(C);
MN=convomat(N);
MP=convomat(P);
ME=convomat(E);


%DFDPsi
DFPsiDPsi=[BC0;Id];
DFEDPsi=[BC1;zeros(K-1,K)];
DFCDPsi=[2*l*para.der2rC0(BC0*[C Psi])*BC0;ZerosKm1K];
DFNDPsi=[2*l*para.der2rN0(BC0*[N Psi])*BC0;ZerosKm1K];
DFPDPsi=[2*l*para.der2rP0(BC0*[P Psi])*BC0;ZerosKm1K];
DFJCDPsi=2*l*para.der2rC1([BC1*[C Psi] V])*BC1;
DFJNDPsi=2*l*para.der2rN1([BC1*[N Psi] V])*BC1;
DFJPDPsi=2*l*para.der2rP1([BC1*[P Psi] V])*BC1;
DFdeltaDPsi=-para.derchi(BC0*Psi)*BC0;
DFlDPsi=Zeros1K;

%DFDE
DFPsiDE=[-2*para.alpha0/l*BC0;Int];
DFEDE=[2*para.alpha1/l*BC1;Id];
DFCDE=[Zeros1K;Int*(-para.zC*MC)];
DFNDE=[Zeros1K;Int*(-para.zN*MN)];
DFPDE=[Zeros1K;Int*(-para.zP*MP)];
DFJCDE=Zeros1K;
DFJNDE=Zeros1K;
DFJPDE=Zeros1K;
DFdeltaDE=Zeros1K;
DFlDE=Zeros1K;

%DFDC
DFPsiDC=zeros(K,K);
DFEDC=[Zeros1K;-para.coupling*l^2/(4*para.lambda^2)*2*Int];
DFCDC=[2*l*para.der1rC0(BC0*[C Psi])*BC0;Id+Int*(-para.zC*ME-para.epsC/2*delta*l*Me)];
DFNDC=zeros(K,K);
DFPDC=zeros(K,K);
DFJCDC=2*l*para.der1rC1([BC1*[C Psi] V])*BC1;
DFJNDC=Zeros1K;
DFJPDC=Zeros1K;
DFdeltaDC=Zeros1K;
DFlDC=Zeros1K;

%DFDN
DFPsiDN=ZerosKK;
DFEDN=[Zeros1K;-para.coupling*l^2/(4*para.lambda^2)*(-1)*Int];
DFCDN=ZerosKK;
DFNDN=[2*l*para.der1rN0(BC0*[N Psi])*BC0;Id+Int*(-para.zN*ME-para.epsN/2*delta*l*Me)];
DFPDN=ZerosKK;
DFJCDN=Zeros1K;
DFJNDN=2*l*para.der1rN1([BC1*[N Psi] V])*BC1;
DFJPDN=Zeros1K;
DFdeltaDN=Zeros1K;
DFlDN=Zeros1K;

%DFDP
DFPsiDP=ZerosKK;
DFEDP=[Zeros1K;-para.coupling*l^2/(4*para.lambda^2)*3*Int];
DFCDP=ZerosKK;
DFNDP=ZerosKK;
DFPDP=[2*l*para.der1rP0(BC0*[P Psi])*BC0;Id+Int*(-para.zP*ME-para.epsP/2*delta*l*Me)];
DFJCDP=Zeros1K;
DFJNDP=Zeros1K;
DFJPDP=2*l*para.der1rP1([BC1*[P Psi] V])*BC1;
DFdeltaDP=Zeros1K;
DFlDP=Zeros1K;

%DFDJC
DFPsiDJC=ZerosK1;
DFEDJC=ZerosK1;
DFCDJC=[1;Int*(-1/4*e)];
DFNDJC=ZerosK1;
DFPDJC=ZerosK1;
DFJCDJC=-1;
DFJNDJC=0;
DFJPDJC=0;
DFdeltaDJC=0;
DFlDJC=para.kappa/(2*delta*para.epsC);

%DFDJN
DFPsiDJN=ZerosK1;
DFEDJN=ZerosK1;
DFCDJN=ZerosK1;
DFNDJN=[1;Int*(-1/4*e)];
DFPDJN=ZerosK1;
DFJCDJN=0;
DFJNDJN=-1;
DFJPDJN=0;
DFdeltaDJN=0;
DFlDJN=0;

%DFDJP
DFPsiDJP=ZerosK1;
DFEDJP=ZerosK1;
DFCDJP=ZerosK1;
DFNDJP=ZerosK1;
DFPDJP=[1;Int*(-1/4*e)];
DFJCDJP=0;
DFJNDJP=0;
DFJPDJP=-1;
DFdeltaDJP=0;
DFlDJP=0;

%DFDdelta
DFPsiDdelta=ZerosK1;
DFEDdelta=ZerosK1;
DFCDdelta=[0;Int*(-para.epsC/2*l*C)];
DFNDdelta=[0;Int*(-para.epsN/2*l*N)];
DFPDdelta=[0;Int*(-para.epsP/2*l*P)];
DFJCDdelta=0;
DFJNDdelta=0;
DFJPDdelta=0;
DFdeltaDdelta=1;
DFlDdelta=-para.kappa*JC/(2*delta^2*para.epsC);

%DFDl
DFPsiDl=[2*para.alpha0/l^2*BC0*E;ZerosKm11];
DFEDl=[-2*para.alpha1/l^2*BC1*E;-para.coupling*l/(2*para.lambda^2)*Int*(3*P-N+2*C+para.rho*e)];
DFCDl=[2*para.rC0(BC0*[C Psi]);Int*(-para.epsC/2*delta*C)];
DFNDl=[2*para.rN0(BC0*[N Psi]);Int*(-para.epsN/2*delta*N)];
DFPDl=[2*para.rP0(BC0*[P Psi]);Int*(-para.epsP/2*delta*P)];
DFJCDl=2*para.rC1([BC1*[C Psi] V]);
DFJNDl=2*para.rN1([BC1*[N Psi] V]);
DFJPDl=2*para.rP1([BC1*[P Psi] V]);
DFdeltaDl=0;
DFlDl=1;


%DFX
DFX=[DFPsiDPsi,DFPsiDE,DFPsiDC,DFPsiDN,DFPsiDP,DFPsiDJC,DFPsiDJN,DFPsiDJP,DFPsiDdelta,DFPsiDl;
     DFEDPsi,DFEDE,DFEDC,DFEDN,DFEDP,DFEDJC,DFEDJN,DFEDJP,DFEDdelta,DFEDl;
     DFCDPsi,DFCDE,DFCDC,DFCDN,DFCDP,DFCDJC,DFCDJN,DFCDJP,DFCDdelta,DFCDl;
     DFNDPsi,DFNDE,DFNDC,DFNDN,DFNDP,DFNDJC,DFNDJN,DFNDJP,DFNDdelta,DFNDl;
     DFPDPsi,DFPDE,DFPDC,DFPDN,DFPDP,DFPDJC,DFPDJN,DFPDJP,DFPDdelta,DFPDl;
     DFJCDPsi,DFJCDE,DFJCDC,DFJCDN,DFJCDP,DFJCDJC,DFJCDJN,DFJCDJP,DFJCDdelta,DFJCDl;
     DFJNDPsi,DFJNDE,DFJNDC,DFJNDN,DFJNDP,DFJNDJC,DFJNDJN,DFJNDJP,DFJNDdelta,DFJNDl;
     DFJPDPsi,DFJPDE,DFJPDC,DFJPDN,DFJPDP,DFJPDJC,DFJPDJN,DFJPDJP,DFJPDdelta,DFJPDl;
     DFdeltaDPsi,DFdeltaDE,DFdeltaDC,DFdeltaDN,DFdeltaDP,DFdeltaDJC,DFdeltaDJN,DFdeltaDJP,DFdeltaDdelta,DFdeltaDl;
     DFlDPsi,DFlDE,DFlDC,DFlDN,DFlDP,DFlDJC,DFlDJN,DFlDJP,DFlDdelta,DFlDl];

