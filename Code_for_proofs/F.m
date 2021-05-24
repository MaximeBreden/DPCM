function FX=F(X,para)

K=(length(X)-5)/5;
V=para.V;
FX=zeros(5*K+5,1);
   
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
if exist('intval','file') && isintval(X(1))
    Int=diag(1./intval(2*(1:K-1)))*spdiags([-ones(K,1),ones(K,1)],[0,2],K-1,K);
    e=intval([1;zeros(K-1,1)]);
    FX=intval(FX);
else
    Int=diag(1./(2*(1:K-1)))*spdiags([-ones(K,1),ones(K,1)],[0,2],K-1,K);
    e=[1;zeros(K-1,1)];
end

%FPsi
FX(1)=BC0*(Psi-2*para.alpha0/l*E)-para.DPsi0;
FX(2:K)=Psi(2:K)+Int*E;

%FE
FX(K+1)=BC1*(Psi+2*para.alpha1/l*E)+para.DPsi1-V;
FX(K+2:2*K)=E(2:K)-para.coupling*l^2/(4*para.lambda^2)*Int*(3*P-N+2*C+para.rho*e);

%FC
FX(2*K+1)=(2*l*para.rC0(BC0*[C Psi])+JC);%/para.kC0;
FX(2*K+2:3*K)=C(2:K)+Int*(-1/4*JC*e-para.zC*convo(C,E)-para.epsC/2*delta*l*C);

%FN
FX(3*K+1)=2*l*para.rN0(BC0*[N Psi])+JN;
FX(3*K+2:4*K)=N(2:K)+Int*(-1/4*JN*e-para.zN*convo(N,E)-para.epsN/2*delta*l*N);

%FP
FX(4*K+1)=(2*l*para.rP0(BC0*[P Psi])+JP);%/para.kP0;
FX(4*K+2:5*K)=P(2:K)+Int*(-1/4*JP*e-para.zP*convo(P,E)-para.epsP/2*delta*l*P);

%FJC
FX(5*K+1)=2*l*para.rC1([BC1*[C Psi] V])-JC;

%FJN
FX(5*K+2)=2*l*para.rN1([BC1*[N Psi] V])-JN;

%FJP
FX(5*K+3)=(2*l*para.rP1([BC1*[P Psi] V])-JP);%/para.kP1;

%Fdelta
FX(5*K+4)=delta-para.chi(BC0*Psi);

%Fl
FX(5*K+5)=l+para.kappa*JC/(2*delta*para.epsC);
