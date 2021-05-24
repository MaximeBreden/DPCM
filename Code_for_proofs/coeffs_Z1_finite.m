function coeffs=coeffs_Z1_finite(X,para,nu,A)

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

ZerosKm11=zeros(K-1,1);
ZerosK1=zeros(K,1);
Zeros2K1=zeros(2*K,1);
Zeros2Km12K=zeros(2*K-1,2*K);
Zeros2K2K=zeros(2*K,2*K);
Zeros12K=zeros(1,2*K);
    
if exist('intval','file') && isintval(X(1))
    Int=diag(1./intval(2*(1:2*K-1)))*spdiags([-ones(2*K,1),ones(2*K,1)],[0,2],2*K-1,2*K);
    ZerosKm11=intval(ZerosKm11);
    ZerosK1=intval(ZerosK1);
    Zeros2K1=intval(Zeros2K1);
    Zeros2Km12K=intval(Zeros2Km12K);
    Zeros2K2K=intval(Zeros2K2K);
    Zeros12K=intval(Zeros12K);
else
    Int=diag(1./(2*(1:2*K-1)))*spdiags([-ones(2*K,1),ones(2*K,1)],[0,2],2*K-1,2*K);
end

e=[1;ZerosKm11];

BC0=[1 2*(-1).^(1:K-1)];
BC1=[1 2*ones(1,K-1)];
BC0_trunc=[zeros(1,K) 2*(-1).^(K:2*K-1)];
BC1_trunc=[zeros(1,K) 2*ones(1,K)];

Trunc_vect=[zeros(K-1,1);ones(K,1)];
Trunc_mat=ones(2*K-1,2*K);
Trunc_mat(1:K-1,1:K)=0;

e_ext=[e;ZerosK1];
E_ext=[E;ZerosK1];
C_ext=[C;ZerosK1];
N_ext=[N;ZerosK1];
P_ext=[P;ZerosK1];
Me_ext=convomat(e_ext);
ME_ext=convomat(E_ext);
MC_ext=convomat(C_ext);
MN_ext=convomat(N_ext);
MP_ext=convomat(P_ext);

%DFDPsi-Adag
DFPsiDPsi=[BC0_trunc;Zeros2Km12K];
DFEDPsi=[BC1_trunc;Zeros2Km12K];
DFCDPsi=[2*l*para.der2rC0(BC0*[C Psi])*BC0_trunc;Zeros2Km12K];
DFNDPsi=[2*l*para.der2rN0(BC0*[N Psi])*BC0_trunc;Zeros2Km12K];
DFPDPsi=[2*l*para.der2rP0(BC0*[P Psi])*BC0_trunc;Zeros2Km12K];
DFJCDPsi=2*l*para.der2rC1([BC1*[C Psi] V])*BC1_trunc;
DFJNDPsi=2*l*para.der2rC1([BC1*[N Psi] V])*BC1_trunc;
DFJPDPsi=2*l*para.der2rC1([BC1*[P Psi] V])*BC1_trunc;
DFdeltaDPsi=-para.derchi(BC0*Psi)*BC0_trunc;
DFlDPsi=Zeros12K;

DFDPsi=[DFPsiDPsi;DFEDPsi;DFCDPsi;DFNDPsi;DFPDPsi;DFJCDPsi;DFJNDPsi;DFJPDPsi;DFdeltaDPsi;DFlDPsi];
ADFDPsi=prod_ext(A,DFDPsi,K);
ADFDPsi_norm=block2norm(ADFDPsi,nu);

%DFDE
DFPsiDE=[-2*para.alpha0/l*BC0_trunc;(Int*Me_ext).*Trunc_mat];
DFEDE=[2*para.alpha1/l*BC1_trunc;Zeros2Km12K];
DFCDE=[Zeros12K;(Int*(-para.zC*MC_ext)).*Trunc_mat];
DFNDE=[Zeros12K;(Int*(-para.zN*MN_ext)).*Trunc_mat];
DFPDE=[Zeros12K;(Int*(-para.zP*MP_ext)).*Trunc_mat];
DFJCDE=Zeros12K;
DFJNDE=Zeros12K;
DFJPDE=Zeros12K;
DFdeltaDE=Zeros12K;
DFlDE=Zeros12K;

DFDE=[DFPsiDE;DFEDE;DFCDE;DFNDE;DFPDE;DFJCDE;DFJNDE;DFJPDE;DFdeltaDE;DFlDE];
ADFDE=prod_ext(A,DFDE,K);
ADFDE_norm=block2norm(ADFDE,nu);

%DFDC
DFPsiDC=Zeros2K2K;
DFEDC=[Zeros12K;(-para.coupling*l^2/(4*para.lambda^2)*Int*2).*Trunc_mat];
DFCDC=[2*l*para.der1rC0(BC0*[C Psi])*BC0_trunc;(Int*(-para.zC*ME_ext-para.epsC/2*delta*l*Me_ext)).*Trunc_mat];
DFNDC=Zeros2K2K;
DFPDC=Zeros2K2K;
DFJCDC=2*l*para.der1rC1([BC1*[C Psi] V])*BC1_trunc;
DFJNDC=Zeros12K;
DFJPDC=Zeros12K;
DFdeltaDC=Zeros12K;
DFlDC=Zeros12K;

DFDC=[DFPsiDC;DFEDC;DFCDC;DFNDC;DFPDC;DFJCDC;DFJNDC;DFJPDC;DFdeltaDC;DFlDC];
ADFDC=prod_ext(A,DFDC,K);
ADFDC_norm=block2norm(ADFDC,nu);

%DFDN
DFPsiDN=Zeros2K2K;
DFEDN=[Zeros12K;(-para.coupling*l^2/(4*para.lambda^2)*Int*(-1)).*Trunc_mat];
DFCDN=Zeros2K2K;
DFNDN=[2*l*para.der1rN0(BC0*[N Psi])*BC0_trunc;(Int*(-para.zN*ME_ext-para.epsN/2*delta*l*Me_ext)).*Trunc_mat];
DFPDN=Zeros2K2K;
DFJCDN=Zeros12K;
DFJNDN=2*l*para.der1rN1([BC1*[N Psi] V])*BC1_trunc;
DFJPDN=Zeros12K;
DFdeltaDN=Zeros12K;
DFlDN=Zeros12K;

DFDN=[DFPsiDN;DFEDN;DFCDN;DFNDN;DFPDN;DFJCDN;DFJNDN;DFJPDN;DFdeltaDN;DFlDN];
ADFDN=prod_ext(A,DFDN,K);
ADFDN_norm=block2norm(ADFDN,nu);

%DFDP
DFPsiDP=Zeros2K2K;
DFEDP=[Zeros12K;(-para.coupling*l^2/(4*para.lambda^2)*Int*3).*Trunc_mat];
DFCDP=Zeros2K2K;
DFNDP=Zeros2K2K;
DFPDP=[2*l*para.der1rP0(BC0*[P Psi])*BC0_trunc;(Int*(-para.zP*ME_ext-para.epsP/2*delta*l*Me_ext)).*Trunc_mat];
DFJCDP=Zeros12K;
DFJNDP=Zeros12K;
DFJPDP=2*l*para.der1rP1([BC1*[P Psi] V])*BC1_trunc;
DFdeltaDP=Zeros12K;
DFlDP=Zeros12K;

DFDP=[DFPsiDP;DFEDP;DFCDP;DFNDP;DFPDP;DFJCDP;DFJNDP;DFJPDP;DFdeltaDP;DFlDP];
ADFDP=prod_ext(A,DFDP,K);
ADFDP_norm=block2norm(ADFDP,nu);

%DFDJC
DFPsiDJC=Zeros2K1;
DFEDJC=Zeros2K1;
DFCDJC=[0;(Int*(-1/4*e_ext)).*Trunc_vect];
DFNDJC=Zeros2K1;
DFPDJC=Zeros2K1;
DFJCDJC=0;
DFJNDJC=0;
DFJPDJC=0;
DFdeltaDJC=0;
DFlDJC=0;

DFDJC=[DFPsiDJC;DFEDJC;DFCDJC;DFNDJC;DFPDJC;DFJCDJC;DFJNDJC;DFJPDJC;DFdeltaDJC;DFlDJC];
ADFDJC=prod_ext(A,DFDJC,K);
ADFDJC_norm=block2norm(ADFDJC,nu);

%DFDJN
DFPsiDJN=Zeros2K1;
DFEDJN=Zeros2K1;
DFCDJN=Zeros2K1;
DFNDJN=[0;(Int*(-1/4*e_ext)).*Trunc_vect];
DFPDJN=Zeros2K1;
DFJCDJN=0;
DFJNDJN=0;
DFJPDJN=0;
DFdeltaDJN=0;
DFlDJN=0;

DFDJN=[DFPsiDJN;DFEDJN;DFCDJN;DFNDJN;DFPDJN;DFJCDJN;DFJNDJN;DFJPDJN;DFdeltaDJN;DFlDJN];
ADFDJN=prod_ext(A,DFDJN,K);
ADFDJN_norm=block2norm(ADFDJN,nu);

%DFDJP
DFPsiDJP=Zeros2K1;
DFEDJP=Zeros2K1;
DFCDJP=Zeros2K1;
DFNDJP=Zeros2K1;
DFPDJP=[0;(Int*(-1/4*e_ext)).*Trunc_vect];
DFJCDJP=0;
DFJNDJP=0;
DFJPDJP=0;
DFdeltaDJP=0;
DFlDJP=0;

DFDJP=[DFPsiDJP;DFEDJP;DFCDJP;DFNDJP;DFPDJP;DFJCDJP;DFJNDJP;DFJPDJP;DFdeltaDJP;DFlDJP];
ADFDJP=prod_ext(A,DFDJP,K);
ADFDJP_norm=block2norm(ADFDJP,nu);

%DFDdelta
DFPsiDdelta=Zeros2K1;
DFEDdelta=Zeros2K1;
DFCDdelta=[0;(Int*(-para.epsC/2*l*C_ext)).*Trunc_vect];
DFNDdelta=[0;(Int*(-para.epsN/2*l*N_ext)).*Trunc_vect];
DFPDdelta=[0;(Int*(-para.epsP/2*l*P_ext)).*Trunc_vect];
DFJCDdelta=0;
DFJNDdelta=0;
DFJPDdelta=0;
DFdeltaDdelta=0;
DFlDdelta=0;

DFDdelta=[DFPsiDdelta;DFEDdelta;DFCDdelta;DFNDdelta;DFPDdelta;DFJCDdelta;DFJNDdelta;DFJPDdelta;DFdeltaDdelta;DFlDdelta];
ADFDdelta=prod_ext(A,DFDdelta,K);
ADFDdelta_norm=block2norm(ADFDdelta,nu);

%DFDl
DFPsiDl=Zeros2K1;
DFEDl=[0;(-para.coupling*l/para.lambda^2*Int*(C_ext+1/2*para.rho*e_ext)).*Trunc_vect];
DFCDl=[0;(Int*(-para.epsC/2*delta*C_ext)).*Trunc_vect];
DFNDl=[0;(Int*(-para.epsN/2*delta*N_ext)).*Trunc_vect];
DFPDl=[0;(Int*(-para.epsP/2*delta*P_ext)).*Trunc_vect];
DFJCDl=0;
DFJNDl=0;
DFJPDl=0;
DFdeltaDl=0;
DFlDl=0;

DFDl=[DFPsiDl;DFEDl;DFCDl;DFNDl;DFPDl;DFJCDl;DFJNDl;DFJPDl;DFdeltaDl;DFlDl];
ADFDl=prod_ext(A,DFDl,K);
ADFDl_norm=block2norm(ADFDl,nu);


coeffs=[ADFDPsi_norm,ADFDE_norm,ADFDC_norm,ADFDN_norm,ADFDP_norm,ADFDJC_norm,ADFDJN_norm,ADFDJP_norm,ADFDdelta_norm,ADFDl_norm];


function v=prod_ext(A,u,K)
K_ext=(size(u,1)-5)/5;
k=size(u,2);
ind_red=[(1:K),K_ext+(1:K),2*K_ext+(1:K),3*K_ext+(1:K),4*K_ext+(1:K),5*K_ext+(1:5)];
u_red=u(ind_red,1:k);
v_red=A*u_red;
v=u;
v(ind_red,1:k)=v_red;
end

function vect=block2norm(v,nu)
K_row=(size(v,1)-5)/5;
K_col=size(v,2);
nu_row=[1 2*nu.^(1:K_row-1)];
nu_col=[1 2*nu.^(1:K_col-1)];

vect=zeros(10,1);
if exist('intval','file') && isintval(nu)
    vect=intval(vect);
end
vect(1)=max((nu_row*abs(v(1:K_row,:)))./nu_col);
vect(2)=max((nu_row*abs(v(K_row+1:2*K_row,:)))./nu_col);
vect(3)=max((nu_row*abs(v(2*K_row+1:3*K_row,:)))./nu_col);
vect(4)=max((nu_row*abs(v(3*K_row+1:4*K_row,:)))./nu_col);
vect(5)=max((nu_row*abs(v(4*K_row+1:5*K_row,:)))./nu_col);
vect(6)=max(abs(v(5*K_row+1,:))./nu_col);
vect(7)=max(abs(v(5*K_row+2,:))./nu_col);
vect(8)=max(abs(v(5*K_row+3,:))./nu_col);
vect(9)=max(abs(v(5*K_row+4,:))./nu_col);
vect(10)=max(abs(v(5*K_row+5,:))./nu_col);
end

end