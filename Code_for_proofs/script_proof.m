clear variables
close all
clc

%"Free" parameters for which we want to compute solutions (and validate them)
pH=7;
Va_ini=0.7;%In physical units (V/NHE)
Va_end=0.0;%In physical units (V/NHE)
nbsteps_V=50;%number of intermediate values in between

%Computation of all the other (nondimensionalized) parameters of the model
para=initialize_para(pH,Va_ini);

K=300;%Number of Chebyshev modes
%Typically, we need way fewer modes to just represent the solution
%accurately, but we will need more modes for the proofs.

%Parameters for Newton's method
tol=10^-8;
it_max=20;

%Parameters used for plotting the solutions
a=-1;
b=1;
nb=500;

%% Getting a solution for the potential Va_ini
nbsteps_coupling=100;%For the numerical continuation, in case we need to compute a solution from scratch
disp=1;%Put to 1 to plot the computed solutions
X=get_sol(Va_ini,para,K,nbsteps_coupling,tol,it_max,disp);
para.coupling=1;

%% Continuation in V
tab_V=linspace(Va_ini/para.resc_V,Va_end/para.resc_V,nbsteps_V);
tab_tJ=0*tab_V;
tab_X=zeros(5*K+5,length(tab_V));
counter=1;
fprintf('\nContinuation in V, from Va=%.4f to Va=%.4f\n',tab_V(1)*para.resc_V,tab_V(end)*para.resc_V)
for V=tab_V
    para.V=V;
    it=0;
    dist=tol+1;
    fprintf('\nCurrent value of Va: %.4f',V*para.resc_V)
    fprintf('\nNewton method, distance between two iterates:\n')
    while it<it_max && dist>tol && dist<Inf && not(isnan(dist))
        Xold=X;
        X=X-DF(X,para)\F(X,para);
        dist=norm(X-Xold,1)
        it=it+1;
    end
    JC=X(5*K+1);
    JN=X(5*K+2);
    JP=X(5*K+3);
    l=X(5*K+5);
    tildeJ=(-3*(JC/(4*para.epsC)+JP)+JN/para.epsN)/(2*l);
    tab_tJ(counter)=tildeJ;
    tab_X(:,counter)=X;
    counter=counter+1;
end

if disp %plot of the last solution 
    Psi=X(1:K);
    C=X(2*K+1:3*K);
    N=X(3*K+1:4*K);
    P=X(4*K+1:5*K);
    delta=X(5*K+4);
    l=X(5*K+5);
    figure
    plotCheb(para.resc_V*Psi,a,b,nb,'k')
    legend('$\Psi$', 'interpreter', 'latex', 'location', 'northwest')
    title(['$V_a$ = ',num2str(para.resc_V*V)], 'interpreter', 'latex')
    set(gca,'FontSize',15)
    figure
    plotCheb(P,a,b,nb,'r')
    hold on
    plotCheb(N,a,b,nb,'b')
    plotCheb(C,a,b,nb,'g')
    legend('P','N','C', 'location', 'northeast')
    title(['$V_a$ = ',num2str(para.resc_V*V)], 'interpreter', 'latex')
    fprintf('\n l = %f(nm), delta = %f(\x03bcm/year)\n',l*para.resc_l,delta*para.resc_delta)
    set(gca,'FontSize',15)
    drawnow
end

tab_V=tab_V*para.resc_V;
tab_tJ=tab_tJ*para.resc_tildeJ;
tab_delta=tab_X(5*K+4,:)*para.resc_delta;
tab_l=tab_X(5*K+5,:)*para.resc_l;

figure
plot(tab_V,tab_l,'k','Linewidth',2)
xlabel({'$V_a \ (V/NHE)$'},'Interpreter','Latex')
ylabel({'$l \ (nm)$'},'Interpreter','Latex')
set(gca,'FontSize',15)
axis([min(tab_V) max(tab_V) min(tab_l)-0.1 max(tab_l)+0.1])

figure
plot(tab_V,tab_delta,'k','Linewidth',2)
xlabel({'$V_a \ (V/NHE)$'},'Interpreter','Latex')
ylabel({'$\delta \ (\mu m/year)$'},'Interpreter','Latex')
set(gca,'FontSize',15)
axis([min(tab_V) max(tab_V) min(tab_delta)-0.001 max(tab_delta)+0.001])

figure
plot(tab_V,tab_tJ,'k','Linewidth',2)
xlabel({'$V_a \ (V/NHE)$'},'Interpreter','Latex')
ylabel({'$\tilde{J} \  (\mu A/cm^2)$'},'Interpreter','Latex')
set(gca,'FontSize',15)
axis([min(tab_V) max(tab_V) min(tab_tJ)-0.001 max(tab_tJ)+0.001])
drawnow


%% Proofs 
% Here we validate a posteriori the approximate solutions obtained above.
% The first "validation" is non rigorous as it is done with usual floatting
% point arithmetic. The second one is truly a rigorous validation, as it
% uses interval arithmetic (thanks to Intlab). If you do not have Intlab,
% this last part will be skipped.

rstar=10^-5;% size of the neighborhood of the numerical solution in which we try to prove the existence of a true solution
nu=1.1;% weight for the l^1 norm on Chebyshev coefficients
criterion_opt_eta='Largest_interval';% Select how we compute the weigths in eta (see optimize_eta.m)

tab_proofs=zeros(14,length(tab_V));
counter=1;
eta0=[0.1, 0.05, 1, 0.3, 1, 0.01, 1, 0.01, 0.01, 0.5];
%eta0=ones(1,10);
fprintf('\n\nVALIDATION\n\n')
for V=tab_V
    Va=V
    para.V=V/para.resc_V;
    X=tab_X(:,counter);
    fprintf('\n"Proof" without interval arithmetic\n')
    [result,rmin,rmax,eta,A]=proof(X,para,nu,rstar,criterion_opt_eta,eta0);% Non rigorous validation (and computation of the weight eta)
    tab_proofs([1,2,5:14],counter)=[result,rmin,eta];
    if result
        eta0=eta;
        if exist('intval','file')
            ipara=para2intval(para);% turns the parameters into intvals
            iX=intval(X);
            inu=intval(nu);
            ieta=intval(eta);
            iA=intval(A);
            fprintf('\nTrue proof with interval arithmetic\n')
            [iresult,irmin,irmax,~,~]=proof(iX,ipara,inu,rstar,'no',[],ieta,iA);% Rigorous validation (using the weight eta computed in the non-rigorous validation)
            tab_proofs(3:4,counter)=[iresult;irmin];
        end
    end
    counter=counter+1;
end

fprintf('\n%d out of %d solutions could be "validated" without interval arithmetic\n',sum(tab_proofs(1,:)),size(tab_proofs,2))
BC0=[1 2*(-1).^(1:K-1)];
if exist('intval','file')
    fprintf('\n%d out of %d solutions could be rigorously validated using interval arithmetic\n',sum(tab_proofs(3,:)),size(tab_proofs,2))
    %Postproccessing the results
    mask=logical(tab_proofs(3,:));
    err_delta=sup(max(tab_proofs(4,mask)./tab_proofs(13,mask))*para.resc_delta);
    err_l=sup(max(tab_proofs(4,mask)./tab_proofs(14,mask))*para.resc_l);
    tab_itJ=intval(zeros(1,length(tab_V)));
    for counter=1:length(tab_V)
        iX=intval(tab_X(:,counter));
        ieta=intval(tab_proofs(5:14,counter));
        iPsi=iX(1:K);
        errPsi=midrad(0,sup(tab_proofs(4,counter)/ieta(1)));
        iN=iX(3*K+1:4*K);
        errN=midrad(0,sup(tab_proofs(4,counter)/ieta(4)));
        iJC=iX(5*K+1);
        iJP=iX(5*K+3);
        il=iX(5*K+5);
        itildeJ=-3*(iJC/(4*ipara.epsC)+iJP)/(2*il)-ipara.rN0(BC0*[iN iPsi]+[errN errPsi])/ipara.epsN;
        tab_itJ(counter)=itildeJ;
    end
    tab_itJ=tab_itJ*para.resc_tildeJ;
    err_tJ=max(rad(tab_itJ(mask)));
    fprintf('\nAmong the rigorously validated solutions:\nthe error on each value of delta is bounded by %e (\x03bcm/year),\nthe error on each value of l is bounded by %e (nm),\nthe error on each value of J is bounded by %e (\x03bc A/cm^2).\n',err_delta,err_l,err_tJ)

    figure
    plot(tab_V(mask),tab_l(mask),'*k','Linewidth',2)
    xlabel({'$V_a \ (V/NHE)$'},'Interpreter','Latex')
    ylabel({'$l \ (nm)$'},'Interpreter','Latex')
    set(gca,'FontSize',15)
    axis([min(tab_V(mask)) max(tab_V(mask)) min(tab_l(mask))-0.1 max(tab_l(mask))+0.1])

    figure
    plot(tab_V(mask),tab_delta(mask),'*k','Linewidth',2)
    xlabel({'$V_a \ (V/NHE)$'},'Interpreter','Latex')
    ylabel({'$\delta \ (\mu m/year)$'},'Interpreter','Latex')
    set(gca,'FontSize',15)
    axis([min(tab_V(mask)) max(tab_V(mask)) min(tab_delta(mask))-0.001 max(tab_delta(mask))+0.001])

    figure
    plot(tab_V(mask),mid(tab_itJ(mask)),'*k','Linewidth',2)
    xlabel({'$V_a \ (V/NHE)$'},'Interpreter','Latex')
    ylabel({'$\tilde{J} \  (\mu A/cm^2)$'},'Interpreter','Latex')
    set(gca,'FontSize',15)
    axis([min(tab_V(mask)) max(tab_V(mask)) min(mid(tab_itJ(mask)))-0.001 max(mid(tab_tJ(mask)))+0.001])
    drawnow
else
    fprintf('\nThe validation was done with floating point arithmetic and is not guaranteed. If you have Intlab, make sure it is in your working directory in order to get certified results.\n')
end
