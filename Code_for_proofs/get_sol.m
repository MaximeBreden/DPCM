function X=get_sol(Va,para,K,nbsteps,tol,it_max,disp)

%Load or compute an approximate solution X for the given Va (applied
%potential in physical units) and the value of pH given in para.
%The output will have K Chebyshev modes. 
%If we need to compute the solution, we first find a solution without
%coupling (i.e. with 0 in the r.h.s. of the poisson equation), and then use
%a (basic) numerical continuation to get up to the model with coupling,
%using nbsteps number of steps.
%Put disp to 1 if you want to plot the computed solutions and to 0
%otherwise.


%Parameters used for plotting the solutions
a=-1;
b=1;
nb=500;

str=['sol_Va_',num2str(Va),'_pH_',num2str(para.pH),'.mat'];
if exist(str,'file')
    %% Loading a precomputed solution
    fprintf('\nLoading a precomputed solution\n')
    load(str,'X')
    para.coupling=1;
    X=change_size(X,K);
    
    if disp %plot
        Psi=X(1:K);
        C=X(2*K+1:3*K);
        N=X(3*K+1:4*K);
        P=X(4*K+1:5*K);
        delta=X(5*K+4);
        l=X(5*K+5);
        figure
        plotCheb(para.resc_V*Psi,a,b,nb,'k')
        legend('$\Psi$', 'interpreter', 'latex', 'location', 'northwest')
        title(['$V_a$ = ',num2str(Va)], 'interpreter', 'latex')
        set(gca,'FontSize',15)
        figure
        plotCheb(P,a,b,nb,'r')
        hold on
        plotCheb(N,a,b,nb,'b')
        plotCheb(C,a,b,nb,'g')
        legend('P','N','C', 'location', 'northeast')
        title(['$V_a$ = ',num2str(Va)], 'interpreter', 'latex')
        fprintf('\n l = %f(nm), delta = %f(\x03bcm/year)\n',l*para.resc_l,delta*para.resc_delta)
        set(gca,'FontSize',15)
        drawnow
    end
    
else
    fprintf('\nNo precomputed solution for this choice of V, lets compute one...\n')
    %% First aproximate solution without coupling
    lmin=10^-1;
    lmax=10^8;
    tol_dicho=10^-2;
    it_max_dicho=10^5;
    fprintf('\nComputation of an initial solution (potentiostatic case, without coupling)\n')
    [l,JC,JN,JP,delta,C_func,N_func,P_func,Psi_0,Psi_1]=find_sol_uncoupled(para,lmin,lmax,tol_dicho,it_max_dicho);

    fprintf('\nConversion of this solution into Chebyshev series\n')
    C=computeCheb(C_func,K);
    N=computeCheb(N_func,K);
    P=computeCheb(P_func,K);

    %% Refinement of the solution (still without coupling) using Newton's method

    X=zeros(5*K+5,1);
    X(1)=(Psi_0+Psi_1)/2;
    X(2)=(Psi_1-Psi_0)/4;
    X(K+1)=(Psi_1-Psi_0)/2;
    X(2*K+(1:K))=C;
    X(3*K+(1:K))=N;
    X(4*K+(1:K))=P;
    X(5*K+(1:5))=[JC;JN;JP;delta;l];

    fprintf('\nRefinement using Newton method, distance between two iterates:\n')
    it=0;
    dist=tol+1;
    while it<it_max && dist>tol && dist<Inf
        Xold=X;
        X=X-DF(X,para)\F(X,para);
        dist=norm(X-Xold,1)
        it=it+1;
    end

    if disp %plot
        Psi=X(1:K);
        C=X(2*K+1:3*K);
        N=X(3*K+1:4*K);
        P=X(4*K+1:5*K);
        delta=X(5*K+4);
        l=X(5*K+5);
        figure
        plotCheb(para.resc_V*Psi,a,b,nb,'k')
        legend('$\Psi$', 'interpreter', 'latex', 'location', 'northwest')
        title(['$V_a$ = ',num2str(Va),' (without coupling)'], 'interpreter', 'latex')
        set(gca,'FontSize',15)
        figure
        plotCheb(P,a,b,nb,'r')
        hold on
        plotCheb(N,a,b,nb,'b')
        plotCheb(C,a,b,nb,'g')
        legend('P','N','C', 'location', 'northeast')
        title(['$V_a$ = ',num2str(Va),' (without coupling)'], 'interpreter', 'latex')
        fprintf('\n l = %f(nm), delta = %f(\x03bcm/year)\n',l*para.resc_l,delta*para.resc_delta)
        set(gca,'FontSize',15)
        drawnow
    end

    %% Continuation to add the coupling
    fprintf('\nNumerical continuation to add the coupling, from s=0 t0 s=1\n')
    tab_s=linspace(0,1,nbsteps);
    for s=tab_s
        para.coupling=s;
        fprintf('\nCurrent value of s: %.4f',s)
        fprintf('\nNewton method, distance between two iterates:\n')
        it=0;
        dist=tol+1;
        while it<it_max && dist>tol && dist<Inf && not(isnan(dist))
            Xold=X;
            X=X-DF(X,para)\F(X,para);
            dist=norm(X-Xold,1)
            it=it+1;
        end
        if dist==Inf || isnan(dist)
            error('Newton did not converge, you should maybe try a larger number of steps (nbsteps)')
        end
    end

    if disp %plot
        Psi=X(1:K);
        C=X(2*K+1:3*K);
        N=X(3*K+1:4*K);
        P=X(4*K+1:5*K);
        delta=X(5*K+4);
        l=X(5*K+5);
        figure
        plotCheb(para.resc_V*Psi,a,b,nb,'k')
        legend('$\Psi$', 'interpreter', 'latex', 'location', 'northwest')
        title(['$V_a$ = ',num2str(Va),' (full coupling)'], 'interpreter', 'latex')
        set(gca,'FontSize',15)
        figure
        plotCheb(P,a,b,nb,'r')
        hold on
        plotCheb(N,a,b,nb,'b')
        plotCheb(C,a,b,nb,'g')
        legend('P','N','C', 'location', 'northeast')
        title(['$V_a$ = ',num2str(Va),' (full coupling)'], 'interpreter', 'latex')
        fprintf('\n l = %f(nm), delta = %f(\x03bcm/year)\n',l*para.resc_l,delta*para.resc_delta)
        set(gca,'FontSize',15)
        drawnow
    end
    
    save(str,'X');
    
    fprintf('\nSolution stored for further usage\n')

end