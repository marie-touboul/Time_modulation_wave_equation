%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%            Time-modulated 1D wave equation problem                                %%%
%%%   BEFORE USING THIS FILE: compute the effective coefficients with Mathematica     %%%
%%%   OUTCOMES:                                                                       %%%
%%%            (i) Maps of the dispersion function                                    %%%
%%%            (ii) Dispersion diagrams                                               %%%
%%%            (iii) Compute error plots and measures                                 %%%
%%%            (iv) Compute measure of non-reciprocity                                %%%
%%%                                                                                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all ; 
clear all ; 
colorplot = inferno;
options=optimset('TolFun', 1.0e-16,'TolX',1.0e-16);             %for the root finding


%% Choose what to plot
map = 'no' ;                                                    %Plot the colormap for the dispersion function
root_finding ='yes' ;                                           %Do root finding along the first branch
modul = "one" ;                                                 %'one' or 'both': Modulating one or both parameters
config = 'beta1' ;                                              %'alpha1' or 'beta1' or 'both': Configuration alpha=1 or beta=1 or both modulated, respectively

%% Physical Parameters to choose 
rho1 = 1000;                    %rho phase A
rho2 = 6000;                    %rho phase B
mu1 = 2.6643e+09;               %mu phase A
mu2 =  2.6643e+09;              %mu phase B
cm = 1800 ;                     % modulation velocity. Has to be outside (c1,c2)
h = 20 ;                        %periodicity
phi = 0.25 ;                    %proportion phase A

%% Effective coefficients computed previously with Mathematica for a given set of physical parameters
cstar = 1392;       % a characteristic speed: same as the one chosen in Mathematica to compute the effective coefficients

%coeff for the case alpha = 1 (homogenization order 2)
C_E0 = 44210500000000/22828871 ;
C_L = 27.3404 ;
C_Nbeta = 1467580901850000000000/42578027471646641 ;


%coeff for the case beta = 1 (homogenization order 2)
rho0_b1 = 1.0667 ;
C_M =  -91.5303  ;
C_Nalpha = -0.0678542 ;


%coeff for the case both parameters varying ( homogenization order 0)
E0 = 0.928709 ;
rho0 = 0.249793 ;
S0 = 0.107279 ;

%% Compute physical quantities and characteristic values
L1 = 0.25*h ;                           %length phase A
L2 = 0.75*h ;                           %length phase B
c1 = sqrt(mu1/rho1);                    %wavespeed phase A
c2 = sqrt(mu2/rho2);                    %wavespeed phase B
Z1 = sqrt(rho1*mu1) ;                   %sqrt impedance phase A (appears in the dispersion function)
Z2 = sqrt(rho2*mu2);                    %sqrt impedance phase A (appears in the dispersion function)

rhobar = phi*rho1+(1-phi)*rho2 ;        %characteristic rho
mubar = 1/(phi/mu1+(1-phi)/mu2);        %characteristic mu 
c_eta = sqrt(mubar/rhobar);             %characteriscit speed. Just to non-dimensionalize the frequency axis


k_lim = 0.69;                           %(-k_lim,k_lim) = interval in wavenumber where we compute dispersion diagrams and errors


%% DECLARE DISPERSION FUNCTIONS

% Homogenization order 2
Dalpha1 = @(k,omega) omega^2-C_E0*k.^2+C_L*k.^2*omega^2-C_Nbeta*omega*k.^3  ;                   %alpha = 1 
Dbeta1 = @(k,omega) -rho0_b1*omega^2+k.^2*cstar^2-C_M*k.^2*omega^2+C_Nalpha*omega^3.*k;         %beta = 1

% Homogenization order 0 
Dalpha_order_zero = @(k,omega) -omega^2+C_E0*k.^2  ;                                            %alpha=1
Dbeta_order_zero = @(k,omega) -rho0_b1*omega^2+k.^2*cstar^2;                                    %beta = 1
Dboth = @(k,omega) -rho0*omega^2+k.^2*E0*cstar^2-2*S0*k.*omega*cstar;                           %Both modulated
Dboth_zero = @(k,omega) -rho0*omega^2+k.^2*E0*cstar^2;                                          %Both modulated - limit S0 = 0 

% Exact Dispersion function in the general case - Obtained by Bloch-Floquet analysis
D_micro = @(k,omega) cos(k*h-cm*(omega-cm*k)*(L1/(c1^2-cm^2)+L2/(c2^2-cm^2)))-cos((omega-cm*k)/(c1^2-cm^2)*c1*L1)*cos((omega-cm*k)/(c2^2-cm^2)*c2*L2)+1/2*(Z1/Z2+Z2/Z1)*sin((omega-cm*k)/(c1^2-cm^2)*c1*L1)*sin((omega-cm*k)/(c2^2-cm^2)*c2*L2);


%% Wavenumber and frequency range
%For the maps
k_map = [(-pi:0.02:-0.008)/h flip(-(-pi:0.02:-0.008)/h)];                                       %k interval
Omega_plot = 0:0.005:1.5*c_eta/h ;                                                              %Omega interval 
%For the root finding 
k_vect = [(-k_lim:0.005:-0.0008)/h flip(-(-k_lim:0.005:-0.0008)/h)];                            %k interval
%Initializations of Omega for the root finding (guess the location of the
%from the maps)
Omega_init_homog2 = 0.9*c_eta/h ;                                                               
Omega_init_homog0 = 0.9*c_eta/h ;
Omega_init_micro = 0.5*c_eta/h ;

%% Plots for one parameter modulated 
if modul == "one"
    if config == "alpha1"
     D_homog0 = Dalpha_order_zero ; 
     D_homog2 = Dalpha1 ; 
    end 
    if config == "beta1" 
        D_homog0 = Dbeta_order_zero ;
        D_homog2 = Dbeta1 ;
    end
    if map == "yes"                                                                         
        %Initialization for zero finding for the homogenized dispersion
        %relations
        %Order 0
        omega_zero = Omega_init_homog0 ; 
        vect_Omega_order0_map = zeros(1,numel(k_map));    
        %Order 2
        omega = Omega_init_homog2 ;
        vect_Omega_map= zeros(1,numel(k_map));                
        %Root finding for homogenized models at order 0 and 2      
        for i = 1:numel(k_map)
            omega = fsolve(@(z) D_homog2(k_map(i),z),omega_zero,options) ;
            vect_Omega_map(1,i) = omega ; 
            omega_zero = fsolve(@(z) D_homog0(k_map(i),z),omega_zero,options) ;
            vect_Omega_order0_map(1,i) = omega_zero ;   
        end
        %Map for the exact dispersion function
        DD_micro = zeros(numel(Omega_plot),numel(k_map));
        for ind = 1:numel(k_map)
            for ind2 = 1:numel(Omega_plot)
                DD_micro(ind2,ind) = D_micro(k_map(ind),Omega_plot(ind2));
            end
        end
        %Plots
        figure ; 
        pcolor(k_map*h,Omega_plot*h/c_eta,log10(abs(DD_micro)))
        colormap(colorplot())
        shading interp ;
        hold on ; 
        plot(k_map(abs(k_map*h)<1.5)*h,vect_Omega_map(abs(k_map*h)<1.5)*h/c_eta,':','Linewidth',2) ; 
        hold on ; 
        plot(k_map*h,vect_Omega_order0_map*h/c_eta,'--k','Linewidth',2) ; 
        colorbar ;
        %ylim([0 1.5])
         ylabel('$\eta(\omega)$','Interpreter','Latex','Fontsize',20) ;
        xlabel('$kh$','Interpreter','Latex','Fontsize',20) ;
        x0=10;
        y0=10;
        width=650;
        height=500;
        set(gcf,'position',[x0,y0,width,height])
        title('$\mathrm{log}_{10}(\mathrm{Disp}(\omega,k))$','Interpreter','Latex','Fontsize',20) ;
        legend('Exact disp. func.','Second-order homog.','Leading-order homog.','Fontsize',11,'Location','Southeast')
    end

    if root_finding == "yes"
        %Initializations
        vect_Omega= zeros(1,numel(k_vect));                %homogenized order 2
        omega = Omega_init_homog2 ;
        vect_Omega_order0 = zeros(1,numel(k_vect));        %homogenized order 0
        omega_zero = Omega_init_homog0 ; 
        vect_Omega_micro = zeros(1,numel(k_vect));         %exact one
        omega_micro = Omega_init_micro ;
       
        %solve for omega from k_i to k_{i+1}
        for i = 1:numel(k_vect)
            omega = fsolve(@(z) D_homog2(k_vect(i),z),omega,options) ;
            vect_Omega(1,i) = omega ; 
            omega_zero = fsolve(@(z) D_homog0(k_vect(i),z),omega_zero,options) ;
            vect_Omega_order0(1,i) = omega_zero ;   
            omega_micro = fsolve(@(z) D_micro(k_vect(i),z),omega_micro,options) ;
            vect_Omega_micro(1,i) = omega_micro ;   
        end
    
        %Plot the zero finding results near the origin 
        figure ; 
        plot(k_vect*h,vect_Omega_micro*h/c_eta,'r','Linewidth',2) ; 
        hold on ; 
        plot(k_vect*h,vect_Omega_order0*h/c_eta,'--k','Linewidth',2) ; 
        hold on ; 
        plot(k_vect*h,vect_Omega*h/c_eta,':','Linewidth',3,'color',[0 0.4470 0.7410]) ; 
        hold on ;         
        ylabel('$\eta(\omega)$','Interpreter','Latex','Fontsize',20) ;
        xlabel('$kh$','Interpreter','Latex','Fontsize',20) ;
        if config == "alpha1"
            title('$\alpha=1$','Interpreter','Latex','Fontsize',20)
        end
        if config == "beta1"
            title('$\beta=1$','Interpreter','Latex','Fontsize',20)
        end  
        legend('Bloch-Floquet solution','Leading- or first-order homogenized model','Second-order homogenized model','Fontsize',11,'Location','North')
        xlim([-k_lim k_lim]); 
        %ylim([0 0.6])
    
        % Plot errors in log-log scale 
        % negative wavenumbers
        ind_neg = k_vect < 0 ; 
        figure ; 
        loglog(abs(k_vect(ind_neg)*h),abs((vect_Omega(ind_neg).^2-vect_Omega_micro(ind_neg).^2)),'-','Linewidth',2) ;
        hold on ; 
        loglog(abs(k_vect(ind_neg)*h),abs((vect_Omega_order0(ind_neg).^2-vect_Omega_micro(ind_neg).^2)),'-','Linewidth',2) ;
        hold on ; 
        loglog(abs(k_vect(ind_neg)*h),abs(10*k_vect(ind_neg)*h).^4,'--','Linewidth',2) ;
        hold on ; 
        loglog(abs(k_vect(ind_pos)*h),abs(k_vect(ind_pos)*h).^6,'--','Linewidth',2) ;
        xlabel('$-kh$','Interpreter','Latex','Fontsize',20)
        ylabel('absolute errors','Interpreter','Latex','Fontsize',20)
        title('negative wavenumbers','Interpreter','Latex','Fontsize',20)
        legend('Second-order homogenized model','Leading- or first-order homogenized model','$\mathcal{O}(k^4)$','$\mathcal{O}(k^6)$','Interpreter','Latex','Fontsize',11)
        set(gca,'TickLabelInterpreter','latex','fontsize',15);   
        % positive wavenumbers
        ind_pos = k_vect > 0 ;
        figure ; 
        loglog(abs(k_vect(ind_pos)*h),abs((vect_Omega(ind_pos).^2-vect_Omega_micro(ind_pos).^2)),'-','Linewidth',2) ;
        hold on ; 
        loglog(abs(k_vect(ind_pos)*h),abs((vect_Omega_order0(ind_pos).^2-vect_Omega_micro(ind_pos).^2)),'-','Linewidth',2) ;
        hold on ; 
        loglog(abs(k_vect(ind_pos)*h),abs(10*k_vect(ind_pos)*h).^4,'--','Linewidth',2) ;
        hold on ; 
        loglog(abs(k_vect(ind_pos)*h),abs(0.5*k_vect(ind_pos)*h).^6,'--','Linewidth',2) ;
        xlabel('$kh$','Interpreter','Latex','Fontsize',20)
        ylabel('absolute errors','Interpreter','Latex','Fontsize',20)
        title('positive wavenumbers','Interpreter','Latex','Fontsize',20)
        legend('Second-order homogenized model','Leading- or first-order homogenized model','$\mathcal{O}(k^4)$','$\mathcal{O}(k^6)$','Interpreter','Latex','Fontsize',11)
        set(gca,'TickLabelInterpreter','latex','fontsize',15);

    
        %Compute "non-reciprocity measures"
        vect_Omega_odd = (vect_Omega-flip(vect_Omega))/2;
        NR_alpha1 = sqrt(trapz(k_vect,vect_Omega_odd.^2)/trapz(k_vect,vect_Omega.^2))*100
    
        %Compute errors between homogenized model and exact relation
        NR_error_homog2 = sqrt(trapz(k_vect,(vect_Omega_micro-vect_Omega).^2)/trapz(k_vect,(vect_Omega_micro).^2))*100
        NR_error_homog0 = sqrt(trapz(k_vect,(vect_Omega_micro-vect_Omega_order0).^2)/trapz(k_vect,(vect_Omega_micro).^2))*100
    end
end

%% Plots for both parameters modulated 
%Plot zeroth order homogenized model
if config == "both"

    %Initialization
    vect_Omega_both = zeros(1,numel(k_vect));               %homogenized modulated both parameters
    omega = 7.5*c_eta/h ;
    vect_Omega_both_zero = zeros(1,numel(k_vect));          %homogenized S0=0, no Willis coupling
    omega_zero = 6*c_eta/h;
    vect_Omega_micro = zeros(1,numel(k_vect));              % exact one
    omega_micro = 150 ;

    %Root finding 
    for i = 1:numel(k_vect)
        omega = fsolve(@(z) Dboth(k_vect(i),z),omega) ;
        vect_Omega_both(1,i) = omega ; 
        omega_zero = fsolve(@(z) Dboth_zero(k_vect(i),z),omega) ;
        vect_Omega_both_zero(1,i) = omega_zero ; 
        omega_micro = fsolve(@(z) D_micro(k_vect(i),z),omega_micro) ;
        vect_Omega_micro(1,i) = omega_micro ;             
    end
    
    %Plots 
    figure;
    plot(k_vect*h,vect_Omega_both*h/c_eta,'Linewidth',3) ; 
    hold on ; 
    plot(k_vect*h,vect_Omega_both_zero*h/c_eta,'--k','Linewidth',3) ; 
%     hold on ; 
%     plot(k_vect,vect_Omega_micro*h/c_eta,'-o','Linewidth',2) ; 
    ylabel('$\eta(\omega)$','Interpreter','Latex','Fontsize',26) ;
    xlabel('$kh$','Interpreter','Latex','Fontsize',26) ;
    legend('$W_0 \neq 0$','$W_0=0$','Interpreter','Latex','Fontsize',22,'Location','SouthEast')
    xlim([-pi pi]); 
    set(gca,'TickLabelInterpreter','latex','fontsize',18);
    %ylim([0 0.6])

    %Compute "non-reciprocity measures"
    Nk = numel(k_vect);
    vect_Omega_both_odd = (vect_Omega_both-flip(vect_Omega_both))/2;
    NR_both = sqrt(trapz(k_vect,vect_Omega_both_odd.^2)/trapz(k_vect,vect_Omega_both.^2))*100
end


