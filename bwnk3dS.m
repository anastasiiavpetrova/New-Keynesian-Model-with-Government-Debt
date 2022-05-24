% bwnk3dS.m
% Find time paths for the New Keynesian model using the backward integration method
% Reference: Brunner, M. and H. Strulik (2002):
% Solution of Perfect Foresight Saddlepoint Problems:
% A Simple Method and Applications
% uses bwnk3ddot2.m and ode45ext.m

% ct = consumption
% bt = government bonds
% pt = inflation
% S0 = steady state
% bs=S(2) - targeted level of government spending (debt)

echo off; clear all;
global r tau teta rho eps phi g M F xi btfinal B0mn B0mx C0mn C0mx P0mn P0mx cmin cmax bmin bmax pmin pmax dir 

% *********************** Initialize Parameters ***********************

% r       		% targeted level of real interest rate
% tau           % targeted level of lump-sum tax
% teta          % adjustment costs
% rho          	% rate of time preference
% eps           % elasticity of substitution for the heterogeneous products
% phi           % inverse of the labor supply elasticity
% g             % government spending
% M             % Taylor coefficient of monetary policy rule
% F             % Taylor coefficient of fiscal policy rule
% xi            % death probability
% It's important:  xi*(xi+rho) = r*(r-rho)/(tau-g)
% theta=teta

r=0.0401; tau=0.168; teta=100; rho=0.04; eps=10; phi=0.5; g=0.15; M=-0.1; F=-0.1; xi=0.005;

% ******************* Steady State ************************************

S0 = [1,(tau-g)/r,0]';

% ******************* Initial Step ************************************

% Eigenvalues ana eigenvectors of the Jacobian matrix at steady state

[V0,D0]=eig(nk3dJacBW(S0, r, tau, teta, rho, eps, phi, g, M, F));
[~, IX0] = sort(real(diag(D0)),'descend');

% We select the MAXIMUM (positive) eigenvalue of the INVERTED Jacobian
% Thus, it's the stable (negative) eigenvalue of the original Jacobian

vU_S0 = real(V0(:, IX0(1)));


% Step Direction

uniteig = vU_S0 / norm(vU_S0);
uniteig = uniteig / sign(uniteig(2));
DD=uniteig;

% Step Length

mu=10^(-4);

%********** Initial Value and Direction of Integration ****************

for dir=0:1
    if dir==0 
        btfinal=0.9*S0(2); 
    else
        btfinal=1.1*S0(2); 
    end
    xstart = S0 + mu * uniteig * (-1)^(dir + 1);
    if dir==0 C0mx=xstart(1); end
    if dir==1 C0mn=xstart(1); end
    if dir==0 B0mx=xstart(2); end
    if dir==1 B0mn=xstart(2); end
    if dir==0 P0mx=xstart(3); end
    if dir==1 P0mn=xstart(3); end

    % ******************* Backward Integration ****************************
   
    [t,x] = ode45modS('bwnk3ddot2',xstart,btfinal,2,10^(-7),1);
    format long
    Z=[t,x];

    % ******************* Calculation of Utility ****************************
   
    TT=Z(:,1);
    Cv=x(:,1);
    Cu=log(x(:,1))-x(:,1).^(1+phi)./(1+phi);
    Exu=exp(-(xi+rho).*Z(:,1));
    a=5; % to ensure a positive utility value
    Ut=Cu.*Exu + a;
    Utf=cumtrapz(Ut);
    utility = Utf(end)
    
    %**************** Reverse time to forward-looking **************************

    t=max(t)-t;

    % Reverse Solution Vector 

    t=t(length(t):-1:1);
    x=x(length(t):-1:1,:);

    % Extract ct, bt, pt

    ct=x(:,1); bt=x(:,2); pt=x(:,3);
    ctrel=ct*100; btrel=bt; ptrel=pt*100;

    % Scale of plots

    % Scale ct

    if S0(1)-min([min(ct),C0mn])<max([max(ct),C0mx])-S0(1) 
    cmin=min([min(ct),C0mn])-max([max(ct),C0mx])+S0(1);
    cmax=max([max(ct),C0mx]);   
    end
    if S0(1)-min([min(ct),C0mn])>max([max(ct),C0mx])-S0(1) 
    cmax=max([max(ct),C0mx])+S0(1)-min([min(ct),C0mn]);
    cmin=min([min(ct),C0mn]);
    end
    
    % Scale bt

    if S0(2)-min([min(bt),B0mn])<max([max(bt),B0mx])-S0(2) 
    bmin=min([min(bt),B0mn])-max([max(bt),B0mx])+S0(2);
    bmax=max([max(bt),B0mx]);  
    end
    if S0(2)-min([min(bt),B0mn])>max([max(bt),B0mx])-S0(2) 
    bmax=max([max(bt),B0mx])+S0(2)-min([min(bt),B0mn]);
    bmin=min([min(bt),B0mn]);
    end
    
    % Scale pt

    if S0(3)-min([min(pt),P0mn])<max([max(pt),P0mx])-S0(3) 
    pmin=min([min(pt),P0mn])-max([max(pt),P0mx])+S0(3);
    pmax=max([max(pt),P0mx]);  
    end
    if S0(3)-min([min(pt),P0mn])>max([max(pt),P0mx])-S0(3) 
    pmax=max([max(pt),P0mx])+S0(3)-min([min(pt),P0mn]);
    pmin=min([min(pt),P0mn]);
    end
    
    % ***************************** Plot the Results **********************

    subplot(2,2,1)
    set(gcf,'PaperUnits', 'centimeters','PaperType','a4letter','Units','centimeters','PaperPosition', [0.5 1.5 15 7.5],'Position',[0.5 1.5 15 7.5],'Color','white');
    plot3(x(:,1),x(:,2),x(:,3),'Color', 'blue')
    hold on
    plot3(S0(1), S0(2), S0(3), '.', 'markersize', 12, 'Color', 'red')
    plot3(xstart(1), xstart(2), xstart(3), '.', 'markersize', 12, 'Color', 'black')
    set(gca,'FontName','Times','FontSize',10);
    axis([min(ct)-0.00001 max(ct)+0.00001 min(bt)-0.01 max(bt)+0.01 min(pt)-0.00001 max(pt)+0.00001])
    xlabel('ct')
    ylabel('bt')
    zlabel('pt')
    grid on;
    
    subplot(2,2,2)
    plot(t,ct,'b-')
    hold on
    plot(t,S0(1),'.','Color', 'red')
    axis([0 max(t) cmin cmax])
    xlabel('time')
    ylabel('ct')
    
    subplot(2,2,3)
    plot(t,bt,'b-')
    hold on
    plot(t,S0(2),'.','Color', 'red')
    axis([0 max(t) bmin bmax])
    xlabel('time')
    ylabel('bt')
    
    subplot(2,2,4)
    plot(t,pt,'b-')
    hold on
    plot(t,S0(3),'.','Color', 'red')
    axis([0 max(t) pmin pmax])
    xlabel('time')
    ylabel('pt')

end

