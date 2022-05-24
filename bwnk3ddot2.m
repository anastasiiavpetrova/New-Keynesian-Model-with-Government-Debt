% Three-dimensional system ofdifferential equations describing the economy

function bwnk3ddot2 = bwnk3ddot2(t,x)

global r tau teta rho eps phi g M F 
ct=x(1);bt=x(2);pt=x(3);

ctdot = (r-rho)*ct+M*pt*ct-r*(r-rho)/(tau-g)*bt;
btdot = (r-F)*bt+M*pt*bt+g-tau+F*(tau-g)/r;
ptdot = (rho+r*(r-rho)/(tau-g)*bt/ct)*pt-((eps-1)/teta)*(ct^(1+phi)-1);

%********** Backward Looking ******************************************

bwnk3ddot2 =-[ctdot;btdot;ptdot];