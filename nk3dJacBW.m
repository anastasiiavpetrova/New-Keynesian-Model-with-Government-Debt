% Jacobian matrix of the system bwnkddot2.m

function jac = nk3dJacBW(x, r, tau, teta, rho, eps, phi, g, M, F)

       %Jacobian of the system - backward: 

       jac = [-(r-rho+M*x(3)), r*(r-rho)/(tau-g), -M*x(1); 
                0, -(r-F+M*x(3)), -F*x(2); 
                (r*(r-rho)/(tau-g)*x(2)*x(3)/x(1)^2)+(eps-1)*(1+phi)*x(1)^phi/teta, -r*(r-rho)/(tau-g)*x(3)/x(1), -r*(r-rho)/(tau-g)*x(2)/x(1)-rho];

end