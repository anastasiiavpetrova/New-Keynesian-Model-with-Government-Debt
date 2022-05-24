function [tout, yout]= ode45modS(ypfun, y0, yfinal, compo, tol, trace)

% ODE45mod.m is ODE45.m modified to stop when a terminal 
% condition yfinal for the integration variable y(compo) is met.
% Accuracy is set with tol2

%	INPUT:
%	F     - String containing name of user-supplied problem description.
%	        Call: yprime = fun(t,y) where F = 'fun'.
%	        t      - Time (scalar).
%	        y      - Solution column-vector.
%	        yprime - Returned derivative column-vector; yprime(i) = dy(i)/dt.
%	y0    - Initial value of y.
%	yfinal- Final value of y(compo), variable in y-vector.
%  compo - Identifies element of y-vector to which the terminal condition
%          yfinal applies.
%	tol   - The desired accuracy. (Default: tol = 1.e-6).
%	trace - If nonzero, each step is printed. (Default: trace = 0).
%
%	OUTPUT:
%	T  - Returned integration time points (column-vector).
%	Y  - Returned solution, one solution column-vector per tout-value.
%  Reverted to forward looking


%	C.B. Moler, 3-25-87, 8-26-91, 9-08-92.
%	Copyright (c) 1984-94 by The MathWorks, Inc.

%New Line: Accuracy in Approximation of final value 

tol2=1.e-4;


% The Fehlberg coefficients:

alpha = [1/4  3/8  12/13  1  1/2]';
beta  = [ [    1      0      0     0      0    0]/4
          [    3      9      0     0      0    0]/32
          [ 1932  -7200   7296     0      0    0]/2197
          [ 8341 -32832  29440  -845      0    0]/4104
          [-6080  41040 -28352  9295  -5643    0]/20520 ]';
gamma = [ [902880  0  3953664  3855735  -1371249  277020]/7618050
          [ -2090  0    22528    21970    -15048  -27360]/752400 ]';
pow = 1/5;
if nargin < 4, tol = 1.e-6; end
if nargin < 5, trace = 0; end

% Initialization

t = 0;
final=0;
hmax = 10000/16;
h = hmax/8;
y = y0(:);
f = zeros(length(y),6);
chunk = 128;
tout = zeros(chunk,1);
yout = zeros(chunk,length(y));
k = 1;
tout(k) = t;
yout(k,:) = y.';

% Determine Direction of Integration

dir=sign(y(compo)-yfinal);

% Integrate until Terminal Condition is Fulfilled

while (dir==1 && y(compo) - yfinal>tol2) || (dir==-1 &&  y(compo) - yfinal<-tol2)
   
   % Compute the slopes

   temp = feval(ypfun,t,y);
   f(:,1) = temp(:);
   for j = 1:5
      temp = feval(ypfun, t+alpha(j)*h, y+h*f*beta(:,j));
      f(:,j+1) = temp(:);
   end

   % Estimate the error and the acceptable error

   delta = norm(h*f*gamma(:,2),'inf');
   tau = tol*max(norm(y,'inf'),1.0);

   % Update the solution only if the error is acceptable

   % Check Approximation of final value

   testy = y + h*f*gamma(:,1);
   if ((dir==1 && testy(compo) - yfinal<0) || (dir==-1 &&  testy(compo) - yfinal>0))
      final=1;
   else
      final=0;
   end 

   disp([t,testy(compo) - yfinal])
   if (delta <= tau && final==0) 
      t = t + h;
      y = y + h*f*gamma(:,1);
      k = k+1;
      if k > length(tout)
         tout = [tout; zeros(chunk,1)];
         yout = [yout; zeros(chunk,length(y))];
      end
      tout(k) = t;
      yout(k,:) = y.';
   end  
   
   % Update the step size

   if delta ~= 0.0

      % Bisection Method for Approximation of final value

      if final==0
      h = min(hmax, 0.8*h*(tau/delta)^pow);
      else
         h=h/2;
      end 
   else
      h=h/2;
        
   end 
   
      
   
end 

tout = tout(1:k);
yout = yout(1:k,:);

