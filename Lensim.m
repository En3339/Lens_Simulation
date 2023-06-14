% The code simulates a Gaussian Beam pass through a lens. the input beam waist is at the center of a focal plane in front of the lens.
% And it calculate's the complex amplitude at the other focal plane.
% It involves both Numerical method and Analytical method. and both are proved to agree with energy conservation law.
% Auther: En Lu
% Date 14/06/2023


% Parameters
lambda = 960e-9;      % Wavelength
k = 2*pi/lambda;      % Wave number
W_0 = 7e-6;           % Waist radius
f=0.03; b = lambda*f;

Nx = 1024;
dx = 0.2e-6;dy=dx;

x = (0:(Nx-1))*dx;x = fftshift(x);x(1:(find(x==0)-1))=x(1:(find(x==0)-1))-x(find(x==0)-1)-x(find(x==0)+1);x=ifftshift(x);
y = x;

fx = (0:(Nx-1))/Nx/diff(x(1:2));
fx = fftshift(fx);fx(1:(find(fx==0)-1))=fx(1:(find(fx==0)-1))-fx(find(fx==0)-1)-fx(find(fx==0)+1);fx=ifftshift(fx);
fy = fx;
u = fx*b;
v = fy*b;
dfx = diff(fx(1:2));

du = diff(u(1:2));dv=du;

% Evaluate complex amplitude at waist
% Distribute at spatial matrix of [X,Y]*b
[X,Y] = meshgrid(x,y); 
U_in = exp(-(((X).^2 + (Y).^2)/W_0^2));   % Complex amplitude at waist

%FFT
H = 1/(1i*lambda*f);
Uf = H*fft2(U_in)*dx*dy;

%Energy before the lens and after, fit energy conservation:
sum(abs(U_in(:)).^2)*dx^2
sum(abs(Uf(:)).^2)*du^2

%Analytic expression for the integral here, from the
%mathematica notebook

[U,V] = meshgrid(u,v); 
Uan = pi*exp(-pi.^2.*(U.^2+V.^2)./(sqrt(1/W_0^4)*b^2))./(b*sqrt(1/W_0^4));

sum(abs(Uan(:)).^2)*du^2

%Analytical intergration result
Interg=pi*W_0^2/2
