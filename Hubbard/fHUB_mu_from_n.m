function [mu,Emu,rho,k0,E0] = fHUB_mu_from_n(U,n,tol,rho,k0,verbose)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3 || isempty(tol),tol=1e-10;end
if nargin < 6 || isempty(verbose),verbose=true;end

if nargin < 5 || isempty(rho) || isempty(k0)
    [E0,k0,rho] = fHUBGS_fixedn(U,n,tol,verbose);
end

Rfun = @(x)(1/pi * integral(@(y)(cos(y.*x)./(1 + exp(0.5*U*y))),0,Inf,'AbsTol',tol,'RelTol',tol,'ArrayValued',true));
Kfun = @(x,y) cos(x).*(Rfun(sin(x) - sin(y)) + Rfun(sin(x) + sin(y))); % symmetric Kernel function, for less numerical effort
Gfun = @(x) rho(k0).*Kfun(x,k0);

sigstr = Fie(1,0,k0,1,Kfun,Gfun,tol,10*tol);
sig = @(x) ntrpFie(sigstr,abs(x));

dEdk = -2*(cos(k0)*rho(k0) + integral(@(x)(cos(x).*sig(x)),0,k0,'AbsTol',tol,'RelTol',tol));
dndk = rho(k0) + integral(sig,0,k0,'AbsTol',tol,'RelTol',tol);

mu = dEdk/dndk - U/2;

Emu = -mu*n;

end

