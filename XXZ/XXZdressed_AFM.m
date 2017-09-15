clear;
close all;
% clc;

verbose=1;

tol = 1e-10;
dig=ceil(-log10(tol));
frmt=['%2.',int2str(dig),'e'];
opts = optimset('TolX',tol,'TolFun',tol);

Delta = -3;
h = 1.2;

phi=acosh(-Delta);
Q=pi/phi;

% Takahashi defines elliptic integrals in terms of u=sqrt(m) (or u^2=m)
% i.e. K(u) = int_0^pi/2 dx (1-u^2 sin^2 x)^(-1/2)
% matlab uses m=u^2
m0 = findm0(Q,tol);
K0 = ellipke(m0);

%% critical critical hc, above which m>0
% Takahashi's eq. (4.35) is more than cryptic... It is very ambiguous as to what u' should be. In any case I
% think the 1/phi is wrong, as it cancels with the phi in the exact expression for rho0 in (4.26)
% Also, he defines H ~ -2hm, so use h-> 2h (i.e. add a multiplicative factor of 2)
hc = 2*sinh(phi)*K0*sqrt(1-m0)/pi;
disp(['hc=',num2str(hc)]);
assert(abs(h)<hc,'|h|<hc');

Ker = @(x,n)((1/(2*pi))*phi*sinh(n*phi)./(cosh(n*phi) - cos(phi*x)));
Kfun = @(x,y)(-Ker(x-y,2)-Ker(x+y,2));
rho0fun = @(x) Ker(x,1);
e0fun = @(x)(h - 2*pi*sinh(phi)/phi*Ker(x,1));

[E0,rho,x0,Edr0,m]=fXXZGS_fixedh(Delta,0,tol,verbose);

sigstr = Fie(1,0,x0,1,Kfun,@(x)(h*ones(size(x))),tol,10*tol);
sig = @(x)(ntrpFie(sigstr,x));

Edr = @(x)(Edr0(x) + sig(x));

xx=linspace(-x0,x0,200);

plot(xx,-Edr(xx),'-r',xx,-Edr0(xx)-h/2,'.b');

