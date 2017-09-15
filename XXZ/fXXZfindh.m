function [h,Eh] = fXXZfindh(Delta,m,rho,B,tol,verbose)
% calculates the magnetic field that corresponds to a ground state with some fixed given magnetization m
% we require the following three integral equations:
%
% rho(x) + int_{-B}^{B} K(x,y) rho(y) dy = rho0(x)          (1)
% E0 = - Delta/4 - h*m + int_{-B}^{B} e0(x) rho(x) dx       (2)
%  m = 1/2 - int_{-B}^{B} rho(x) dx                         (3)
%
% Here E0 is the ground state energy, m is the magnetization, rho(x) is the Bethe root density function, 
% e0(x) are the bare energies (without the magnetic field interaction contribution) and B is the fermi repidity.
%
% We determine h from exploiting that (2) gives the lowest energy for a given h (or m), so dE0/dm =  - h + dI/dm = 0, 
% where I is the integral in (2). We thus have h = dI/dm.
% Using the chain rule we can further write dI/dm = dI/dB (dm/dB)^{-1} in terms of the fermi rapidity. These two
% derivatives we can actually calculate. For this we need to additionally calculate the derivative of rho w.r.t B:
% sigma(x) := drho(x)/dB
% This give rise to the following integral equation for sigma(x):
%
% sigma(x) + int_{-B}^{B} K(x,y) sigma(y) dy = - rho(B)*[K(x,B) + K(x,-B)]  (4)
%
% with (4) solved for sigma(x) we can then calculate
%
% dI/dB =  2*e0(B)*rho(B) + int_{-B}^{B} e0(x) sigma(x) dx
% dm/dB = -2*rho(B) - int_{-B}^{B} sigma(x) dx 
%
% and h = dI/dB (dm/dB)^{-1}

assert(Delta<1,'Delta<1');
assert(-0.5<m && m<0.5,'-0,5 < m < 0.5!');

if nargin<6 || isempty(verbose),verbose=false;end
if nargin<5 || isempty(tol),tol=1e-10;end
if nargin<4 || isempty(rho) || isempty(B)
    [~,rho,B]=fXXZGS_fixedm(Delta,m,tol,verbose);
else
    mtmp = 0.5-2*integral(rho,0,B,'AbsTol',tol,'RelTol',tol);
    if abs(mtmp-m)>tol,warning('MATLAB:fXXZfindh',...
            ['m from rho & x0 and given m are not consistent (',num2str(mtmp),' and ',num2str(m),', diff=',num2str(abs(mtmp-m),'%2.10e'),')']);end
end

if Delta == -1 % istropic gapless AFM
    Ker = @(x,n)(n./(pi*(x.^2 + n^2)));
    Kfun = @(x,y) (-Ker(x-y,2)- Ker(x+y,2));% formulate a symmetric kernel function, to only integrate from 0 to x
    e0fun = @(x) (-2*pi*Ker(x,1));
elseif Delta < -1 % anisotropic AFM (gapless for |m|>0)
    phi = acosh(-Delta);
    Ker = @(x,n) (phi*sinh(n*phi)./(2*pi*(cosh(n*phi) - cos(phi*x))));
    Kfun = @(x,y) (-Ker(x-y,2)-Ker(x+y,2)); % formulate a symmetric kernel function, to only integrate from 0 to x
    e0fun = @(x)(-2*pi*sinh(phi)*Ker(x,1)/phi);
elseif Delta > -1 && Delta < 1 % gapless anisotropic Luttinger Liquid phase
    gamma = acos(-Delta);
    Ker = @(x,n)(gamma*sin(n*gamma)./(2*pi*(cosh(gamma*x) -cos(n*gamma) )));
    Kfun = @(x,y)(-Ker(x-y,2)-Ker(x+y,2)); % formulate a symmetric kernel function, to only integrate from 0 to x
    e0fun = @(x)(-2*pi*sin(gamma)/gamma*Ker(x,1));
end

if m==0
    if Delta < -1 % determine critical hc for gapped AFM (above which |m|>0 and the system becomes gapless)
        m0 = findm0(pi/phi,tol);
        K0 = ellipke(m0);
        % Takahashi's eq. (4.35) is more than cryptic... It is very ambiguous as to what u' should be. In any case I
        % think the 1/phi is wrong, as it cancels with the phi in the exact expression for rho in (4.26)
        % Also, he defines H ~ -2hm, so use h-> 2h (i.e. add additional 2)
        h = 2*sinh(phi)*K0*sqrt(1-m0)/pi; % this is really the critical hc, above which m>0
    else
        h = 0;
    end
else
    drhodBstr = Fie(1,0,B,1,Kfun,@(x)(Kfun(x,B).*rho(B)),tol,10*tol);
    drhodB = @(x)(ntrpFie(drhodBstr,x));
    
    dedB = e0fun(B)*rho(B) + integral(@(x)(e0fun(x).*drhodB(x)),0,B,'RelTol',tol,'AbsTol',tol);
    dmdB = rho(B) + integral(drhodB,0,B,'RelTol',tol,'AbsTol',tol);
    
    h = -dedB/dmdB;
end
    
Eh = -h*m;
end

function [m0,F]=findm0(Q,tol)
fun = @(m)(Q*ellipke(1-m) - ellipke(m));

% ellipke(1) = Inf, therefor search interval must be [eps,1-eps] with eps small (use tol)
right = 1-tol;

% we know that fun is monotonically decreasing, so try left=0.5
% and if negative there, decrease until positive
left = 0.5;
while fun(left)<0,left=left/2;end

m0 = fzero(fun,[left,right]); % u is between 0 and 1
F = fun(m0);
end
