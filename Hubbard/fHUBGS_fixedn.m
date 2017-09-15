function [E0,k0,rho_c] = fHUBGS_fixedn(U,n,tol,verbose)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates Ground state energy for the fermionic Hubbard model in the TD
% limit for a given particle density at zero magnetic field
%
% use particle-hole symmetric Hamiltonian (n -> n-1/2):
% (nu-1/2)(nd-1/2) = nu*nd - 1/2(nu+nd) + 1/4
%
% H = - T + U \sum_j (nu_j - 1/2)(nd_j - 1/2) - mu \sum_j (nu_j + nd_j) - h/2 \sum_j (nu_j - nd_j)
%   = - T + U \sum_j (nu_j)(nd_j) - (mu + U/2) \sum_j (nu_j + nd_j) + L*U/4
%
% where T is the hopping term [T = - \sum_(js) (c*_(sj)c_(sj+1) + c*_(sj+1)c_(sj) ], n = nu + nd 
% and L->infty
%
% therefore: correct energies to E0 -> E0 - U/2*n + U/4
%
% input: U   ... interaction strength
%        n  ... particle density
%        tol ... tolerance for numerical integrations
%
%   Refs.: [1] E. Lieb, F. Wu, PRL 20, 1445 (1968)
%          [2] H. Frahm, V. Korepin, PRB 42, 10553 (1990)
%          [3] F. Essler et al.: The One-Dimensional Hubbard Model, Cambridge (2005)
%
% Valentin Stauber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3||isempty(tol),tol=1e-10;end;
if nargin<4||isempty(verbose),verbose=false;end;

assert(n>=0 && n<=2,'0<=n<= 2!');
frmt=['%2.',int2str(ceil(-log10(tol))),'e'];

if n==1
    %% HALF FILLING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % n=1 <==> k0 = pi -> Bethe Equations can be solved by Fourier transform (from [1])
    
    efun = @(x)(besselj(0,x).*besselj(1,x)./(x.*(1+exp(0.5*x*U))));
    rhofun = @(x,k)(cos(x.*sin(k)).*besselj(0,x)./(1 + exp(0.5*x*U)));
    
    E0 = -4*integral(efun,0,Inf,'AbsTol',tol,'RelTol',tol) - U/4;
    k0 = pi;
%     rho_c = @(k)(1/pi*(1/2 + cos(k).*integral(@(x)rhofun(x,k),0,Inf,'AbsTol',tol,'RelTol',tol,'ArrayValued',true)));
    rho_c = @(k)(1/pi*(1/2 + cos(k).*integral(@(x)rhofun(x,k),0,Inf,'AbsTol',tol,'RelTol',tol)));
elseif n>1
    % due to symmetries it is enough to consider nu<=0.5 and nd<=0.5, i.e. n<=1
    % (c.f. [1] & [3])
    [E0,k0,rho_c] = fHUBGS_fixedn(U,2-n,tol,verbose);
elseif n<1 
    % general n, but zero magnetization: Essler 5.5.4
    Rfun = @(x)(1/pi * integral(@(y)(cos(y.*x)./(1 + exp(0.5*U*y))),0,Inf,'AbsTol',tol,'RelTol',tol,'ArrayValued',true));
    Kfun = @(x,y) cos(x).*(Rfun(sin(x) - sin(y)) + Rfun(sin(x) + sin(y))); % symmetric Kernel function, for less numerical effort
%     Rfun = @(x)(1/pi * integral(@(y)(cos(y.*x)./(1 + exp(0.5*U*y))),0,Inf,'AbsTol',tol,'RelTol',tol));
%     Kfun = @(x,y) cos(x).*(arrayfun(Rfun,sin(x) - sin(y)) + arrayfun(Rfun,sin(x) + sin(y))); % symmetric Kernel function, for less numerical effort
    Gfun = @(x)(ones(size(x))/(2*pi));
    
    kold = 0;
    knew = pi*n;
    
    err.Ferr=0;
    err.Ierr=0;
    
    ct = 1;
    while abs(knew - kold)>tol
        kold = knew;
        if verbose,fprintf('iteration %u: ',ct);end
        % Fredholm
        if verbose,fprintf('Fredholm: ');end
        [rhocstr,err.Ferr] = Fie(1,0,kold,1,Kfun,Gfun,tol,10*tol);
        if verbose,fprintf('done, zero search: ');end
        % zero search
        rho_c = @(x)(ntrpFie(rhocstr,abs(x)));
        [knew,err.Ierr] = fzero(@(k)(n - 2*integral(rho_c,0,k)),kold);
        if verbose,disp(['done, k(',int2str(ct),'): ',num2str(knew,frmt),', dk: ',num2str(knew-kold,frmt)]);end
        ct=ct+1;
    end
    k0 = knew;
    
    efun = @(x)(-2*cos(x).*rho_c(x));
    E0 = 2*integral(efun,0,k0,'AbsTol',tol,'RelTol',tol) - U/2*n + U/4;
end


if verbose,disp(['E0(U=',num2str(U),';n=',num2str(n),'): ',num2str(E0,frmt)]);end
end

