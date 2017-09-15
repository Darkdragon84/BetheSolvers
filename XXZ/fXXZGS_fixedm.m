function [E0,rho,B,Edr]=fXXZGS_fixedm(Delta,m,tol,verbose)
% solves Bethe ansatz equations in the TD limit for XXZ Hamiltonian
%
%   H = -\sum_j Sx_{j} Sx_{j+1} +  Sy_{j} Sy_{j+1} + Delta Sz_{j} Sz_{j+1} 
%
% with fixed magnetization -1/2 < m < 1/2 and Delta < 1 (gapped FM trivial)
%
% E0   ... ground state energy density
% rho  ... Bethe root distribution function. rho0 is such that m = 1/2 - int_{-B}^{B} rho(x) dx
% B    ... fermi rapidity: max. Bethe spectral parameter, i.e. max. occupied Bethe root (basically the integration limit for all integral equations)
% Edr  ... dressed energy, the elementary excitations are then given by -Edr(x) and have momentum p(x) = 2*pi*int_{x}^{B} rho(x) dx
%
% we use the following two integral equations:
% 
% rho(x) + int_{-B}^{B} K(x,y) rho(y) dy = rho0(x)  (1)
% m = 1/2 - int_{-B}^{B} rho(x) dx                  (2)
% 
% the task is to calculate rho(x) and B from these two equations, and we do so by starting with some initial B and
% alternately solve (1) for rho(x) and (2) for B.
%
% With rho(x) and B known we can calculate the ground state energy
% E0 = -Delta/4 + int_{-B}^{B} e0(x) rho(x) dx  (hwew e0(x) is the bare energy without the magnetic field contribution)
%
% We can also calculate the dressed energy Edr(x), which satisfies
%
%   Edr(x) + int_{-B}^{B} K(x,y) Edr(y) dy = e0(x)
%
% Watch out however that e0(x) is missing the correct magnetic field, so Edr(x) is only correct up to some constant
% offset!

if nargin<2||isempty(m),m=0;end;
if nargin<3||isempty(tol),tol=1e-10;end;
if nargin<4||isempty(verbose),verbose=false;end;

assert(-0.5<m && m<0.5,'-0.5 < m < 0.5');

frmt=['%2.',int2str(ceil(-log10(tol))),'e'];

%% calculate ground state energy and root density
if Delta == -1 % istropic gapless AFM (from Takahashi)
    
    % watch out: Fie solves lam*rho(x) - integral(Kfun(x,y)*rho(y),-x0,x0) = Gfun(x), beware of the minus
    Ker = @(x,n)(n./(pi*(x.^2 + n^2)));
    Kfun = @(x,y) (-Ker(x-y,2)-Ker(x+y,2));% formulate a symmetric kernel function, to only integrate from 0 to x
    Gfun = @(x) Ker(x,1);
    e0fun = @(x) (-2*pi*Ker(x,1));
    
    if m==0
        E0 = 0.25 - log(2);
        rho = @(x)(sech(pi*x/2)/4);
        B = Inf;
    else
        xold = 0;
        xnew = 1;
        
        err.Ferr=0;
        err.Ierr=0;
        
        ct = 1;
        while abs(xnew - xold)>tol
            xold = xnew;
            if verbose,fprintf('iteration %u: ',ct);end
            % Fredholm
            if verbose,fprintf('Fredholm: ');end
            [rhostr,err.Ferr] = Fie(1,0,xold,1,Kfun,Gfun,tol,10*tol);
            if verbose,fprintf('done, zero search: ');end
            % zero search
            rho = @(x)(ntrpFie(rhostr,abs(x)));
            [xnew,err.Ierr] = fzero(@(x)(m - 0.5 + 2*integral(rho,0,x)),xold);
            if verbose,disp(['done, x(',int2str(ct),'): ',num2str(xnew,frmt),', dx: ',num2str(xnew-xold,frmt)]);end
            ct=ct+1;
        end
        B = xnew;
        
        E0 = 1/4 + 2*integral(@(x)(e0fun(x).*rho(x)),0,B);
    end
elseif Delta<-1 % gapped AFM
    phi=acosh(-Delta);
    Q=pi/phi;
    
    % watch out: Fie solves lam*rho(x) - integral(Kfun(x,y)*rho(y),-x0,x0) = Gfun(x), beware of the minus
    Ker = @(x,n) (phi*sinh(n*phi)./(2*pi*(cosh(n*phi) - cos(phi*x))));
    Kfun = @(x,y) (-Ker(x-y,2)-Ker(x+y,2)); % formulate a symmetric kernel function, to only integrate from 0 to x
    Gfun = @(x) Ker(x,1);
    e0fun = @(x)(-2*pi*sinh(phi)*Ker(x,1)/phi);
    
    if m==0
        %% ground state energy (from Takahashi)
        sumfac = 0;
        fac=2*sinh(phi);
        nn = 1;
        summand=1/(exp(2*nn*phi)+1);
        while abs(fac*summand)>tol
            sumfac = sumfac + fac*summand;
            nn = nn + 1;
            summand=1/(exp(2*nn*phi)+1);
        end
        E0 = -Delta/4 - sinh(phi)/2 - sumfac;
        %% root density
        % Takahashi defines elliptic integrals in terms of u=sqrt(m) (or u^2=m)
        % i.e. K(u) = int_0^pi/2 dx (1-u^2 sin^2 x)^(-1/2)
        % matlab uses m=u^2
        m0 = findm0(Q,tol); % u0 is between 0 and 1
        K0 = ellipke(m0);
        rho = @(x)(rho_h0_fun(x,K0,m0,Q));
        B = Q;
    else
        xold = 0;
        xnew = (1-2*m)*Q;
        
        err.Ferr=0;
        err.Ierr=0;
        
        ct = 1;
        while abs(xnew - xold)>tol
            xold = xnew;
            if verbose,fprintf('iteration %u: ',ct);end
            % Fredholm
            if verbose,fprintf('Fredholm: ');end
            [rhostr,err.Ferr] = Fie(1,0,xold,1,Kfun,Gfun,tol,10*tol);
            if verbose,fprintf('done, zero search: ');end
            % zero search
            rho = @(x)(ntrpFie(rhostr,abs(x)));
            [xnew,err.Ierr] = fzero(@(x)(m - 0.5 + 2*integral(rho,0,x)),xold);
            if verbose,disp(['done, x(',int2str(ct),'): ',num2str(xnew,frmt),', dx: ',num2str(xnew-xold,frmt)]);end
            ct=ct+1;
        end
        B = xnew;
        
        E0 = -Delta/4 + 2*integral(@(x)(e0fun(x).*rho(x)),0,B,'RelTol',tol,'AbsTol',tol);
    end
elseif -1<Delta && Delta<1 % Luttinger Liquid (gapless) phase, both AFM and FM
    gamma = acos(-Delta);
    p0 = pi/gamma;
    
    % watch out: Fie solves lam*rho(x) - integral(Kfun(x,y)*rho(y),-x0,x0) = Gfun(x), beware of the minus
    Ker = @(x,n)(gamma*sin(n*gamma)./(2*pi*(cosh(gamma*x) -cos(n*gamma) )));
    Kfun = @(x,y)(-Ker(x-y,2)-Ker(x+y,2)); % formulate a symmetric kernel function, to only integrate from 0 to x
    Gfun = @(x)(Ker(x,1));
    e0fun = @(x)(-2*pi*sin(gamma)/gamma*Ker(x,1));
    
    if m==0
        %% ground state energy and root density (from Cloiseaux/Gaudin and Takahashi)

        % from Takahashi or Cloiseaux/Gaudin (their (53) can be misread, the cosh is in the denominator!!)
        rho = @(x)(sech(pi*x/2)/4); 
        % from Cloiseaux/Gaudin (their Theta=gamma, rho=-Delta, epsilon = e + Delta/4)
        % Takahashi's (4.44) is wrong, the integral must be divided by 2!
        efun = @(x)(1-tanh(x*gamma)./tanh(x*pi)); 
        
        E0 = -Delta/4 - sin(gamma)*integral(efun,0,Inf,'RelTol',tol,'AbsTol',tol); % integral can deal with infinite integration limits!!
        B = Inf;
    else
        xold = 0;
        xnew = p0;
        
        err.Ferr=0;
        err.Ierr=0;
        
        ct = 1;
        while abs(xnew - xold)>tol
            xold = xnew;
            if verbose,fprintf('iteration %u: ',ct);end
            % Fredholm
            if verbose,fprintf('Fredholm: ');end
            [rhostr,err.Ferr] = Fie(1,0,xold,1,Kfun,Gfun,tol,10*tol);
            if verbose,fprintf('done, zero search: ');end
            % zero search
            rho = @(x)(ntrpFie(rhostr,abs(x)));
            [xnew,err.Ierr] = fzero(@(x)(m - 0.5 + 2*integral(rho,0,x)),xold);
            if verbose,disp(['done, x(',int2str(ct),'): ',num2str(xnew,frmt),', dx: ',num2str(xnew-xold,frmt)]);end
            ct=ct+1;
        end
        B = xnew;
        
        E0 = -Delta/4 + 2*integral(@(x)(e0fun(x).*rho(x)),0,B,'RelTol',tol,'AbsTol',tol);
    end
else
    error('not implemented yet');
end

Edrstr = Fie(1,0,B,1,Kfun,e0fun,tol,10*tol);
Edr = @(x)(ntrpFie(Edrstr,abs(x)));

end

function rho=rho_h0_fun(x,K0,m0,Q)
[~,~,dn]=ellipj(K0.*x./Q,m0); % this is fine and checked now!!
rho=K0./(2*pi*Q).*dn;  % this is fine and checked now!!
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


