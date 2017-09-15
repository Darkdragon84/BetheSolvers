clear;
close all;
clc;

%% parameters

U = 5;

n = 1;
% mu = -2;

dp = 0.01;

tol = 1e-10;
verbose = 1;
% usesave = false;
usesave = true;
datafldr = 'data';

opts = optimset('TolX',tol);
prec = ceil(-log10(tol));
frmt = sprintf('%%2.%ie',prec);

fixedn = false;
if exist('n','var')==1
    if exist('mu','var')==1,warning('n and mu given, using n');end
    fixedn = true;
    name = ['U',num2str(U),'_n',num2str(n),'_dp',num2str(dp),'_tol',num2str(tol),'.mat'];
    disp(['using fixed n=',num2str(n)]);
elseif exist('mu','var')==1
    name = ['U',num2str(U),'_mu',num2str(mu),'_dp',num2str(dp),'_tol',num2str(tol),'.mat'];
    disp(['using fixed mu=',num2str(mu)]);
else
    warning('no n or mu given, assuming half filling (n=1, mu=0)');
    mu = 0;
    name = ['U',num2str(U),'_mu0_dp',num2str(dp),'_tol',num2str(tol),'.mat'];
end

hf = false;
if fixedn
    if n<=0 || n>=2,error('0<n<2');end
    if n==1
        warning('half filled band, charge excitations gapped');
        hf = true;
    end
else
    muminus = 2-U/2-4*integral(@(x)(besselj(1,x)./(x.*(1+exp(0.5*x*U)))),0,Inf,'AbsTol',tol,'RelTol',tol);
    musat = -2-U/2;
    if mu < musat,error('completely empty');end
    if mu > abs(musat),error('completely filled');end
    if abs(mu) < abs(muminus)
        warning('half filled band, charge excitations gapped');
        hf = true;
    end
end

%% calculate excitations

pv = 0:dp:1;
Np = length(pv);

if usesave && exist([datafldr,'/',name],'file')==2
    F = load([datafldr,'/',name]);
    E0 = F.E0;
    k0 = F.k0;
    n = F.n;
    mu = F.mu;
    
    Edr_c = F.Edr_c;
    Edr_s = F.Edr_s;
    rho_c = F.rho_c;
    rho_s = F.rho_s;
    
    E_c = F.E_c;
    kv_c = F.kv_c;
    pv_c = F.pv_c;
    pfc = F.pfc;
    
    E_s = F.E_s;
    kv_s = F.kv_s;
    pv_s = F.pv_s;
    pfs = F.pfs;
    
else
    disp('ground state');
    if fixedn
        if ~hf
            [mu,Emu,~,k0,En] = fHUB_mu_from_n(U,n,tol,[],[],verbose);
            [E0,ntmp,k0,Edr_c,Edr_s,rho_c,rho_s] = fHUBGS_fixedmu(U,mu,tol/10,verbose,k0);
            disp('check energies:');
            disp(['E(n)  = ',num2str(En+Emu,frmt)]);
            disp(['E(mu) = ',num2str(E0,frmt)]);
            disp('check filling:');
            disp(['n0 = ',num2str(n,frmt)]);
            disp(['n  = ',num2str(ntmp,frmt)]);
        else
            [E0,~,k0,Edr_c,Edr_s,rho_c,rho_s] = fHUBGS_fixedmu(U,0,tol,verbose);
            mu = 0;
            n = 1;
        end
    else
        [E0,n,k0,Edr_c,Edr_s,rho_c,rho_s] = fHUBGS_fixedmu(U,mu,tol,verbose);
    end
    
    disp('done');
    
    % DRESSED MOMENTA (Essler 5.5.3)
    %
    % pc(k) = 2*pi*int_{0}^{k} rho_c(k') dk'
    % ps(Lam) = 2*pi*int_{0}^{Lam} rho_s(Lam') dLam'
    %
    % for ps it is actually slightly faster to evaluate (5.106).
    % as a side remark: z1(Inf) in (5.99) turns out to be z1(Inf) = pi*(n - ndown) for T=0.
    % For zero field we have ndown = n/2, so z1(Inf) = pi*n/2. 
    % Further, with int_{-Inf}^{Inf} sig1(Lam) dLam = ndown (and sig1 symmetric) we have 
    % int_{\Lam}^{Inf} sig1(Lam) dLam = ndown/2 - int_{0}^{\Lam} sig1(Lam) dLam and so finally
    % ps(Lam) = z1(Inf) - 2*pi*int_{\Lam}^{Inf} sig1(Lam) dLam = 2*pi*int_{0}^{Lam} rho_s(Lam') dLam'
    % I have no clue why (5.99) is so complicated when in fact it is formally exactly the same expression as pc(k)...
    % 
    
    if hf
%         pc = @(k)1/pi*(k + 2*integral(@(x) besselj(0,x).*sin(x.*sin(k))./(x.*(1 + exp(0.5*U*x))),0,Inf,'AbsTol',tol,'RelTol',tol,'ArrayValued',true));
%         ps = @(Lam)1/pi*integral(@(x) besselj(0,x).*sin(x.*Lam)./(x.*cosh(0.25*U*x)),0,Inf,'AbsTol',tol,'RelTol',tol,'ArrayValued',true);
        pc = @(k)1/pi*(k + 2*integral(@(x) besselj(0,x).*sin(x.*sin(k))./(x.*(1 + exp(0.5*U*x))),0,Inf,'AbsTol',tol,'RelTol',tol));
        ps = @(Lam)1/pi*integral(@(x) besselj(0,x).*sin(x.*Lam)./(x.*cosh(0.25*U*x)),0,Inf,'AbsTol',tol,'RelTol',tol);
        
        pfc = NaN; % charge excitations are gapped, there is no fermi momentum
        pfs = 0.5; % spin fermi momentum: ps(Inf)
    else
        pc = @(k) 2*integral(rho_c,0,k,'AbsTol',tol,'RelTol',tol);
        % for zero field the dressed spin momentum ps can also be written as (5.106), where the integral over Lam was
        % performed first
        ps = @(Lam) (n/2 - 2/pi*integral(@(k) rho_c(k).*(atan(exp(-2*pi/U*(Lam-sin(k)))) + atan(exp(-2*pi/U*(Lam+sin(k))))),0,k0,'AbsTol',tol,'RelTol',tol));
%         ps = @(Lam) 2/U*integral(@(k) rho_c(k).*integral(@(x) sech(2*pi/U*(x-sin(k))),0,Lam,'AbsTol',tol,'RelTol',tol,'ArrayValued',true),-k0,k0,'AbsTol',tol,'RelTol',tol);
        
        pfc = n; % charge fermi momentum: pc(k0)
        pfs = n/2; % spin fermi momentum: ps(Inf)
    end
    dpcdk = @(k) 2*rho_c(k);
    dpsdL = @(Lam) 2*rho_s(Lam);
    
    %% CHARGE EXCITATIONS
    disp('charge excitations');
    
    pv_c = pv;
    NPc = Np;
    
    E_c = zeros(1,NPc);
    kv_c = zeros(1,NPc);
    E_c(1) = Edr_c(0);
    kv_c(1) = 0;
    
    ct = 0;
    for kk = 1:NPc-1
        kold = kv_c(kk);
        knew = kold + dp/dpcdk(kold);
        knew = fzero(@(k)(pc(k) - pv_c(kk+1)),knew,opts);
        
        E_c(kk+1) = Edr_c(knew);
        kv_c(kk+1) = knew;
        
        if floor(10*kk/(NPc-1)) > ct
            disp([num2str(100*kk/(NPc-1),'%2.2f'),' % done']);
            ct = ct + 1;
        end
    end
    
    %% SPIN EXCITATIONS
    disp('spin excitations');
    
    pv_s = 0:dp:(pfs-tol); % to prevent that pv_s(end) is exactly pfs = ps(Inf), subtract tol
    NPs = length(pv_s);
    
    E_s = zeros(1,NPs);
    kv_s = zeros(1,NPs);
    E_s(1) = Edr_s(0);
    kv_s(1) = 0;
    
    ct = 0;
    for kk = 1:NPs-1
        kold = kv_s(kk);
        knew = kold + dp/dpsdL(kold);
        knew = fzero(@(k)(ps(k) - pv_s(kk+1)),knew,opts);
        
        E_s(kk+1) = Edr_s(knew);
        kv_s(kk+1) = knew;
        
        if floor(10*kk/(NPs-1)) > ct
            disp([num2str(100*kk/(NPs-1),'%2.2f'),' % done']);
            ct = ct + 1;
        end
    end
    
    if usesave
        save([datafldr,'/',name],'E0','k0','n','mu','Edr_c','Edr_s','rho_c','rho_s','E_c','kv_c','pv_c','pfc','E_s','kv_s','pv_s','pfs');
    end
end



%% elementary charge
if hf
%     pch = pv_c;
%     Ech = -E_c;
    
    pch = [-fliplr(pv_c),pv_c(2:end)];
    Ech = [-fliplr(E_c),-E_c(2:end)];
else
    pfcind = find(pv_c<pfc,1,'last');
    
%     pch = pv_c(1:pfcind);
%     Ech = -E_c(1:pfcind);
%     pcp = pv_c(pfcind+1:end);
%     Ecp =  E_c(pfcind+1:end);

    pch = [-pfc,-fliplr(pv_c(1:pfcind)),pv_c(2:pfcind),pfc];
    Ech = [0,-fliplr(E_c(1:pfcind)),-E_c(2:pfcind),0];
    
    pcp = [-fliplr(pv_c(pfcind+1:end)),-pfc,pfc,pv_c(pfcind+1:end)];
    Ecp = [fliplr(E_c(pfcind+1:end)),0,0,E_c(pfcind+1:end)];
end

%% elementary spin

if abs(pv_s(end)-pfs)>tol
    pse = [pv_s,pfs];
    Ese = -[E_s,0];
else
    pse = pv_s;
    Ese = -E_s;
end

pse = [-fliplr(pse(2:end)),pse];
Ese = [fliplr(Ese(2:end)),Ese];

%% spin-spin

[psstmp1,psstmp2] = ndgrid(pse,pse);
[Esstmp1,Esstmp2] = ndgrid(Ese,Ese);

% pss1 = mod(psstmp1 + psstmp2 + n + 1,2) - 1;
% pss2 = mod(psstmp1 + psstmp2 - n + 1,2) - 1;
pss = mod(psstmp1 + psstmp2 + 1,2) - 1;
Ess = Esstmp1 + Esstmp2;

% pss1 = pss1(:);
% pss2 = pss2(:);
pss = pss(:);
Ess = Ess(:);

%% hole-hole

[phhtmp1,phhtmp2] = ndgrid(pch,pch);
[Ehhtmp1,Ehhtmp2] = ndgrid(Ech,Ech);

phh = mod(phhtmp1 + phhtmp2 + 1,2) - 1;
Ehh = Ehhtmp1 + Ehhtmp2;

phh = phh(:);
Ehh = Ehh(:);

%% spin-hole

[pshtmp1,pshtmp2] = ndgrid(pse,pch);
[Eshtmp1,Eshtmp2] = ndgrid(Ese,Ech);

% psh1 = mod(pshtmp1 + pshtmp2 + n + 1,2) - 1;
% psh2 = mod(pshtmp1 + pshtmp2 - n + 1,2) - 1;
psh = mod(pshtmp1 + pshtmp2 + 1,2) - 1;
Esh = Eshtmp1 + Eshtmp2;

% psh1 = psh1(:);
% psh2 = psh2(:);
psh = psh(:);
Esh = Esh(:);

% the following excitations are only present away from half filling (where there are particle excitations)
if ~hf
    %% particle-particle
    [ppptmp1,ppptmp2] = ndgrid(pcp,pcp);
    [Epptmp1,Epptmp2] = ndgrid(Ecp,Ecp);
    
    ppp = mod(ppptmp1 + ppptmp2 + 1,2) - 1;
    Epp = Epptmp1 + Epptmp2;
    
    ppp = ppp(:);
    Epp = Epp(:);
    
    %% particle-hole
    [pphtmp1,pphtmp2] = ndgrid(pcp,pch);
    [Ephtmp1,Ephtmp2] = ndgrid(Ecp,Ech);
    
    pph = mod(pphtmp1 + pphtmp2 + 1,2) - 1;
    Eph = Ephtmp1 + Ephtmp2;
    
    pph = pph(:);
    Eph = Eph(:);
    %% spin-particle
    
    [psptmp1,psptmp2] = ndgrid(pse,pcp);
    [Esptmp1,Esptmp2] = ndgrid(Ese,Ecp);
    
    psp = mod(-psptmp1 + psptmp2 + 1,2) - 1;
    Esp = Esptmp1 + Esptmp2;
    
    psp = psp(:);
    Esp = Esp(:);
end
%% plot results

fh = figure('position',[200,200,1400,800]);
ah1 = axes;
ah2 = axes;
ah3 = axes;
ah4 = axes;

lhs = line(pse,Ese,'linestyle','none','marker','.','color','k','parent',ah1);
lhch = line(pch,Ech,'linestyle','none','marker','.','color','r','parent',ah1);
lg1 = {'spin','hole'};

lhchh = line(phh,Ehh,'linestyle','none','marker','.','color','r','parent',ah2);
lg2 = {'hole-hole'};

% lhss1 = line(pss1,Ess,'linestyle','none','marker','.','color','r','parent',ah3);
% lhss2 = line(pss2,Ess,'linestyle','none','marker','.','color','b','parent',ah3);
% lg3 = {'spin-spin 1','spin-spin 2'};
lhss1 = line(pss,Ess,'linestyle','none','marker','.','color','r','parent',ah3);
lg3 = {'spin-spin'};


% lhsh1 = line(psh1,Esh,'linestyle','none','marker','.','color','r','parent',ah4);
% lhsh2 = line(psh2,Esh,'linestyle','none','marker','.','color','b','parent',ah4);
% lg4 = {'spin-hole 1','spin-hole 2'};
lhsh1 = line(psh,Esh,'linestyle','none','marker','.','color','r','parent',ah4);
lg4 = {'spin-hole'};

if ~hf
    lhcp = line(pcp,Ecp,'linestyle','none','marker','.','color','b','parent',ah1);
    lg1 = [lg1,{'particle'}];
    
    lhcpp = line(ppp,Epp,'linestyle','none','marker','.','color','b','parent',ah2);
    lhcph = line(pph,Eph,'linestyle','none','marker','.','color','k','parent',ah2);
    lg2 = [lg2,{'particle-particle','particle-hole'}];
    
    lhsp = line(psp,Esp,'linestyle','none','marker','.','color','k','parent',ah4);
    lg4 = [lg4,{'spin-particle'}];
end
%% plot fixes
lm = 0.04;
bm = 0.04;
wi = 0.45;
he = 0.45;

y1 = get(ah1,'ylim');
y2 = get(ah2,'ylim');
y3 = get(ah3,'ylim');
y4 = get(ah4,'ylim');
ylim(ah1,[0,y1(2)]);
ylim(ah2,[0,y2(2)]);
ylim(ah3,[0,y3(2)]);
ylim(ah4,[0,y4(2)]);


set(ah1,'position',[lm,0.5+bm,wi,he]);
set(ah2,'position',[0.5+lm,0.5+bm,wi,he]);
set(ah3,'position',[lm,bm,wi,he]);
set(ah4,'position',[0.5+lm,bm,wi,he]);

legend(ah1,lg1);
legend(ah2,lg2);
legend(ah3,lg3);
legend(ah4,lg4);