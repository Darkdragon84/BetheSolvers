clear;
clc;
close all;


%% parameters
Delta = -3;
h = 0;
tol = 1e-12;
make4 = false;
% make4 = true;
dP = 0.02;

assert(abs(h)<1-Delta,'|h|<1-Delta! otherwise trivial FM!');

opts = optimset('TolX',tol);
dig=ceil(-log10(tol));
frmt=['%2.',int2str(dig),'e'];


[E0,rho,B,Edr,m0]=fXXZGS_fixedh(Delta,h,tol,1);
disp(['E0=',num2str(E0,frmt),', m0=',num2str(m0,frmt)]);


pf = 2*integral(rho,0,B);
% disp([pf,0.5-m0,0.5-m0-pf]);
%% calculate energies and momenta

pmax = pf;
p0=(0:dP:pmax).';
% p0 = linspace(0,pmax,16).';
NP=numel(p0);
x=zeros(NP,1);

ct=NP/10;
for kk=2:NP
    xg = x(kk-1) + dP/(2*rho(x(kk-1))); % guess for next x, dk/dx = -2*rho(x)
    
    [xtmp,fval] = fzero(@(y)(p0(kk) - p0(kk-1) - 2*integral(rho,x(kk-1),y)),xg);
    x(kk)=xtmp;
    if kk>=ct
        disp([num2str(100*kk/NP),' %']);
        ct = ct + NP/10;
    end
end

e0 = -Edr(x);
p0 = [-pf;-flipud(p0(2:end));p0;pf];
% e0 = [0;flipud(e0(2:end));e0;0];
e0 = [-Edr(-B);flipud(e0(2:end));e0;-Edr(B)];

%% single particle
K1 = p0;
E1 = e0;

%% two particle

[Ktmp1,Ktmp2] = ndgrid(p0,p0);
[Etmp1,Etmp2] = ndgrid(e0,e0);

K2 = mod(Ktmp1 + Ktmp2 + 1,2) - 1;
E2 = Etmp1 + Etmp2;

K2 = K2(:);
E2 = E2(:);

%% three particle
[Ktmp1,Ktmp2] = ndgrid(K2,p0);
[Etmp1,Etmp2] = ndgrid(E2,e0);

K3 = mod(Ktmp1 + Ktmp2 + 1,2) - 1;
E3 = Etmp1 + Etmp2;

K3 = K3(:);
E3 = E3(:);

%% four particle
% make4 = numel(p0)^4<1e9;
if make4
    [Ktmp1,Ktmp2] = ndgrid(K3,p0);
    [Etmp1,Etmp2] = ndgrid(E3,e0);
    
    K4 = mod(Ktmp1 + Ktmp2 + 1,2) - 1;
    E4 = Etmp1 + Etmp2;
    
    K4 = K4(:);
    E4 = E4(:);
end


%%
fh1 = figure('color','w','position',[50,50,1200,800]);

ah11 = axes('position',[0.05,0.55,0.42,0.42]);
title(ah11,'elementary');
lh1 = line(K1,E1,'parent',ah11,'color','b');
xlabel('k/\pi');
ylabel('\DeltaE','rotation',0);
ylim(ah11,[0,1.1*max(E1)]);
xlim(ah11,[-1,1]);

ah12 = axes('position',[0.55,0.55,0.42,0.42]);
title(ah12,'two particle');
lh2 = line(K2,E2,'parent',ah12,'color','r','linestyle','none','marker','.');
xlabel('k/\pi');
ylabel('\DeltaE','rotation',0);
ylim(ah12,[0,1.1*max(E2)]);
xlim(ah12,[-1,1]);

ah21 = axes('position',[0.05,0.05,0.42,0.42]);
title(ah21,'three particle');
lh3 = line(K3,E3,'parent',ah21,'color','g','linestyle','none','marker','.');
xlabel('k/\pi');
ylabel('\DeltaE','rotation',0);
ylim(ah21,[0,1.1*max(E3)]);
xlim(ah21,[-1,1]);

if make4
    ah22 = axes('position',[0.55,0.05,0.42,0.42]);
    title(ah22,'four particle');
    lh4 = line(K4,E4,'parent',ah22,'color','k','linestyle','none','marker','.');
    xlabel('k/\pi');
    ylabel('\DeltaE','rotation',0);
    ylim(ah22,[0,1.1*max(E4)]);
    xlim(ah22,[-1,1]);
end




