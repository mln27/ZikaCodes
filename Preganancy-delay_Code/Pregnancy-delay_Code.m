%Zika pregnancy delay model
%Fixed parameter values:
global Const;

bH = 40*7;
HR = 1/3;
Const.BH2 = 1/((1-HR)*bH);
Const.BH1 = 1/(HR*bH);
Const.DHM = 1/(364*11.8);
Const.DHW = 1/(29.8*364);
Const.Gratio = 0.492;
Const.TransM1 = 1/(364*25);
Const.TransM2 = 1/(35*364);
Const.TransW1 = 1/(21.7*364);
Const.TransW2 = 1/(26.8*364);
Const.P = 1.9*Const.TransW2;

 N0 = [0.20 0.15 0.24 (1-Const.P*bH)*0.238 0.238*Const.P*bH*(1-HR) 0.238*Const.P*bH*HR 0.052 0.12];
Const.Data = ([54 38 18	24 35 51 64	137	140	175	276	353	587	1130 1232 1340 1699 1652 2002 1758 1687 2575 3102 4277 5721 6228 5396 4670 3855 3758 3223 2682 2397 2504])';


%Vector-Host transmission model
function dy = ZikadiffMass(t,y,par)

BetaM = par(1); BetaH = par(2); C = par(3); Alpha = par(4); Gamma = par(5); TauM = par(6); DO = par(7); 
Theta = par(end-1);
global Const;

%Seasonal variation of mosquitoes birth rate
 TPY = DO*(1+(1/4)*sin((t)*pi/(182)+1.0472));

%Model variables Mosquitoes
OS = y(1); OE = y(2); OI = y(3);
%Model variables humans
%Pre-reproductive age men
SM1 = y(4); EM1 = y(5); IM1 = y(6); RM1 = y(7);
%Pre-reproductive age women
SW1 = y(8); EW1 = y(9); IW1 = y(10); RW1 = y(11); 
%Reproductive age men
SM2 = y(12); EM2 = y(13); IM2 = y(14); RM2 = y(15);
%Reproductive age women, non pregnant
SW2 = y(16); EW2 = y(17); IW2 = y(18); RW2 = y(19);
%Reproductive age women, first trimester pregnancy
SWPR2 = y(20); EWPR2 = y(21); IWPR2 = y(22); RWPR2 = y(23);
%Reproductive age women, late pregnancy
SWP2 = y(24); EWP2 = y(25); IWP2 = y(26); RWP2 = y(27);
%Post-reproductive age men
SM3 = y(28); EM3 = y(29); IM3 = y(30); RM3 = y(31);
%Post-reproductive age women
SW3 = y(32); EW3 = y(33); IW3 = y(34); RW3 = y(35); 
%comunaltive infections
CI = y(36);
%comunaltive women infections
CW = y(37);
%comunaltive reproductive age infections
CR = y(38);

%total population size
NH = SM1+SW1+SM2+SW2+SWP2+SWPR2+SM3+SW3+...
    EM1+EW1+EM2+EW2+EWP2+EWPR2+EM3+EW3+...
    IM1+IW1+IM2+IW2+IWP2+IWPR2+IM3+IW3+...
    RM1+RW1+RM2+RW2+RWP2+RWPR2+RM3+RW3;

%humans force ofinfection
LambdaM = BetaM.*sum(C.*(IM1+IW1+IM2+Theta*(IW2+IWP2+IWPR2)+IM3+IW3))./sum(NH);
%mosquitoes force of infection
Lambda1 = C.*(BetaH)*(OI);

dx1 = TPY.*(OS+OE+OI)-(LambdaM+DO).*OS;
dx2 = LambdaM.*OS-(TauM+DO).*OE;
dx3 = TauM.*OE-(DO).*OI;
%pre-reproductive males and females
dx4 = Const.BH2*Const.Gratio*(SWP2+EWP2+IWP2+RWP2)-(Const.TransM1+Lambda1)*SM1;
dx5 = -Const.TransM1*EM1+(Lambda1)*SM1-Alpha*EM1;
dx6 = -Const.TransM1*IM1+Alpha*EM1-Gamma*IM1;
dx7 = -Const.TransM1*RM1+Gamma*IM1;
dx8 = Const.BH2*(1-Const.Gratio)*(SWP2+EWP2+IWP2+RWP2)-(Const.TransW1+Lambda1)*SW1;
dx9 = -Const.TransW1*EW1+(Lambda1)*SW1-Alpha*EW1;
dx10 = -Const.TransW1.*IW1+Alpha*EW1-Gamma*IW1;
dx11 = -Const.TransW1*RW1+Gamma*IW1;
%Reproductive males
dx12 = Const.TransM1*SM1-(Const.TransM2+Lambda1)*SM2;
dx13 = Const.TransM1*EM1-Const.TransM2*EM2+(Lambda1)*SM2-Alpha*EM2;
dx14 = Const.TransM1*IM1-Const.TransM2*IM2+Alpha*EM2-Gamma*IM2;
dx15 = Const.TransM1*RM1-Const.TransM2*RM2+Gamma.*IM2;
%Reproductive females not pregnant
dx16 = Const.TransW1*SW1-(Const.TransW2+Theta*Lambda1)*SW2+Const.BH2*SWP2-Const.P*SW2;
dx17 = Const.TransW1*EW1-Const.TransW2*EW2+(Theta*Lambda1)*SW2-Alpha*EW2+Const.BH2*EWP2-Const.P*EW2;
dx18 = Const.TransW1*IW1-Const.TransW2*IW2+Alpha*EW2-Gamma.*IW2+Const.BH2*IWP2-Const.P*IW2;
dx19 = Const.TransW1*RW1-Const.TransW2*RW2+Gamma.*IW2+Const.BH2*RWP2-Const.P*RW2;
%Reproductive females pregnant at risk
dx20 = Const.P*SW2-(Const.BH1+Theta*Lambda1)*SWPR2;
dx21 = Const.P*EW2+(Theta*Lambda1)*SWPR2-(Const.BH1+Alpha)*EWPR2;
dx22 = Const.P*IW2-Const.BH1*IWPR2+Alpha*EWPR2-Gamma*IWPR2;
dx23 = Const.P*RW2-Const.BH1*RWPR2+Gamma*IWPR2;
%Reproductive females pregnant not at risk
dx24 = Const.BH1*SWPR2-(Const.BH2+Theta*Lambda1)*SWP2;
dx25 = Const.BH1*EWPR2+(Theta*Lambda1)*SWP2-(Alpha+Const.BH2)*EWP2;
dx26 = Const.BH1*IWPR2+Alpha*EWP2-Gamma*IWP2-Const.BH2*IWP2;
dx27 = Const.BH1*RWPR2+Gamma*IWP2-Const.BH2*RWP2;
%Post-reprodutive males and females
dx28 = Const.TransM2*SM2-(Const.DHM+Lambda1)*SM3;
dx29 = Const.TransM2*EM2-Const.DHM*EM3+(Lambda1)*SM3-Alpha*EM3;
dx30 = Const.TransM2*IM2-Const.DHM*IM3+Alpha*EM3-Gamma*IM3;
dx31 = Const.TransM2*RM2-Const.DHM*RM3+Gamma*IM3;
dx32 = Const.TransW2*SW2-(Const.DHW+Lambda1)*SW3;
dx33 = Const.TransW2*EW2-Const.DHW*EW3+(Lambda1)*SW3-Alpha*EW3;
dx34 = Const.TransW2*IW2-Const.DHW*IW3+Alpha*EW3-Gamma*IW3;
dx35 = Const.TransW2*RW2-Const.DHW*RW3+Gamma*IW3;

%comulative infected cases
dx36 = (Alpha)*(EM1+EW1+EM2+EW2+EM3+EW3+EWP2+EWPR2);
%comulative infected cases among women
dx37 = Alpha*(EW1+EW2+EW3+EWP2+EWPR2);
%cumulative infected cases among people of reproductive age
dx38 = Alpha*(EW2+EWP2+EWPR2+EM2);

dy = [dx1 dx2 dx3 dx4 dx5 dx6 dx7 dx8 dx9 dx10 dx11 dx12...
    dx13 dx14 dx15 dx16 dx17 dx18 dx19 dx20 dx21 dx22 dx23...
    dx24 dx25 dx26 dx27 dx28 dx29 dx30 dx31 dx32 dx33...
    dx34 dx35 dx36 dx37 dx38]';

%Prior distributions
function F = ZikainitsNoAge(n)

F = zeros(n,14);

for i = 1:n

BetaM = rand(1)*(1-0.001)+0.001;
BetaH  = rand(1)*(0.75-0.001)+0.001;
C = rand(1)*(1-0.33)+0.33;
Gamma = 1/(rand(1)*(22-3)+3);
Alpha = 1/(rand(1)*(17-6)+6);
TauM = 1/gamrnd(12,1,1);
rep = rand(1)*(0.27-0.09)+0.09;
ART = (rand(1)*(0.77-0.5)+0.5);
DO = 1/(rand(1)*(42-8)+8);
disp = rand(1)*1000;
Ih0 = randi(20);
Iv0 = rand(1)*0.002;
Theta = rand(1)*(10-1)+1;
Tmax = (randi(7)+38)*7;

G(i,:) = [BetaM,BetaH,C,Alpha,Gamma,TauM,DO,rep,ART,disp,Ih0,Iv0,Theta,Tmax];

end
F = G;

%Likelihood function
function [F] = ZikalikelihoodMass(para,c0)

global Const;

tmax=para(end);
t= [0:tmax];

curve=ode15s(@ZikadiffMass,t,c0,[],para);
curve=deval(curve,t)';
CI = (curve(:,36));
CW = curve(:,37);
CR = curve(:,38);

I = tmax/7;
Zikainc = zeros(1,I);
for j = 1:I
Zikainc(j) = CI((j)*7)-CI((j-1)*7+1);
end

Zikaincdata = Const.Data';

L1=length(Zikaincdata);
LHI=zeros(1,L1);

for i=1:L1
  if(Zikainc(end-(i-1))>Zikaincdata(end-(i-1)))
    LHI(1,i)= nbinpdf(Zikaincdata(end-(i-1)),para(10),para(10)/(para(10)+ceil(para(8).*Zikainc(end-(i-1)))));
    else
         LHI(1,i)=0;
     end
end
F{1} = prod(LHI(1,:))*betapdf(CW(end-7)/CI(end-7),12615,7224)*betapdf(CR(end-7)/CI(end-7),12459,7380);


%%Bayesian mdeling