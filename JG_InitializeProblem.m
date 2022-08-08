function [t,x] = JG_InitializeProblem (Lconc,Rconc,tf_input,calcium,cBCL2,newXIAP)

% Function Definitions
% k#f is the forward binding rate constant
% k#b is the backwards binding rate constant
% Lconc is the concentration of ligand (TRAIL)
% Rconc is the concentration of receptor (DR4)
% tf_input is the end time of integration
%calcium in uM

% t is the time vector for numerical integration
% x is the matrix of reactant values, where each column is a different
% reactant

%% Initial Conditions

% Non-zero initial conditions (in molecules per cell)

calcium = calcium*10^-6;       %convert calcium to M
calcium = calcium*(10E12/17);  %convert calcium to #/CC

flip  = 1E2;  % Flip
pC8   = 2E4;  % procaspase-8 (pro-C8)
Bar   = 1E3;  % Bifunctional apoptosis regulator
pC3   = 1E4;  % procaspase-3 (pro-C3)
pC6   = 1E4;  % procaspase-6 (pro-C6)
PARP  = 1E6;  % C3* substrate
Bid   = 4E4;  % Bid
Bax   = 1E5;  % Bax
M     = 5E5;  % mitochondrial binding sites for activated Bax
CytoC = 5E5;  % cytochrome c
Smac  = 1E5;  % Smac
pC9   = 1E5;  % procaspase-9 (pro-C9)
Apaf  = 1E5;  % Apaf-1
CP    = 1E5;  % calpain
pC12  = 0;    % procaspase-12 (pro-C12)
CS    = 1E5;  % calpastatin
Bcl2c = cBCL2;  % cytosolic Bcl-2
Bcl2  = 2E4;  % mitochondrial Bcl-2
XIAP  = newXIAP;  % XIAP


transloc=.01;% rate of translocation between the cytosolic and mitochondrial compartments
v= .1; % mitochondria compartment volume/cell volume

% Initialize the full vector of initial conditions (IC)
IC=zeros(75,1);
IC(1) =Lconc;    IC(2) =Rconc;    IC(5) =flip; IC(7) =pC8;  IC(10)=Bar;
IC(12)=pC3;  IC(15)=pC6;  IC(19)=XIAP; IC(21)=PARP; IC(24)=Bid;
IC(27)=Bcl2c; IC(29)=Bax;  IC(33)=Bcl2; IC(39)=M;    IC(42)=CytoC;
IC(45)=Smac; IC(49)=Apaf; IC(52)=pC9; IC(59) = calcium; IC(60) = CP;
IC(62)=CS; IC(68)=pC12;

%% Rate Constants

% L + R <--> L:R  ---> R*
k(1)=4E-6; k_(1)= 5E-3; kc(1)=1E-5;

% flip + DISC <-->  flip:DISC
k(2)=5E-7; k_(2)=5E-4;

% pC8 + DISC <--> DISC:pC8 --> C8 + DISC
k(3)=5E-7; k_(3)=5E-4; kc(3)=1;

% C8 + BAR <--> BAR:C8
k(4)=5E-8; k_(4)=5E-4;

% pC3 + C8 <--> pC3:C8 --> C3 + C8
k(5)=5E-8; k_(5)=5E-4; kc(5)=1;

% pC6 + C3 <--> pC6:C3 --> C6 + C3
k(6)=5E-7; k_(6)=5E-4; kc(6)=1;

% pC8 + C6 <--> pC8:C6 --> C8 + C6
k(7)=6E-9; k_(7)=5E-4; kc(7)=1;

% XIAP + C3 <--> XIAP:C3 --> XIAP + C3_U
k(8)=2E-7; k_(8)=1E-4;  kc(8)=1E-4;

% PARP + C3 <--> PARP:C3 --> CPARP + C3
k(9)=1E-7; k_(9)=1E-3; kc(9)=1;

% Bid + C8 <--> Bid:C8 --> tBid + C8
k(10)=5E-7; k_(10)=1E-4; kc(10)=1;

% tBid + Bcl2c <-->  tBid:Bcl2c
k(11)=1E-8; k_(11)=1E-4;

% Bax + tBid <--> Bax:tBid --> aBax + tBid
k(12)=1E-9; k_(12)=1E-4; kc(12)=1;

% aBax <-->  MBax
k(13)=transloc; k_(13)=transloc;

% MBax + Bcl2 <-->  MBax:Bcl2
k(14)=1E-7; k_(14)=1E-4;

% MBax + MBax <-->  MBax:MBax == Bax2
k(15)=1E-7; k_(15)=1E-4;

% Bax2 + Bcl2 <-->  MBax2:Bcl2
k(16)=1E-7; k_(16)=1E-4;

% Bax2 + Bax2 <-->  Bax2:Bax2 == Bax4
k(17)=1E-7; k_(17)=1E-4;

% Bax4 + Bcl2 <-->  MBax4:Bcl2
k(18)=1E-7; k_(18)=1E-4;

% Bax4 + Mit0 <-->  Bax4:Mito -->  AMito
k(19)=1E-7; k_(19)=1E-4; kc(19)=1;

% AMit0 + mCtoC <-->  AMito:mCytoC --> AMito + ACytoC
k(20)=2E-7; k_(20)=1E-4; kc(20)=10;

% AMit0 + mSMac <-->  AMito:mSmac --> AMito + ASMAC
k(21)=2E-7; k_(21)=1E-4; kc(21)=10;

% ACytoC <-->  cCytoC
k(22)=transloc; k_(22)=transloc;

% Apaf + cCytoC <-->  Apaf:cCytoC --> Apaf*
k(23)=5E-8; k_(23)=1E-4; kc(23)=1;

% Apaf* + Procasp9 <-->  Apoptosome
k(24)=5E-9; k_(24)=1E-4;

% Apop + pCasp3 <-->  Apop:cCasp3 --> Apop + Casp3
k(25)=5E-8; k_(25)=1E-4; kc(25)=1;

% ASmac <-->  cSmac
k(26)=transloc; k_(26)=transloc;

% Apop + XIAP <-->  Apop:XIAP
k(27)=2E-7; k_(27)=1E-4;

% cSmac + XIAP <-->  cSmac:XIAP
k(28)=7E-7; k_(28)=1E-4;

% calcium + calpain <--> calpainC
k(29)=5E-10; k_(29)=1E-4;

% calcium + calpastatin <--> calpastatin*
k(30)=5E-10; k_(30)=1E-4;

% calpastatin* + calpainC <--> calpastatin*:calpainC --> calpain_blocked
k(31)=3E-7; k_(31)=1E-4; kc(31)=1;

% C3 + calpain_blocked <--> C3:calpain_blocked --> C3 + calpain*
k(32)=1E-7; k_(32)=1E-4; kc(32)=1;

%calpain* + pCasp12 <--> calpain*:pCasp12 --> calpain* + Casp12
k(33)=1E-7; k_(33)=1E-4; kc(33)=1;

%Casp12 + pCasp9 <--> Casp12:pCasp9 --> Casp12 + Casp9
k(34)=5E-9; k_(34)=1E-4; kc(34)=1;

%calpain* + Bid <--> calpain*:Bid --> calpain* + tBid
k(35)=5E-7; k_(35)=1E-4; kc(35)=1;

%Casp12 + XIAP <--> Casp12:XIAP
k(36)=2E-7; k_(36)=1E-4;

% calpain* + Bcl2c <--> calpain*:Bcl2c --> cBcl2c
k(37)=2E-7; k_(37)=1E-4; kc(37)=.1;

%% Numerical Integration

%set timespan
tt=[0:60:tf_input]; %Calculate every 60 seconds

%set initial conditions
x0=IC;

%run ODE
[t, x]=ode15s(@rhs,tt,x0);

%% Subfunctions

    function xp=rhs(t,x)
        % ODE with dimensions
        xp=double(x);
       xp(1) = -k(1)*x(1)*x(2) +k_(1)*x(3) ; % Ligand 

xp(2) = -k(1)*x(1)*x(2) +k_(1)*x(3) ; % R

xp(3) =  k(1)*x(1)*x(2) -k_(1)*x(3) -kc(1)*x(3) ; % L:R complex

xp(4) = kc(1)*x(3) +...
    -k(2)*x(4)*x(5) +k_(2)*x(6) +...
    -k(3)*x(4)*x(7) +k_(3)*x(8) +kc(3)*x(8); % R*

xp(5)= -k(2)*x(4)*x(5) +k_(2)*x(6) ; % flip

xp(6)=  k(2)*x(4)*x(5) -k_(2)*x(6) ; % flip:R*

xp(7) = -k(3)*x(4)*x(7) +k_(3)*x(8) +...
    -k(7)*x(7)*x(17) +k_(7)*x(18) ; % pC8 

xp(8) =  k(3)*x(4)*x(7) -k_(3)*x(8) -kc(3)*x(8) ; % R*:pC8

xp(9) =  kc(3)*x(8) +...
    -k(4)*x(9)*x(10) +k_(4)*x(11) +...
    -k(5)*x(9)*x(12) +k_(5)*x(13) +kc(5)*x(13) +...
    +kc(7)*x(18) +...
    -k(10)*x(9)*x(24) +k_(10)*x(25) +kc(10)*x(25) ; % C8

xp(10) = -k(4)*x(9)*x(10) +k_(4)*x(11) ; % Bar

xp(11) =  k(4)*x(9)*x(10) -k_(4)*x(11) ; % Bar:C8

xp(12)= -k(5)*x(9)*x(12) +k_(5)*x(13) +...
    -k(25)*x(12)*x(53) +k_(25)*x(54) ; % pC3

xp(13)=  k(5)*x(9)*x(12) -k_(5)*x(13) -kc(5)*x(13) ; % C8:pC3

xp(14)=  kc(5)*x(13) +...
    -k(6)*x(14)*x(15) +k_(6)*x(16) +kc(6)*x(16) +...
    -k(8)*x(14)*x(19) +k_(8)*x(20) +...
    -k(9)*x(14)*x(21) +k_(9)*x(22) +kc(9)*x(22) +...
    +kc(25)*x(54)...
    -k(32)*x(65)*x(14)+k_(32)*x(66)+kc(32)*x(66); % C3 

xp(15)= -k(6)*x(14)*x(15) +k_(6)*x(16) ; % pC6

xp(16)=  k(6)*x(14)*x(15) -k_(6)*x(16) -kc(6)*x(16) ; % C3:pC6

xp(17)=  kc(6)*x(16) +...
    -k(7)*x(7)*x(17) +k_(7)*x(18) +kc(7)*x(18) ; % C6

xp(18)=  k(7)*x(7)*x(17) -k_(7)*x(18) -kc(7)*x(18) ; % C6:pC8

xp(19)= -k(8)*x(14)*x(19) +k_(8)*x(20) +kc(8)*x(20) +... 
    -k(27)*x(19)*x(53) +k_(27)*x(56) +...
    -k(28)*x(19)*x(55) +k_(28)*x(57) ...
    -k(36)*x(70)*x(19)+k_(36)*x(73); % XIAP

xp(20)=  k(8)*x(14)*x(19) -k_(8)*x(20) -kc(8)*x(20) ; % XIAP:C3

xp(21)= -k(9)*x(14)*x(21) +k_(9)*x(22) ; % PARP

xp(22)=  k(9)*x(14)*x(21) -k_(9)*x(22) -kc(9)*x(22) ; % C3:PARP

xp(23)= kc(9)*x(22) ; % CPARP

xp(24)= -k(10)*x(9)*x(24) +k_(10)*x(25)....
            -k(35)*x(24)*x(67)+k_(35)*x(72); % Bid

xp(25)=  k(10)*x(9)*x(24) -k_(10)*x(25) -kc(10)*x(25) ; % C8:Bid

xp(26)=  kc(10)*x(25) +...
    -k(11)*x(26)*x(27) +k_(11)*x(28) +...
    -k(12)*x(26)*x(29) +k_(12)*x(30) + kc(12)*x(30)...
    +kc(34)*x(69)*x(9)+kc(35)*x(72); % tBid

xp(27)= -k(11)*x(26)*x(27) +k_(11)*x(28)...
        -k(37)*x(67)*x(27)+k_(37)*x(74); % Bcl2c

xp(28)= +k(11)*x(26)*x(27) -k_(11)*x(28) ; % Bcl2c:tBid

xp(29)= -k(12)*x(26)*x(29) +k_(12)*x(30) ; % Bax

xp(30)=  k(12)*x(26)*x(29) -k_(12)*x(30) - kc(12)*x(30) ; % tBid:Bax

xp(31)=  kc(12)*x(30) +...
    -k(13)*x(31) + k_(13)*x(32) ; % Bax*

xp(32)=  k(13)*x(31) - k_(13)*x(32) +...
    -1/v*k(14)*x(32)*x(33) +k_(14)*x(34) +...
    -1/v*2*k(15)*x(32)^2 +2*k_(15)*x(35) ; % Baxm

xp(33)= -1/v*k(14)*x(32)*x(33) +k_(14)*x(34) +...
    -1/v*k(16)*x(33)*x(35) +k_(16)*x(36) +...
    -1/v*k(18)*x(33)*x(37) +k_(18)*x(38) ; % Bcl2 

xp(34)=  1/v*k(14)*x(32)*x(33) -k_(14)*x(34) ; % Baxm:Bcl2

xp(35)=  1/v*k(15)*x(32)^2 -k_(15)*x(35) +...
    -1/v*k(16)*x(33)*x(35) +k_(16)*x(36) +...
    -2/v*k(17)*x(35)^2 +2*k_(17)*x(37) ; % Bax2

xp(36)=  1/v*k(16)*x(33)*x(35) -k_(16)*x(36) ; % Bax2:Bcl2

xp(37)= 1/v*k(17)*x(35)^2 -k_(17)*x(37)+... 
    -1/v*k(18)*x(33)*x(37) +k_(18)*x(38) +...
    -1/v*k(19)*x(39)*x(37) +k_(19)*x(40) ; % Bax4 

xp(38)= 1/v*k(18)*x(33)*x(37) -k_(18)*x(38) ; % Bax4:Bcl2

xp(39)= -1/v*k(19)*x(39)*x(37) +k_(19)*x(40); % M

xp(40)=  1/v*k(19)*x(39)*x(37) -k_(19)*x(40) -kc(19)*x(40) ; % Bax4:M

xp(41)=  kc(19)*x(40) +...
    -1/v*k(20)*x(41)*x(42) +k_(20)*x(43) +kc(20)*x(43) +...
    -1/v*k(21)*x(41)*x(45) +k_(21)*x(46) +kc(21)*x(46) ; % M*

xp(42)= -1/v*k(20)*x(41)*x(42) +k_(20)*x(43) ; % CytoCm

xp(43)=  1/v*k(20)*x(41)*x(42) -k_(20)*x(43) -kc(20)*x(43) ; % M*:CytoCm

xp(44)=  kc(20)*x(43) +...
    -k(22)*x(44) +k_(22)*x(48) ; % CytoCr

xp(45)= -1/v*k(21)*x(41)*x(45) +k_(21)*x(46) ; % Smacm

xp(46)=  1/v*k(21)*x(41)*x(45) -k_(21)*x(46) -kc(21)*x(46) ; % M*:Smacm

xp(47)=  kc(21)*x(46) +...
    -k(26)*x(47) +k_(26)*x(55) ; % Smacr

xp(48)=  k(22)*x(44) -k_(22)*x(48) +...
    -k(23)*x(48)*x(49) +k_(23)*x(50) +kc(23)*x(50) ; % CytoC

xp(49)= -k(23)*x(48)*x(49) +k_(23)*x(50) ; % Apaf

xp(50)=  k(23)*x(48)*x(49) -k_(23)*x(50) -kc(23)*x(50) ; % Apaf:CytoC

xp(51)= kc(23)*x(50) +...
    -k(24)*x(51)*x(52) +k_(24)*x(53); % Apaf*

xp(52)= -k(24)*x(51)*x(52) +k_(24)*x(53) -k(33)*x(52)*x(68)...
            -k(34)*x(70)*x(52)+k_(34)*x(71); % pC9

xp(53)=  k(24)*x(51)*x(52) -k_(24)*x(53) +...
    -k(25)*x(12)*x(53) +k_(25)*x(54) +kc(25)*x(54) +...
    -k(27)*x(19)*x(53) +k_(27)*x(56) +kc(34)*x(71) ; % Apop

xp(54)=  k(25)*x(12)*x(53) -k_(25)*x(54) -kc(25)*x(54) ; % Apop:pC3

xp(55)=  k(26)*x(47) -k_(26)*x(55) +...
    -k(28)*x(19)*x(55) +k_(28)*x(57) ; % Smac

xp(56)=  k(27)*x(19)*x(53) -k_(27)*x(56) ; % Apop:XIAP 

xp(57)=  k(28)*x(19)*x(55) -k_(28)*x(57) ; % Smac:XIAP

xp(58)=  kc(8)*x(20); % C3_Ub

xp(59)=  -k(29)*x(59)*x(60)+k_(29)*x(61)...
        -k(30)*x(59)*x(62)+k_(29)*x(63); %calcium

xp(60)=  -k(29)*x(59)*x(60)+k_(29)*x(61); %calpain

xp(61)=  k(29)*x(59)*x(60)-k_(29)*x(61)...
            -k(31)*x(61)*x(63)+k_(31)*x(64);%calpainC
    
xp(62)=  -k(30)*x(59)*x(62)+k_(30)*x(63);  %calpastatin

xp(63)=  k(30)*x(59)*x(62)-k_(30)*x(63)...
            -k(31)*x(61)*x(63)+k_(31)*x(64); %calpastatin*

xp(64)=   k(31)*x(61)*x(63)-k_(31)*x(64)-kc(31)*x(64); %calpastatin*:calpainC

xp(65)=   kc(31)*x(64)-k(32)*x(65)*x(14)+k_(32)*x(66);  %calpain_blocked

xp(66)=   k(32)*x(65)*x(14)-k_(32)*x(66)-kc(32)*x(66);  %C3:calpain_blocked

xp(67)=   kc(32)*x(66)...
            -k(33)*x(67)*x(68)+k_(33)*x(69)+kc(33)*x(69)...
            -k(35)*x(24)*x(67)+k_(35)*x(72)+kc(35)*x(72)...
            -k(37)*x(67)*x(27)+k_(37)*x(74)+kc(37)*x(74);  %calpain*
    
xp(68)=   -k(33)*x(67)*x(68)+k_(33)*x(69);   %pC12

xp(69)=   k(33)*x(67)*x(68)-k_(33)*x(69)-kc(33)*x(69);   %pC12:calpain*
        
xp(70)=   kc(33)*x(69)...
            -k(34)*x(70)*x(52)+k_(34)*x(71)+kc(34)*x(71)...
            -k(36)*x(70)*x(19)+k_(36)*x(73); %C12
       
xp(71)=   k(34)*x(70)*x(52)-k_(34)*x(71)-kc(34)*x(71); %C12:pC9
        
xp(72)=   k(35)*x(24)*x(67)-k_(35)*x(72)-kc(35)*x(72);  %calpain*:Bid

xp(73)=   k(36)*x(70)*x(19)-k_(36)*x(73);  %C12:XIAP

xp(74)=   k(37)*x(67)*x(27)-k_(37)*x(74)-kc(37)*x(74); %calpain*:Bcl2c

xp(75)=   kc(37)*x(74);     %cleaved Bcl2c
    end
end