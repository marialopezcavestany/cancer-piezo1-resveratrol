close all
clear

%TRAIL exposure duration [s]
dur = 24*3600;

%End time [s]
endTime = 24*3600;
timePoint = endTime;

%TRAIL concentration (molecules/cell) -> Choose one
%Lconc = 0;        % No TRAIL treatment
%Lconc = 30;       % 0.05ng/mL TRAIL
Lconc = 3000;     % 50ng/mL TRAIL
%Lconc = 12000;    % 200ng/mL TRAIL

%Receptor concentration (molecules/cell)
Rconc = 200;

%Cytosolic Bcl-2 concentration
cBCL2 = 2E6;  % cytosolic Bcl-2
newXIAP = 1E5; %XIAP

%Calcium concentration in uM -> choose one
%calcium = 0.116;      %Calcium free buffer
%calcium = 0.162;      %Control
%calcium = 0.134;      %RSV
calcium = 1;          %Yoda1
%calcium = 1.760;      %Yoda1 + RSV

%Runs the ODEs
[xnew,tnew] = JG_Duration(Lconc,Rconc,dur,calcium,cBCL2,newXIAP);

%Plot cPARP
figure;
plot(tnew, xnew(:,23));

%Plot cytosolic Smac
figure;
plot(tnew, xnew(:,55));
