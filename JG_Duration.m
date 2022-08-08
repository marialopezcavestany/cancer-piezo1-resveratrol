function [xnew,tnew] = JG_Duration (Lconc,Rconc,dur,calcium,cBCL2,newXIAP)

% this function is used to simulate the continuous duration of TRAIL exposure with death
% receptor for 3D binding kinetics

% Lconc is the concentration of ligand binding receptor [#/cc]
% Rconc is the concentration of receptor binding ligand [#/cc]

% xnew is the matrix of reagents (each collumn is a different reagent, and
% each row is a given time.
% tnew is the vector of time corresponding to xnew.

xnew = [];
tnew = [];

[t_col,x_col] = JG_InitializeProblem(Lconc,Rconc,dur,calcium,cBCL2,newXIAP);

xnew = [xnew;x_col];
tnew = [tnew;t_col];

clear t_col x_col

end