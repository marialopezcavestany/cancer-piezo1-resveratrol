close all
clear

%TRAIL exposure duration [s]
dur = 24*3600;

%End time [s]
endTime = 24*3600;
timePoint = endTime;

%TRAIL concentration (molecules/cell)
Lconc = 3000;     % 50ng/mL TRAIL

%Receptor concentration (molecules/cell)
Rconc = 200;

%calcium concentration in uM -> choose one
%calcium = 0.116; %Calcium free buffer
%calcium = 0.162; %Control
%calcium = 0.134; %RSV
calcium = 1; %Yoda1
%calcium = 1.760; %Yoda1 + RSV

%Set up the random normal Bcl-2 and XIAP calculation
n=500;
randBCL2 = 2E6 + randn (1,n)*1E6;
randXIAP = 1E5 + randn (1,n)*1E5;
storeValues = [];
storeCPARP = [];
countCells = 0;
countDeadCells = 0;

%This loop runs both the Duration and InitializeProblem scripts for each of
%the randomly generated Bcl-2 and XIAP values
for i=1:n
    
    cBCL2 = randBCL2(1,i);
    newXIAP = randXIAP (1,i);
    
    %This helps keep track of where you are in the loop
    %disp(i);

    if cBCL2 > 0.01
         
        if newXIAP > 0.01
        % Run the simulation using the random Bcl-2 and XIAP concentrations
         [xnew,tnew] = JG_Duration(Lconc,Rconc,dur,calcium,cBCL2,newXIAP);
    
        %Set up a matrix that records all the model outputs and also just
        %the cPARP values to make graphing at the end simpler
            storeValues = [storeValues xnew];
            storeCPARP = [storeCPARP xnew(:,23)];
        
            countCells = countCells + 1;

        %This adds up the number of dead cells (number of cells where the 
        % cPARP concentration is above the apoptosis threshold)
        if xnew(1441,23) >= 5E5
            countDeadCells = countDeadCells + 1;
            
        end
      end
    end

end

%Calculate the cell viability based on the number of cells below the
%apoptosis threshold
CellViability = ((countCells - countDeadCells) / countCells) * 100;

%Add cell apoptosis line and make the figure
Apoptosis = ones(1441,1)* 5E5;
graphCPARP = [storeCPARP Apoptosis];

%Calculate random variable statistics and histograms
 stats.meanBCL = mean(randBCL2);
 stats.medianBCL = median(randBCL2);
 stats.stdBCL = std(randBCL2);
 stats.meanXIAP = mean(randXIAP);
 stats.medianXIAP = median(randXIAP);
 stats.stdXIAP = std(randXIAP);
 %figure;
 %histogram(randBCL2);
 %figure;
 %histogram(randXIAP);

 %Once this script has finished running, run the FigureCustomization script
 %to make the figure
 