% simulator - Main file for DFBA batch simulations
% 
% This is the file that needs to be run in order to performe the DFBA
% simulations.
%  
%     
% 
% Inputs:
%    input1 - y0, the initial concentration and conditions vector
%    input2 - tspan, the time span of the batch process
% 
% Outputs:
%    output1 - y, time profile of all the specified metabolites 
%    output2 - t, vector containing the duration of the simulated batch process
% 
% Other m-files required: FBA.m, DFBA.m
%
% Author: Rafael David and Willians Santos
% Affiliation: University of São Paulo (USP)
% email: willians_soad@gmail.com
% Created: 03-Mar-2023 ; Last revision: 03-Mar-2023 

% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% User inputs:

% Initial concentrations vector y0. Enter the values for biomass (g/L), glucose or glycerol or xylose (mmol/L), 
% ammonia (mmol/L), PHB (mmol/L), acetate (mmol/L), lactate (mmol/L), ethanol (mmol/L), 
% formate (mmol/L), and succinate (mmol/L), respectively. Metabolites can be removed or other
% metabolites can be added, but that will require changes in other parts of
% the scripts as well, as will be pointed out in each specific section.

% For the present work, an initial molar concentration of 138.76 mmol/L was used
% when glucose was the carbon source, an initial molar concentration of 262.96 mmol/L was used
% when glycerol was the carbon source, and an initial molar concentration of 166.52 mmol/L was used
% when xylose was the carbon source. That is so that the same inital concentration of
% carbon source in terms of mass (25 g/L), was used for all simulations.
% Initial biomass concentration used was 0.25 g/L.
% initial ammonia concentration used was 166.67 mmol/L.


y0=[0.25 138.76 166.67 0.0 0.0 0.0 0.0 0.0 0.0];  

% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% the time span for the batch simulation, tspan. Indicates how long will
% the batch operation last and the time step used. For this study, a time
% step of 0.05 was used. As for the batch duration, the times chosen for
% each simulation were such as to allow all glucose to be consumed.

h = 0.05;
tspan=[0.0:h:6];

% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% Setting the approach.

% Direct approach (DA) for DFBA simulation. This is the approach used in this study.

[t,y] = ode15s(@(t,y)dFBA(t,y),tspan,y0); % ODE15s to solve system of differential equations
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% plotting the graph for the metabolites as a function of time.
% If the user removes or adds metabolites, the plot function has to be
% changed to account for that.

figure(1)
subplot(2,1,1)
plot(t,y(:,2))
hold on
plot(t,y(:,3))
hold on
plot(t,y(:,4))
hold on
plot(t,y(:,5))
hold on
plot(t,y(:,6))
hold on
plot(t,y(:,7))
hold on
plot(t,y(:,8))
hold on
plot(t,y(:,9))
xlabel('Time [h]')
ylabel('Concentrations [mmol/L]')
legend('Carbon source','Ammonia', 'PHB', 'Acetate', 'Lactate', 'Ethanol', 'Formate', 'Succinate', 'Location', 'eastoutside','FontSize',8)


%plotting the graph for the biomass as a function of time.

subplot(2,1,2)
plot(t,y(:,1))
xlabel('Time [h]')
ylabel('Biomass [g/L]')
% ------------------------------------------------------------------------
