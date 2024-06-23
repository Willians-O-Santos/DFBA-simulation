tic
% simulator - Main file for DFBA fed-batch simulations
% 
% This is the file that needs to be run in order to performe the DFBA
% simulations.
%  
%     
% 
% Inputs:
%    input1 - y0, the initial concentration and conditions vector
%    input2 - tend, final time of the simulation
%    input3 - ne, the number of elements to be used in the simulation
%    input4 - teta.vg_max, maximum glucose uptake rate
%    input5 - teta.Kg, the glucose saturation constant
%    input6 - teta.Kig, the glucose inhibition constant
%    input7 - teta.Sfeed, the feed glucose concentration
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
% formate (mmol/L), and succinate (mmol/L), and initial volume (L), respectively. 
% Metabolites can be removed or other metabolites can be added, but that will require 
% changes in other parts of the scripts as well, as will be pointed out in each specific section.

% For the present work, an initial molar concentration of 83.26 mmol/L was used
% when glucose was the carbon source, an initial molar concentration of 157.78 mmol/L was used
% when glycerol was the carbon source, and an initial molar concentration of 99.91 mmol/L was used
% when xylose was the carbon source. That is so that the same inital concentration of
% carbon source in terms of mass (15 g/L), was used for all simulations.
% Initial biomass concentration used was 0.25 g/L.
% initial ammonia concentration used was 550 mmol/L.

y0=[0.25 83.26 550 0.0 0.0 0.0 0.0 0.0 0.0 200000];
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% time step:

tend=8; % Enter final time of the simulation, (h)

ne=160; % Enter a number of elements, which in turn defines the number of steps

h=tend/ne; % number of steps
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% Kinetic parameters, just uncomment the ones for the model that will be
% used in the simulations.


% Maximum carbon source uptake rate:

teta.vg_max=10.5; % E. coli glucose under aerobic conditions,(mmol/gCDW.h)
% teta.vg_max=18.5; % E. coli glucose under anaerobic conditions,(mmol/gCDW.h)
% teta.vg_max=13; % E. coli glycerol under aerobic conditions,(mmol/gCDW.h)
% teta.vg_max=7.5; % E. coli xylose under aerobic conditions,(mmol/gCDW.h)



% teta.vg_max=3.0; % C. necator glucose under aerobic conditions,(mmol/gCDW.h)
% teta.vg_max=8; % C. necator glycerol under aerobic conditions,(mmol/gCDW.h)


% teta.vg_max=22.5; % S. cerevisiae glucose under aerobic conditions,(mmol/gCDW.h)
% teta.vg_max=32; % S. cerevisiae xylose under aerobic conditions,(mmol/gCDW.h)


% saturation constant:

teta.Kg=0.015; % E. coli glucose saturation constant,(mmol/L)
% teta.Kg=0.019; % E. coli glycerol saturation constant,(mmol/L)
% teta.Kg=0.11; % E. coli xylose saturation constant,(mmol/L)


% teta.Kg=3.36; % C. necator glucose saturation constant,(mmol/L)
% teta.Kg=11.04; % C. necator glycerol saturation constant,(mmol/L)


% teta.Kg=4.884; % S. cerevisiae glucose saturation constant,(mmol/L)
% teta.Kg=80.33; % S. cerevisiae xylose saturation constant,(mmol/L)


% carbon source inhibition constant:

teta.Kig=2140.12; % E. coli glucose,(mmol/L)
% teta.Kig=11370.6; % E. coli glycerol,(mmol/L)
% teta.Kig=24042.6; % E. coli xylose,(mmol/L)


% teta.Kig=9738.47; % C. necator glucose,(mmol/L)
% teta.Kig=47237.1; % C. necator glycerol,(mmol/L)


% teta.Kig=27102.575; % S. cerevisiae glucose,(mmol/L)
% teta.Kig=32523.81; % S. cerevisiae xylose,(mmol/L)


% acetate inhibition constant:

teta.Kiac=100; % E. coli,(mmol/L)


% teta.Kiac= C. necator does not produce acetate; 


% teta.Kiac=107.143; % S. cerevisiae,(mmol/L)


% ethanol inhibition constant:

teta.Kieth=434.12; % E. coli,(mmol/L)


% teta.Kieth= C. necator does not produce ethanol;


% teta.Kieth=15.35; % S. cerevisiae,(mmol/L)

% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% feed concentration:

teta.Sfeed=2775.31; % Feed glucose concentration used in this work (mmol/L)

%teta.Sfeed=5259.17; % Feed glycerol concentration used in this work (mmol/L)

%teta.Sfeed=3330.45; % Feed xylose concentration used in this work (mmol/L)

% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% Preallocating variables.

tspan=[0.0 h]; % time span of the simulations
y = zeros(ne+1,10); % preallocating concentrations matrix to improve performance
y(1,:) = y0; % Adding the initial conditions in the concentration matrix
t = zeros(ne+1,1); % preallocating time vector to improve performance, (h)
Fin = zeros(ne+1,1); % preallocating feed vector to improve performance, (L/h)
k=1;
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% Setting kinetic profiles for the first iteration.

% Glucose uptake kinetic, (mmol/gCDW.h):

vg = (teta.vg_max*y0(2)/(teta.Kg + y0(2) + ((y0(2)^2)/teta.Kig)))*(1/(1 + (y0(5)/teta.Kiac)))*(1/(1 + (y0(7)/teta.Kieth))); % for E. coli and S. cerevisiae
% vg = (teta.vg_max*y0(2)/(teta.Kg + y0(2) + ((y0(2)^2)/teta.Kig))); % for C. necator

% Feed expression for exponential feeding profile, (L/h). If the user
% intends to use a different feeding profile in the simulation, the feeding
% expressions have to be changed.

Fin(1)=(vg*y0(1)*y0(10)) / (teta.Sfeed-y0(2)); 



% For no feeding, change the feeding expressions to zero.
% Fin(1)= 0; % No feed
% ------------------------------------------------------------------------


for i=0:h:tend-h

% ------------------------------------------------------------------------
%Direct approach (DA)
    [ta,ya] = ode15s(@(t,y)dFBA(t,y,Fin(k),teta),tspan,y0); % Calling DFBA function
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% Updating variables:
    y0=ya(end,:); % Updating initial conditions for next iteration
    y(k+1,:)=ya(end,:); % Updating concentrations matrix
    t(k+1)=h*k; % Updating time vector
    
    vg = (teta.vg_max*y0(2)/(teta.Kg + y0(2) + ((y0(2)^2)/teta.Kig)))*(1/(1 + (y0(5)/teta.Kiac)))*(1/(1 + (y0(7)/teta.Kieth))); % Updating uptake, for E. coli and S. cerevisiae
    %vg = (teta.vg_max*y0(2)/(teta.Kg + y0(2) + ((y0(2)^2)/teta.Kig))); % Updating uptake,for C. necator
    Fin(k+1)=(vg*y0(1)*y0(10)) / (teta.Sfeed-y0(2)); % Updating feed
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% Stops feeding after reaching a defined maximum volume.
% If user does not intend to do so, just comment this block of code.
if y0(10) < 204000 % Where 235000 L is the defined maximum volume in this study
    Fin(k+1)= (vg*y0(1)*y0(10)) / (teta.Sfeed-y0(2));
else
    Fin(k+1)= 0;
end
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% Uncomment if the intention is to not have feeding.
%Fin(k+1)= 0; %
% ------------------------------------------------------------------------


k=k+1;


    
end

% ------------------------------------------------------------------------
% plotting the graph for the metabolites as a function of time.
% If the user removes or adds metabolites, the plot function has to be
% changed to account for that.

figure(2)
subplot(2,2,1)
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

subplot(2,2,2)
plot(t,y(:,1))
xlabel('Time [h]')
ylabel('Biomass [g/L]')



%plotting the graph for the volume as a function of time.

subplot(2,2,3)
plot(t,y(:,10))
xlabel('Time [h]')
ylabel('Volume [L]')



%plotting the graph for the feeding as a function of time.

subplot(2,2,4)
plot(t,Fin)
xlabel('Time [h]')
ylabel('F [L/h]')

toc