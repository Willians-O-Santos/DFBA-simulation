% dFBA - function that carries out the DFBA simulation
% 
% This function is called by the simulator.m file in order to performe the DFBA
% simulation for each iteration
%  
% syntax:
%    dydt = dFBA(t,y)
% 
% Inputs:
%    input1 - t, time at each iteration
%    input2 - y, concentration vector
% 
% Outputs:
%    dydt = the system of ordinary differential equations to be solved
% 
% 
% Other m-files required: simulator.m, FBA.m
%
% Author: Rafael David and Willians Santos
% Affiliation: University of São Paulo (USP)
% email: willians_soad@gmail.com
% Created: 03-Mar-2023 ; Last revision: 03-Mar-2023 

% ------------------------------------------------------------------------ 


function dydt = dFBA(t,y)

% ------------------------------------------------------------------------
% Maximum oxygen uptake:

% vo2 = 0.0; % For anaerobic simulations
vo2 = 15.0; % E. coli under aerobic conditions
% vo2 = 5.0; % C. necator under aerobic conditions
% vo2 = 1.5; % S. cerevisiae under aerobic conditions
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% Glucose uptake kinetics:

if y(2) <=0.0 %stops consuming carbon source if concentration falls to zero
  vg=0.0;
else
    
% Uncomment the uptake expression for the model and condition of your
% simulation:

vg = (10.5*y(2)/(0.015+y(2)+((y(2)^2)/2140.12)))*(1/(1 + (y(5)/100)))*(1/(1 + (y(7)/434.12))); % E. coli glucose under aerobic conditions
% vg = (18.5*y(2)/(0.015+y(2)+((y(2)^2)/2140.12)))*(1/(1 + (y(5)/100)))*(1/(1 + (y(7)/434.12))); % E. coli glucose under anaerobic conditions
% vg = (13*y(2)/(0.019+y(2)+((y(2)^2)/11370.6)))*(1/(1 + (y(5)/100)))*(1/(1 + (y(7)/434.12))); % E. coli glycerol under aerobic conditions
% vg = (7.5*y(2)/(0.11+y(2)+((y(2)^2)/24042.6)))*(1/(1 + (y(5)/100)))*(1/(1 + (y(7)/434.12))); % E. coli xylose under aerobic conditions


% vg = (3.0*y(2)/(3.36+y(2)+((y(2)^2)/9738.47))); % C. necator glucose under aerobic conditions, and does not produce ac and eth
% vg = (8*y(2)/(11.04+y(2)+((y(2)^2)/47237.1))); % C. necator glycerol under aerobic conditions, and does not produce ac and eth

    
% vg = (22.5*y(2)/(4.884+y(2)+((y(2)^2)/27102.575)))*(1/(1 + (y(5)/107.143)))*(1/(1 + (y(7)/15.35))); % S. cerevisiae glucose under aerobic conditions
% vg = (32*y(2)/(80.33+y(2)+((y(2)^2)/32523.81)))*(1/(1 + (y(5)/107.143)))*(1/(1 + (y(7)/15.35))); % S. cerevisiae xylose under aerobic conditions


end
% ------------------------------------------------------------------------


 
% ------------------------------------------------------------------------
% calling the FBA.m function.
% If the user removes or adds metabolites to be tracked, 
% the output of the function has to be changed to account for the new
% metabolites.

[mi,vam,vphb,vac,vlac,veth,vform,vsucc,exitflag] = FBA(vg,vo2);
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
%avoids negative concentrations:

if exitflag == -2
 vg=0.0;
 vo2=0.0;

else
end
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% system of ordinary differential equations to be solved. If the user wants
% to remove or add other metabolites to be tracked, the differential
% equation of this metabolites needs to be added.

dydt(1,1) = mi*y(1); % biomass concentration as a function of time
dydt(2,1) = -vg*y(1); % carbon source concentration as a function of time
dydt(3,1) = vam*y(1); % ammonia concentration as a function of time
dydt(4,1) = vphb*y(1); % phb concentration as a function of time
dydt(5,1) = vac*y(1); % acetate concentration as a function of time
dydt(6,1) = vlac*y(1); % lactate concentration as a function of time
dydt(7,1) = veth*y(1); % ethanol concentration as a function of time
dydt(8,1) = vform*y(1); % formate concentration as a function of time
dydt(9,1) = vsucc*y(1); % succinate concentration as a function of time
% ------------------------------------------------------------------------

end