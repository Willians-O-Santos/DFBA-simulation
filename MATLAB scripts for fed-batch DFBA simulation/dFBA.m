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
%    input3 - Fin, feeding uptake rate (L/h)
%    input4 - teta, kinetic parameters
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


function dydt = dFBA(t,y,Fin,teta)

% ------------------------------------------------------------------------
% Maximum oxygen uptake:

% vo2 = 0.0; % For anaerobic simulations
vo2 = 15.0; % E. coli under aerobic conditions
% vo2 = 5.0; % C. necator under aerobic conditions
% vo2 = 1.5; % S. cerevisiae under aerobic conditions
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% Carbon source uptake kinetics:

vg = (teta.vg_max*y(2)/(teta.Kg + y(2) + ((y(2)^2)/teta.Kig)))*(1/(1 + (y(5)/teta.Kiac)))*(1/(1 + (y(7)/teta.Kieth))); % For E. coli and S. cerevisiae
% vg = (teta.vg_max*y(2)/(teta.Kg + y(2) + ((y(2)^2)/teta.Kig))); %For C. necator

% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% Feed concentration:

Gf=teta.Sfeed; %mmol/L
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
 vo=0.0;

else
end
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% system of ordinary differential equations to be solved. If the user wants
% to remove or add other metabolites to be tracked, the differential
% equation of this metabolites needs to be added.

dydt(1,1) = (mi*y(1) - ((Fin*y(1))/y(10))); % biomass concentration as a function of time
dydt(2,1) = ((Fin/y(10))*(Gf - y(2)) - vg*y(1)); % carbon source concentration as a function of time
dydt(3,1) = (vam*y(1) - ((Fin*y(3))/y(10))); % ammonia concentration as a function of time
dydt(4,1) = (vphb*y(1) - ((Fin*y(4))/y(10))); % phb concentration as a function of time
dydt(5,1) = (vac*y(1) - ((Fin*y(5))/y(10))); % acetate concentration as a function of time
dydt(6,1) = (vlac*y(1) - ((Fin*y(6))/y(10))); % lactate concentration as a function of time
dydt(7,1) = (veth*y(1) - ((Fin*y(7))/y(10))); % ethanol concentration as a function of time
dydt(8,1) = (vform*y(1) - ((Fin*y(8))/y(10))); % formate concentration as a function of time
dydt(9,1) = (vsucc*y(1) - ((Fin*y(9))/y(10))); % succinate concentration as a function of time
dydt(10,1) = Fin; % volume as a function of time

end