% FBA - function that carries out the FBA simulation
% 
% This function is called by the DFBA.m file in order to performe the FBA
% simulation for each iteration.
%  
% syntax:
%    [mi,vam,vphb,vac,vlac,veth,vfor,vsucc,exitflag] = FBA(vg,vo2)
% 
% Inputs:
%    input1 - vg, the glucose uptake rate
%    input2 - vo2, the oxygen uptake rate
% 
% Outputs:
%    output1 - mi, the resulting cell growth rate of the current iteration 
%    output2 - vam, the resulting ammonia uptake flux of the current iteration
%    output3 - vphb, the resulting PHB production flux of the current iteration
%    output4 - vac, the resulting acetate production flux of the current iteration
%    output5 - veth, the resulting ethanol production flux of the current iteration
%    output6 - vform, the resulting formate production flux of the current iteration
%    output7 - vac, the resulting acetate production flux of the current iteration
%    output8 - vsucc, the resulting succinate production flux of the current iteration
%    output9 - exitflag, returns a value that describes the exit condition
% 
% 
% Other m-files required: simulator.m, DFBA.m
%
% Author: Rafael David and Willians Santos
% Affiliation: University of São Paulo (USP)
% email: willians_soad@gmail.com
% Created: 03-Mar-2023 ; Last revision: 03-Mar-2023 

% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% Main FBA function. If the user removes or adds metabolites to be tracked, 
% the output of the function has to be changed to account for the new
% metabolites.

function [mi,vam,vphb,vac,vlac,veth,vform,vsucc,exitflag] = FBA(vg,vo2)

% -----------------------------------------------------------------------



% ------------------------------------------------------------------------
% Uncomment the Indexes for the model you want to use:

% E. coli iML1515 model indexes of key reactions:

glc = 181;
gly = 100; %glycerol exchange reaction
xyl = 212; %xylose exchange reaction
o2 = 1982;
ac = 107;
co2 = 46;
am = 1501; %nh4 ammonia
obj = 2669; %biomass
phb = 2717;
eth = 1526; 
form = 217; 
succ = 113; 
lac = 348;  



% C. necator H16 RehMBEL1391 model indexes of key reactions: 

% glc = 1541;
% gly = 110; %glycerol exchange reaction
% o2 = 143;
% ac = 10;
% co2 = 133;
% am = 125; % NH4 ammonia
% obj = 21; %biomass
% phb = 60;
% lac = 152;
% succ = 70;
% eth = 162;
% form= 119;



% S. cerevisiae iMM904 model indexes of key reactions: 

% glc = 508;
% xyl = 638; %xylose exchange reaction
% o2 = 599;
% ac = 502;
% co2 = 547;
% am = 596; %NH4 ammonia
% obj = 1521; %biomass
% phb = 1581;
% lac = 566;
% succ = 594;
% eth = 473;
% form= 487;
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------

% Reading the model, uncomment the block related to the model you want to use:


% E. coli iML1515 model:

model = load('iML1515_PHB_NADPH.mat'); % with NADPH-dependent PHB pathway
% model = load('iML1515_PHB_NADH.mat'); % with NADH-dependent PHB pathway
model= model.model;



% C. necator H16 RehMBEL1391 model:

% c_necator = load('Cupriavidus_necator_H16_bounds30_with_glucose_no_EXac_EXeth_EXlac.mat');
% model= c_necator.model;



% S. cerevisiae iMM904 model:

% model = load('iMM904_with_NADPH_PHB.mat'); 
% model = load('iMM904_with_NADH_PHB.mat'); 
% model= model.model;
% ------------------------------------------------------------------------




% ------------------------------------------------------------------------
% Setting up the linprog function:

objvector = model.c; % variable with the objective vector from the model
f=-objvector; % sets the objective for the linprog function
A=[]; % stoichiometric matrix for inequality constraints, null in the case of FBA
b=[]; % right hand side of inequality constraints, null in the case of FBA
Aeq=model.S; % stoichiometric matrix of the model
beq=model.b; % right hand side of equality constraints
lb=model.lb; % variable with the lower bounds from the model
ub=model.ub; % variable with the upper bounds from the model
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% Constraints for the linprog function:

%For E. coli:
% set glucose as the carbon source
lb(glc)=-vg;
lb(gly)=0;
lb(xyl)=0; 

% set glycerol as the carbon source
% lb(glc)=0;
% lb(gly)=-vg;
% lb(xyl)=0; 

% set xylose as the carbon source
% lb(glc)=0;
% lb(gly)=0;
% lb(xyl)=-vg;

% set oxygen uptake
% lb(o2)=-vo2;

% set PHB flux
% Fix a chosen flux to PHB synthesis, by setting the upper and lower bound, in order to explore the trade-off between biomass and PHB formation.
% ub(phb)=0;
% lb(phb)=0; 





%For C. necator:

% set glucose as the carbon source
% lb(glc)=-vg;
% lb(gly)=0;
 
% set glycerol as the carbon source
% lb(glc)=0;
% lb(gly)=-vg;

% set oxygen uptake
%  lb(o2)=-vo2;

% set PHB flux
% Fix a chosen flux to PHB synthesis, by setting the upper and lower bound, in order to explore the trade-off between biomass and PHB formation.
% ub(phb)=0;
% lb(phb)=0;




%For S. cerevisiae:

% set glucose as the carbon source
% lb(glc)=-vg;
% lb(xyl)=0; 

% set xylose as the carbon source
% lb(glc)=0;
% lb(xyl)=-vg;

% set oxygen uptake
% lb(o2)=-vo2;

% set PHB flux
% Fix a chosen flux to PHB synthesis, by setting the upper and lower bound, in order to explore the trade-off between biomass and PHB formation.
% ub(phb)=0;
% lb(phb)=0;

% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
% Changing the objective function to PHB synthesis for the case of
% non-growth associated production simulations, just uncomment the block of 
% code for the model you want to use:


% E. coli iML1515 model:

% objvector = model.c; % variable with the objective vector from the model
% objvector(2669) = 0; % growth is no longer the objective function
% objvector(2717) = 1; % PHB synthesis as the objective function



% C. necator H16 RehMBEL1391 model:

% objvector = model.c; % variable with the objective vector from the model
% objvector(21) = 0; % growth is no longer the objective function
% objvector(60) = 1; % sets PHB synthesis as the objective function



% S. cerevisiae iMM904 model:

% objvector = model.c; % variable with the objective vector from the model
% objvector(1521) = 0; % growth is no longer the objective function
% objvector(1581) = 1; % sets PHB synthesis as the objective function



% Completing the objective function change, uncomment the following line:

f=-objvector; %changes the objective to PHB synthesis
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% Running the FBA simulation:

% FBA simulation for linear objective function:
    [v,~,exitflag] = linprog(f,A,b,Aeq,beq,lb,ub);
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% Avoids negative concentrations.
% If the user removes or adds metabolites to be tracked, the respective
% removals and additions have to be also done here.

if exitflag == -2
 mi=0.0;
 vam=0.0;
 vphb=0.0;
 vac=0.0;
 vlac=0.0;
 veth=0.0;
 vform=0.0;
 vsucc=0.0;

else
mi=v(obj);
vam=v(am);
vphb=v(phb);
vac=v(ac);
vlac=v(lac);
veth=v(eth);
vform=v(form);
vsucc=v(succ);
end
% ------------------------------------------------------------------------
end