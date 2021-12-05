%% driver_read_simulation
% Read equations and species from the equations file, then set values. 
% Will run and plot simulation results. 
%% Initializing

clear all;
% Default parameters
kon_standard = 1e-1; % 1/(min * M)
koff_standard = 0.1; % min^-1

%% Set species indices and initial conditions
% e.g., sp.R1=1, means the value for receptor 1 will be stored as the first
% index in our initial conditions below (y0), the first entry in
% our differential equations (dydt in 'eqns_read'), and the first column of
% our simulation results (Y).

% read file as char
eqns_text = fileread("eqns_read.m");

% find species names in equations file
species_all = regexp(eqns_text,'(?<=sp\.)\w*','match');
species_all=unique(species_all);
%   Documentation and examples for using regular expressions: 
%      https://www.mathworks.com/help/matlab/ref/regexp.html
%      I also looked at examples on stackoverflow to understand the usage.
%   Translating '(?<=sp\.)\w*':
%      I'm using the regexp look-around assertion '(?<=test)expr'
%      This looks for an expression (expr) that's behind a defined value (test)
%      '?<=' = look for something behind...
%      'sp\.' = 'sp.' (what I'm looking behind)
%               I need the '\' to make sure the '.' is literally
%               interpretted as a period.
%      '\w*' = Matches a word of any length (includeing underscores and #'s)
%              This command will stop at spaces and other symbols. 

% store species values (the indices) as a struct
for spi = 1:length(species_all)
    sp.(species_all{spi})=spi;
end

%% Set initial conditions

% Set initial conditions all to 0 
y0 = zeros(1,max(cell2mat(struct2cell(sp))));

% Now setting some species to be non-zero initially
y0(sp.La) = 10; %M (molar)
y0(sp.Lb) = 10; %M
y0(sp.Lc) = 10; %M
y0(sp.R1) = 10; %M 
y0(sp.R2) = 10; %M

%% Find parameter values

% find parameter names in equations file
params_all = regexp(eqns_text,'(?<=pa\.)\w*','match');
params_all=unique(params_all);
% see note on regexp in species section

%% Set parameter values
% for simplicity, I am setting all kon and koff values for each receptor &
% ligand to be equivalent, but you can change this 
koff_pnames = params_all(contains(params_all,"koff"));
for kind = 1:length(koff_pnames)
    pa.(koff_pnames{kind}) = koff_standard;
end

kon_pnames = params_all(contains(params_all,"kon"));
for kind = 1:length(kon_pnames)
    pa.(kon_pnames{kind}) = kon_standard;
end

% UNCOMMENT OR ADD LINES TO SET A SPECIFIC PARAMETER VALUES
 % pa.koff_R1_La = 0.5; 


%% assignment check
% Checking if you have assigned all the parameters in the equations file
if length(intersect(fieldnames(pa),params_all')) ~= length(params_all')
   error("Error. You have not set: \n%s",...
          strjoin(setdiff(params_all',fieldnames(pa)),'\n'));  
end

% Checking if you have assigned all the species in the equations file
if length(fieldnames(sp)) ~= length(y0)
   error("Error. \n Missing %i initial conditions",...
          length(fieldnames(sp))-length(y0) );  
end

%% Run simulation
tspan = [0 10]; %minutes
options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
[T,Y] = ode45(@eqns_read,tspan,y0,options,sp,pa);


%% Plot results
figure(); set(gcf, 'Color', 'w'); 
tiledlayout('flow')
for spi = 1:size(Y,2) 
    nexttile
    plot(T,Y(:,spi),'LineWidth',2)
    xlabel("minutes")
    ylabel("molar")
    ylim([0, max(Y,[],'all')]);
    title(species_all{spi},'Interpreter','none')
    % also useful for reformatting for interpreter (try this out) 
 %  title(regexprep(regexprep(species_all{spi},'_(\w*)','_{$1}'),"\{(\w+)_(\w+)\}",'\{$1,$2\}'),'Interpreter','tex')
    % % Notes: 
    % %  inner regexprep adds curly braces around text following the first underscore,
    % %  outer regexprep replaces any underscores that are between curly braces with a comma 
end
sgtitle("Concentrations of Species");
