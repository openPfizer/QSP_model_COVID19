% Drive file
% hr time units


% define simulation times
function [T,X] = function_run_model_noplots(data_dictionary,virus_inoculation,IC_dAT)
	
	TSTART = 0;
	TSTOP = 500*24; % n days to hours
	Ts = 0.1;

	TSIM = (TSTART:Ts:TSTOP);
	
	data_dictionary.parameters.hAT1 = data_dictionary.initial_condition(2);
	data_dictionary.parameters.hAT2 = data_dictionary.initial_condition(3);
	data_dictionary.parameters.hVI = 0;

	[T,X] = SolveBalances(TSTART,TSTOP,Ts,data_dictionary);
	data_dictionary.initial_condition = X(end,:);
	data_dictionary.initial_condition(1) = virus_inoculation; % load virus
	data_dictionary.parameters.hAT1 = data_dictionary.initial_condition(2);
	data_dictionary.parameters.hAT2 = data_dictionary.initial_condition(3);
	data_dictionary.parameters.hVI = 1;
	
	switch data_dictionary.pknab

		case 'REGN'
			data_dictionary.initial_condition(40) = data_dictionary.nAb0; % REGN initial condition casirivimab
			data_dictionary.initial_condition(53) = data_dictionary.nAb0; % REGN initial condition imdevimab
		case 'LY'
			data_dictionary.initial_condition(40) = data_dictionary.nAb0;	% Lilly initial condition bamlanivimab
			data_dictionary.initial_condition(53) = data_dictionary.nAb0;	% Lilly initial condition etesevimab
		otherwise
			data_dictionary.initial_condition(40) = 0;	% no initial condition
			data_dictionary.initial_condition(53) = 0;	% no initial condition
	end

	if IC_dAT ~= 0
		data_dictionary.initial_condition(5) = IC_dAT; % load virus
		data_dictionary.initial_condition(6) = IC_dAT; % load virus
	end
	% % 
	% % call solver
	TSTART = 0;
	TSTOP2 = 75*24;
	T2=0;
	 [T2,X2] = SolveBalances(TSTART,TSTOP2,Ts,data_dictionary);
	% % 

	
	T_basal = T((TSTOP-25*24)/Ts:end)-TSTOP;
	X_basal = X((TSTOP-25*24)/Ts:end,:);
	T_virus = [T_basal;T2];
	X_virus = [X_basal;X2];
	T = [T_virus];
	X = [X_virus];


end