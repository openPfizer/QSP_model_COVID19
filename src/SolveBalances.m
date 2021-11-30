function [t,x] = SolveBalances(TSTART,TSTOP,Ts,data_dictionary)

  % Load DataFile struct -
  TSIM = TSTART:Ts:TSTOP;
  initial_condition_vector = data_dictionary.initial_condition;
  % opts = odeset('NonNegative',1:length(initial_condition_vector),'Events', @covidEventsFcn);
  opts = odeset('Events', @covidEventsFcn);
  
  avfun = @(tv) 0;
  args = {avfun};
  if find(ismember(fieldnames(data_dictionary),'pkav'))    
    % mean clinical PK approximation
    tpkav = data_dictionary.pkcont + data_dictionary.parameters.delay*24;
    avpk = data_dictionary.pkcon;
    avfun = @(tv) interp1(tpkav, avpk, tv);
    args = {avfun};
  else
    avfun = @(tv) 0;
    args = {avfun};
    data_dictionary.parameters.IC50_av = 1e5;
  end

  % Call the ODE solver -
  [t,x,te,ye,ie] = ode15s(@(t,x) covid19_dxdt(t,x,data_dictionary,args{:}),TSIM,initial_condition_vector, opts); 

  t1 = [];
  x1 = [];
    if (t(end)<TSTOP) && ~isempty(ye)
      % disp('enter')
      TSIM = t(end):Ts:TSTOP+2*Ts;
      opts = odeset('Events', []);
      data_dictionary.initial_condition = x(end,:);
      initial_condition_vector = data_dictionary.initial_condition;
      tmp_dictionary = data_dictionary;
      tmp_dictionary.parameters.k_int = 0; % load virus - 100/uL -> 2e5
      % tmp_dictionary.parameters.A_V = 0;
      [t1,x1] = ode15s(@(t,x) covid19_dxdt(t,x,tmp_dictionary,args{:}),TSIM,initial_condition_vector, opts);
      % disp('exit')
    end

  t = [t;t1(2:end)];
  x = [x;x1(2:end,:)];

end
