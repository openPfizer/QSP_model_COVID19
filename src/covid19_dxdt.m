function dxdt =  covid19_dxdt(t,x,data_dictionary, varargin)
% SECTION - Define parameters and species 	
  	% unpack plasma and lung volume
  	vol_plasma = data_dictionary.vol_plasma;
  	vol_alv_ml = data_dictionary.vol_alv_ml;

	p = data_dictionary.parameters;
	for species_index = 1:length(data_dictionary.initial_condition)
	  	s.(data_dictionary.species_names(species_index,2)) = x(species_index);
	  	IC.(data_dictionary.species_names(species_index,2)) = data_dictionary.initial_condition(species_index);
	end
	% unpack species array into meaningful variables	
	
	V = s.V;
	AT1 = s.AT1;
	AT2 = s.AT2;
	I = s.I;
	dAT1 = s.dAT1;
	dAT2 = s.dAT2;
	pDC = s.pDC;
	M1 = s.M1;
	N = s.N;
	Th1 = s.Th1;
	Th17 = s.Th17;
	CTL = s.CTL;
	Treg = s.Treg;
	TNFa = s.TNFa;
	IL6 = s.IL6;
	IL1b = s.IL1b;
	IFNb = s.IFNb;
	IFNg = s.IFNg;
	IL2 = s.IL2;
	IL12 = s.IL12;
	IL17 = s.IL17;
	IL10 = s.IL10;
	TGFb = s.TGFb;
	GMCSF = s.GMCSF;
	SPD = s.SPD;
	FER = s.FER;
	TNFa_c = s.TNFa_c;
	IL6_c = s.IL6_c;
	IL1b_c = s.IL1b_c;
	IFNb_c = s.IFNb_c;
	IFNg_c = s.IFNg_c;
	IL2_c = s.IL2_c;
	IL12_c = s.IL12_c;
	IL17_c = s.IL17_c;
	IL10_c = s.IL10_c;
	TGFb_c = s.TGFb_c;
	GMCSF_c = s.GMCSF_c;
	SPD_c = s.SPD_c;
	FER_c = s.FER_c;
	Ab = s.Ab;
	CRPExtracellular = s.CRPExtracellular;
	Blood_CRP = s.Blood_CRP;
	pDC1 = s.pDC1;	
	Ab_Igx = s.Ab_Igx;
	Ab12 = s.Ab12;
	Ab21 = s.Ab21;
	Ab22 = s.Ab22;

	
	% unpack parameters
	A_V = p.A_V;
	k_int = p.k_int;
	km_int_IFNb = p.km_int_IFNb;
	
	b_V = p.b_V;
	mu_AT2 = p.mu_AT2;
	k_ROS_AT2 = p.k_ROS_AT2;
	
	k_AT1_AT2 = p.k_AT1_AT2;
	km_AT1_AT2 = p.km_AT1_AT2;
	k_IFNb_kill = p.k_IFNb_kill;
	k_kill = p.k_kill;
	km_ROS_AT2 = p.km_ROS_AT2;
	b_AT2 = p.b_AT2;
	b_I = p.b_I;
	mu_AT1 = p.mu_AT1;
	
	b_AT1 = p.b_AT1;
	b_dAT1 = p.b_dAT1;
	k_damage_TNFa = p.k_damage_TNFa;
	km_damage_TNFa = p.km_damage_TNFa;
	k_damage_IL6 = p.k_damage_IL6;
	km_damage_IL6 = p.km_damage_IL6;
	k_damage_IL1b = p.k_damage_IL1b;
	km_damage_IL1b = p.km_damage_IL1b;
	k_damage_IFNg = p.k_damage_IFNg;
	km_damage_IFNg = p.km_damage_IFNg;
	k_damage_cyt = p.k_damage_cyt;
	a_DC = p.a_DC;
	kbasal_DC = p.kbasal_DC;
	b_DC = p.b_DC;
	lb_DC = p.lb_DC;
	k_DC_TNFa = p.k_DC_TNFa;
	km_DC_TNFa = p.km_DC_TNFa;
	k_DC_IFNg = p.k_DC_IFNg;
	km_DC_IFNg = p.km_DC_IFNg;
	k_DC_IL6 = p.k_DC_IL6;
	km_DC_IL6 = p.km_DC_IL6;
	km_DC_IL10 = p.km_DC_IL10;
	a_M1 = p.a_M1;
	kbasal_M1 = p.kbasal_M1;
	k_v = p.k_v;
	k_I = p.k_I;


	b_M1 = p.b_M1;
	lb_M1 = p.pb_M1;
	k_M1_TNFa = p.k_M1_TNFa;
	km_M1_TNFa = p.km_M1_TNFa;
	k_M1_GMCSF = p.k_M1_GMCSF;
	km_M1_GMCSF = p.km_M1_GMCSF;
	k_M1_IFNg = p.k_M1_IFNg;
	km_M1_IFNg = p.km_M1_IFNg;
	km_M1_IL10 = p.km_M1_IL10;
	a_N = p.a_N;
	k_N_IFNg = p.k_N_IFNg;
	km_N_IFNg = p.km_N_IFNg;
	k_N_TNFa = p.k_N_TNFa;
	km_N_TNFa = p.km_N_TNFa;
	k_N_GMCSF = p.k_N_GMCSF;
	km_N_GMCSF = p.km_N_GMCSF;
	k_N_IL17c = p.k_N_IL17c;
	km_N_IL17c = p.km_N_IL17c;
	b_N = p.b_N;
	a_Th1 = p.a_Th1;
	b_Th1 = p.b_Th1;
	k_Th1_IL2 = p.k_Th1_IL2;
	K_Th1_IL12 = p.K_Th1_IL12;
	k_Th1_IL12IL2 = p.k_Th1_IL12IL2;
	K_Th1_IL12IL2 = p.K_Th1_IL12IL2;
	K_Th1_IL10 = p.K_Th1_IL10;
	K_Th1_TGFb = p.K_Th1_TGFb;
	k_Th1_IFNg = p.k_Th1_IFNg;
	K_IFNg_Th1 = p.K_IFNg_Th1;
	K_Th1_IL6 = p.K_Th1_IL6;
	k_Th1_Th17 = p.k_Th1_Th17;
	K_Th1_Th17 = p.K_Th1_Th17;
	k_Th1_Treg = p.k_Th1_Treg;
	K_Th1_Treg = p.K_Th1_Treg;
	a_Th17 = p.a_Th17;
	b_Th17 = p.b_Th17;
	k_Th17_TGFb = p.k_Th17_TGFb;
	K_Th17_TGFb = p.K_Th17_TGFb;
	K_Th17_IL2 = p.K_Th17_IL2;
	K_Th17_IFNg = p.K_Th17_IFNg;
	K_Th17_IL10 = p.K_Th17_IL10;
	k_Th17_IL6 = p.k_Th17_IL6;
	km_Th17_IL6 = p.km_Th17_IL6;
	k_Th17_IL1b = p.k_Th17_IL1b;
	km_Th17_IL1b = p.km_Th17_IL1b;
	a_CTL = p.a_CTL;
	b_CTL = p.b_CTL;
	k_CTL_IL2 = p.k_CTL_IL2;
	K_CTL_IL12 = p.K_CTL_IL12;
	k_CTL_IL12IL2 = p.k_CTL_IL12IL2;
	K_CTL_IL12IL2 = p.K_CTL_IL12IL2;
	K_CTL_IL10 = p.K_CTL_IL10;
	K_CTL_TGFb = p.K_CTL_TGFb;
	K_CTL_IL6	 = p.K_CTL_IL6	;
	k_CTL_IFNg = p.k_CTL_IFNg;
	K_CTL_IFNg = p.K_CTL_IFNg;
	kmax_MHC1 = p.kmax_MHC1;
	km_MHC1_IFNb = p.km_MHC1_IFNb;
	a_Treg = p.a_Treg;
	b_Treg = p.b_Treg;
	k_Treg_IL2 = p.k_Treg_IL2;
	K_Treg_IL2 = p.K_Treg_IL2;
	K_Treg_IL17 = p.K_Treg_IL17;
	K_Treg_IL6 = p.K_Treg_IL6;
	k_Treg_TGFb = p.k_Treg_TGFb;
	K_Treg_TGFb = p.K_Treg_TGFb;
	kbasal_SPD = p.kbasal_SPD;
	a_SPD = p.a_SPD;
	b_SPD = p.b_SPD;
	kbasal_FER = p.kbasal_FER;
	a_FER = p.a_FER;
	b_FER = p.b_FER;
	a_tnf = p.a_tnf;
	a_tnf_at1 = p.a_tnf_at1;
	a_tnf_i = p.a_tnf_i;
	a_tnf_at2 = p.a_tnf_at2;
	a_tnf_m1 = p.a_tnf_m1;
	a_tnf_th1 = p.a_tnf_th1;
	a_tnf_th17 = p.a_tnf_th17;
	b_tnf = p.b_tnf;
	a_il6 = p.a_il6;
	b_il6 = p.b_il6;
	a_il6_at1 = p.a_il6_at1;
	a_il6_i = p.a_il6_i;
	a_il6_at2 = p.a_il6_at2;
	a_il6_m1 = p.a_il6_m1;
	a_il6_th17 = p.a_il6_th17;
	a_il6_neu = p.a_il6_neu;
	a_ifng = p.a_ifng;
	b_ifng = p.b_ifng;
	a_ifng_dc = p.a_ifng_dc;
	a_ifng_th1 = p.a_ifng_th1;
	a_ifng_ctl = p.a_ifng_ctl;
	a_ifnb = p.a_ifnb;
	b_ifnb = p.b_ifnb;
	a_ifnb_at1 = p.a_ifnb_at1;
	a_ifnb_i = p.a_ifnb_i;
	a_ifnb_d = p.a_ifnb_d;
	a_ifnb_dc = p.a_ifnb_dc;
	a_il2_dc = p.a_il2_dc;
	a_il2_th1 = p.a_il2_th1;
	b_il2 = p.b_il2;
	a_il2 = p.a_il2;
	a_il12_m1 = p.a_il12_m1;
	a_il12_dc = p.a_il12_dc;
	b_il12 = p.b_il12;
	a_il12 = p.a_il12;
	a_il17_th17 = p.a_il17_th17;
	a_il17_ctl = p.a_il17_ctl;
	b_il17 = p.b_il17;
	a_il17 = p.a_il17;
	a_il10_treg = p.a_il10_treg;
	b_il10 = p.b_il10;
	a_il10 = p.a_il10;
	a_tgfb_th17 = p.a_tgfb_th17;
	a_tgfb_treg = p.a_tgfb_treg;
	b_tgfb = p.b_tgfb;
	a_tgfb = p.a_tgfb;
	a_gmcsf_m1 = p.a_gmcsf_m1;
	a_gmcsf_th1 = p.a_gmcsf_th1;
	a_gmcsf_th17 = p.a_gmcsf_th17;
	b_gmcsf = p.b_gmcsf;
	a_gmcsf = p.a_gmcsf;
	a_il1b = p.a_il1b;
	b_il1b = p.b_il1b;
	a_il1b_at1 = p.a_il1b_at1;
	a_il1b_i = p.a_il1b_i;
	a_il1b_at2 = p.a_il1b_at2;
	a_il1b_m1 = p.a_il1b_m1;
	a_il1b_dc = p.a_il1b_dc;
	kbasal_ROS = p.kbasal_ROS;
	b_ROS = p.b_ROS;
	ktr_TNFa = p.ktr_TNFa;
	ktr_IL6 = p.ktr_IL6;
	ktr_IL1b = p.ktr_IL1b;
	ktr_IFNb = p.ktr_IFNb;
	ktr_IFNg = p.ktr_IFNg;
	ktr_IL2 = p.ktr_IL2;
	ktr_IL12 = p.ktr_IL12;
	ktr_IL17 = p.ktr_IL17;
	ktr_IL10 = p.ktr_IL10;
	ktr_TGFb = p.ktr_TGFb;
	ktr_GMCSF = p.ktr_GMCSF;
	ktr_SPD = p.ktr_SPD;
	ktr_FER = p.ktr_FER;
	k_CTL_I_SPD = p.k_CTL_I_SPD;
	k_CTL_I_Fer = p.k_CTL_I_Fer;

	%%%% LUNG TO CENTRAL TRANSPORT RATES
	tr_TNFa = ktr_TNFa*TNFa;
	tr_IL6 = ktr_IL6*IL6;
	tr_IL1b = ktr_IL1b*IL1b;
	tr_IFNb = ktr_IFNb*IFNb;
	tr_IFNg = ktr_IFNg*IFNg;
	tr_IL2 = ktr_IL2*IL2;
	tr_IL12 = ktr_IL12*IL12;
	tr_IL17 = ktr_IL17*IL17;
	tr_IL10 = ktr_IL10*IL10;
	tr_TGFb = ktr_TGFb*TGFb;
	tr_GMCSF = ktr_GMCSF*GMCSF;
	tr_SPD = ktr_SPD*SPD;
	tr_FER = ktr_FER*FER;


	tr_pDC = p.ktr_pDC*(s.pDC/vol_alv_ml -s.pDC_c);
	tr_M1 =  p.ktr_M1*(s.M1/vol_alv_ml -s.M1_c);
	tr_N = p.ktr_N*(s.N/vol_alv_ml -s.N_c);
	tr_Th1 = p.ktr_Th1*(s.Th1/vol_alv_ml -s.Th1_c);
	tr_Th17 = p.ktr_Th17*(s.Th17/vol_alv_ml -s.Th17_c);
	tr_CTL = p.ktr_CTL*(s.CTL/vol_alv_ml -s.CTL_c);
	tr_Treg = p.ktr_Treg*(s.Treg/vol_alv_ml -s.Treg_c);
% !SECTION - end Define parameters and species 

% 	NOTE - hVI to ensure virus is not shed during healthy subject simulation (due to numerical error) 
	hVI = p.hVI;
%% SECTION - Initializing antiviral (simplistic infusion or custom PK)  
	Cav = varargin{1};
	Cavint = Cav(t);
	Cavint(isnan(Cavint) || Cavint<=0) = 0;
	IC50_av = p.IC50_av;

% !SECTION - end Initializing antiviral (simplistic infusion or custom PK)

%% SECTION - Initialize endogenou Ab and exogenous nAb with Lilly as default and Imax = 0
	Imaxbnab = 0; % default Imax for nAbs
	IC50_bnAb1 = 0.03;
	IC50_bnAb2 = 0.0900;
	CL_Aba = 0.27/24; % convert units to d 
	Qa =0.375/24; % convert units to d
	V_Ab1a = 2.87;
	V_Ab2a = 2.71;
	CL_Abb = 0.128/24; % convert units to d
	Qb =0.514/24; % convert units to d
	V_Ab1b = 2.38;
	V_Ab2b =1.98;
	ab_clearance = 0;
	if find(ismember(fieldnames(data_dictionary),'pknab'))
		switch data_dictionary.pknab

			case 'REGN'
				% REGN-COV2 neuralizing antibody
				CL_Aba = 0.22/24;
				Qa = 0;
				V_Ab1a = 4.5;
				V_Ab2a = 4.5;
				CL_Abb = 0.128/24; % convert units to d
				Qb = 0;
				V_Ab1b = 4.5;
				V_Ab2b = 1.98;
				IC50_bnAb2 = 0.0561; %ug/mL
				IC50_bnAb1 = 0.0631; %ug/mL
				Imaxbnab = 0.28; % REGN Imax
			case 'LY'

				% % Lilly neutralizing antibody parameters
				CL_Aba = 0.27/24; % convert units to d
				Qa =0.375/24; % convert units to d
				V_Ab1a = 2.87;
				V_Ab2a = 2.71;
				CL_Abb = 0.128/24; % convert units to d
				Qb =0.514/24; % convert units to d
				V_Ab1b = 2.38;
				V_Ab2b =1.98;
				IC50_bnAb1 = 0.3; %ug/mL
				IC50_bnAb2 = 0.09; %ug/mL
				Imaxbnab = 0.30; % Lilly Imax
		end
	end
	
	% SECTION - exogenous neutralizing antibody dynamics
	dAbdt = 0;
	dAb12dt = 0;
	dAb21dt = 0;
	dAb22dt = 0;
	CONC1 = 0;
	CONC2 = 0;

	if t>=p.tau_Ab && hVI
		CONC1 = Ab/V_Ab1a;
		CONC2 = Ab21/V_Ab1b;
		dAbdt =  - Ab*CL_Aba/V_Ab1a - Qa/V_Ab1a*Ab + Qa/V_Ab2a*Ab12;
		dAb12dt = Qa/V_Ab1a*Ab - Qa/V_Ab2a*Ab12;
		dAb21dt = - Ab21*CL_Abb/V_Ab1b - Qb/V_Ab1b*Ab21 + Qb/V_Ab2b*Ab22;
		dAb22dt = Qb/V_Ab1b*Ab21 - Qb/V_Ab2b*Ab22;
	end	

	% !SECTION - end exogenous neutralizing antibody dynamics

	% SECTION endogenous antibody dynamics
		dAb_Igxdt = 0;
		IC50_Igx = 0.03;
		ImaxIgx = 1;
		thalf_Igx = 51; % h
		CL_Ab_Igx = 0.693/thalf_Igx;
		V_Ab_Igx = 3.5;
		if t>=p.tau_Ab_Igx && hVI
			k_prod_igx = 2; % endogenous Ab production
			dAb_Igxdt =  k_prod_igx - Ab_Igx*CL_Ab_Igx/V_Ab_Igx;
			ImaxIgx = 1; % endogenous Ab imax 
			IC50_Igx = 3; % endogenous Ab IC50
		end	
		CONC_Igx = Ab_Igx/V_Ab_Igx;
	% !SECTION - end endogenous antibody dynamics
% !SECTION - end Initialize endogenou Ab and exogenous nAb

%% SECTION - V - Viral dynamics (# viral mRNA copies/mL)
	% p.nhill = 1;
	nhillab = 1;
	prod_virus_shedding = A_V*I*(1-p.Imax*Cavint^p.nhill/(IC50_av^p.nhill+Cavint^p.nhill));
	virus_endocytosis = k_int*AT2*V*(km_int_IFNb/(km_int_IFNb+IFNb))*(1-Imaxbnab*(CONC1^nhillab/(CONC1^nhillab + IC50_bnAb1^nhillab)+CONC2^nhillab/(CONC2^nhillab+IC50_bnAb2^nhillab)))*(1-ImaxIgx*CONC_Igx/(CONC_Igx+IC50_Igx));
	deg_virus = b_V*V;
	
	dVdt = hVI*(prod_virus_shedding  - deg_virus - ab_clearance);
% !SECTION - end

%% SECTION - AT2 - healthy alveolar type 2 cells (AT2) (# cells)
	ICAT1 = p.hAT1; % Homeostatic level of healthy AT1 cells
	ICAT2 = p.hAT2; % Homeostatic level of susceptible AT2 cells
    damage_cyt_AT = k_damage_cyt*(k_damage_TNFa*(TNFa/(km_damage_TNFa+TNFa)) +  k_damage_IL6*(IL6/(km_damage_IL6+IL6)) ...
        + k_damage_IL1b*(IL1b/(km_damage_IL1b + IL1b)) + k_damage_IFNg*(IFNg/(km_damage_IFNg+IFNg)));
	damage_ROS_AT2 = k_ROS_AT2*N/(km_ROS_AT2+N);	
	growth_AT2 = mu_AT2*(1 + p.k_mu_AT2*(max((ICAT1+ICAT2)-(AT2+AT1),0)/(km_AT1_AT2*(ICAT1+ICAT2) + max((ICAT1+ICAT2)-(AT2+AT1),0))))*AT2; % induced differentiation and growth induced by AT1 decrease from baseline
	deg_AT2 = b_AT2*AT2;
	diff_AT2 = k_AT1_AT2 * AT2 * (1 + p.k_diff_AT1*max(ICAT1-AT1,0)/(km_AT1_AT2*ICAT1 + max(ICAT1-AT1,0)));
	dAT2dt = growth_AT2 - virus_endocytosis - damage_ROS_AT2*AT2 - deg_AT2 - damage_cyt_AT*AT2 - diff_AT2;
% !SECTION - end

%% SECTION I - Infected AT2 (# cells)
	kill_CTL_I = k_kill*I*(CTL);
	deg_I = b_I*I*(1+k_IFNb_kill*IFNb); %/(1*km_kill+IFNb)

	dIdt = virus_endocytosis - deg_I - kill_CTL_I - damage_ROS_AT2*I;
% !SECTION - end

%% SECTION - AT1 - health alveolar type 1 cells (AT1) (# cells)
	growth_AT1 = mu_AT1*(ICAT1-AT1); % no growth, no division, source is from AT2 
	damage_ROS_AT1 = damage_ROS_AT2*AT1;
	deg_AT1 = b_AT1*AT1;

	dAT1dt = growth_AT1 - damage_ROS_AT1 - deg_AT1 - damage_cyt_AT*AT1 + diff_AT2 ;
% !SECTION - end

%% SECTION - dAT1 - damaged AT1 cells (# cells)
	deg_dAT1 = b_dAT1*dAT1;
	ddAT1dt = damage_ROS_AT1 - deg_dAT1 + damage_cyt_AT*AT1;
% !SECTION - end

%% SECTION - dAT2 - damaged AT2 cells (# cells)
	deg_dAT2 = b_dAT1*dAT2;
	ddAT2dt = damage_ROS_AT2*(I+AT2)  - deg_dAT2 + damage_cyt_AT*AT2;
% !SECTION - end



% % % % % Immune cells % % % % % 

%% SECTION - pDC - lung pulmonary dendritic cells (# cells)
    act_pDC_TNFa = k_DC_TNFa*(TNFa/(km_DC_TNFa+TNFa));
	act_pDC_IFNg = k_DC_IFNg*(IFNg/(km_DC_IFNg+IFNg));
    act_pDC_GMCSF = k_M1_GMCSF*(GMCSF/(km_M1_GMCSF+GMCSF));
	inh_pDC_IL10 = (km_DC_IL10/(km_DC_IL10+IL10));

    dpDCdt = a_DC*(kbasal_DC+ k_v*log(1+V) + k_I*log(1+I)+ p.k_dAT*log((dAT1+dAT2)))*(1+act_pDC_TNFa + act_pDC_IFNg + act_pDC_GMCSF)* inh_pDC_IL10 - lb_DC*b_DC*pDC - tr_pDC*vol_alv_ml;  
% !SECTION - end

%% SECTION - M1 - lung M1 macrophages (# cells)
	act_M1_TNFa = k_M1_TNFa*(TNFa/(km_M1_TNFa+TNFa));
	act_M1_GMCSF = k_M1_GMCSF*(GMCSF/(km_M1_GMCSF+GMCSF));
	act_M1_IFNg = k_M1_IFNg*(IFNg/(km_M1_IFNg+IFNg));
	inh_M1_IL10 = (km_M1_IL10/(km_M1_IL10+IL10));

	dM1dt = a_M1*(kbasal_M1+ k_v*log(1+V)+ k_I*log(1+I) + p.k_dAT*log((dAT1+dAT2)))*(1+ act_M1_TNFa + act_M1_GMCSF + act_M1_IFNg)*inh_M1_IL10 - b_M1*M1 - tr_M1*vol_alv_ml;
% !SECTION - end

%% SECTION - N - lung Neutrophils (# cells)
    act_N_IFNg = k_N_IFNg*IFNg/(IFNg+km_N_IFNg); 
    act_N_TNFa= k_N_TNFa*TNFa/(TNFa+km_N_TNFa); 
    act_N_GMCSF = k_N_GMCSF*GMCSF/(GMCSF+km_N_GMCSF);
	rec_N_IL17c = k_N_IL17c*IL17_c/(IL17_c + km_N_IL17c);
% !SECTION - end

%% SECTION - pDC1 - Delay in activation of adaptive immune cells by pDCs
% NOTE - pDCs are primary APCs in model 
	mtr = 0.0500; % transit compartment
	kdie = mtr;
	dpDC1dt = mtr*(pDC - kdie*pDC1);
% !SECTION - end

    dNdt = a_N*(p.kbasal_N+  k_v*log(1+V) +  k_I*log(1+I) + p.k_dAT*log((dAT1+dAT2)))*(1+act_N_IFNg + act_N_TNFa + act_N_GMCSF) + rec_N_IL17c - b_N*N - tr_N*vol_alv_ml ;

%% SECTION - Th1 - lung Th1 Cells (# cells)
	act_Th1_IL12 = k_Th1_IL2*(IL12/(K_Th1_IL12+IL12))*(1+k_Th1_IL12IL2*(IL2/(K_Th1_IL12IL2+IL2)));
	act_Th1_IFNg = k_Th1_IFNg*(IFNg/(K_IFNg_Th1 + IFNg)) * (K_Th1_IL6/(K_Th1_IL6+IL6));
	inh_Th1_IL10_TGFb = (K_Th1_IL10/(K_Th1_IL10+IL10)) * (K_Th1_TGFb/(K_Th1_TGFb+TGFb));
    diff_Th1_Th17 = k_Th1_Th17*Th17*(IL12/(K_Th1_Th17+IL12))*(K_Th1_TGFb/(K_Th1_TGFb+TGFb));
	diff_Th1_Treg = k_Th1_Treg*Treg*(IL12/(K_Th1_Treg+IL12));

	dTh1dt = a_Th1*pDC1*(act_Th1_IL12 + act_Th1_IFNg)*inh_Th1_IL10_TGFb + diff_Th1_Th17 + diff_Th1_Treg - b_Th1*Th1 - tr_Th1*vol_alv_ml;% 
% !SECTION - end

%% SECTION - Th17 - lung TH17 cells (# cells)
	inh_Th17 = (K_Th17_IL2/(K_Th17_IL2 + IL2)) * (K_Th17_IFNg/(K_Th17_IFNg + IFNg)) * (K_Th17_IL10/(K_Th17_IL10 + IL10));
    act_TH17_TGFb = k_Th17_TGFb*(TGFb/(K_Th17_TGFb + TGFb)) * inh_Th17;
	act_Th17_IL6 = k_Th17_IL6*(IL6/(km_Th17_IL6 + IL6)) * inh_Th17;
	act_Th17_IL1b = k_Th17_IL1b*(IL1b/(km_Th17_IL1b + IL1b)) * inh_Th17;
	dTh17dt = a_Th17*pDC1*(act_TH17_TGFb + act_Th17_IL6 + act_Th17_IL1b) - diff_Th1_Th17 - b_Th17*Th17 - tr_Th17*vol_alv_ml;
% !SECTION - end

%% SECTION - CTL - lung CTL (# cells)
	act_CTL_IFNb = kmax_MHC1*(IFNb/(km_MHC1_IFNb+IFNb));
	act_CTL_IL12 = k_CTL_IL2*(IL12/(K_CTL_IL12+IL12))*(1+k_CTL_IL12IL2*(IL12/(K_CTL_IL12IL2+IL12))) * (K_CTL_IL10/(K_CTL_IL10+IL10)) * (K_CTL_TGFb/(K_CTL_TGFb+TGFb));
	act_CTL_IFNg = k_CTL_IFNg*(IFNg/(K_CTL_IFNg + IFNg)) * (K_CTL_IL10/(K_CTL_IL10+IL10)) * (K_CTL_TGFb/(K_CTL_TGFb+TGFb)) * (K_CTL_IL6/(K_CTL_IL6+IL6));

	dCTLdt = a_CTL*pDC1*(1+act_CTL_IFNb)*(1+ act_CTL_IL12 + act_CTL_IFNg)  - b_CTL*CTL - tr_CTL*vol_alv_ml;
% !SECTION - end

%% SECTION - Treg - lung Treg (# cells)
	act_Treg_IL2 = k_Treg_IL2*(IL2/(K_Treg_IL2 + IL2)) * (K_Treg_IL17/(K_Treg_IL17 + IL17)) * (K_Treg_IL6/(K_Treg_IL6 + IL6));
	act_Treg_TGFb = k_Treg_TGFb*(TGFb/(K_Treg_TGFb + TGFb)) * (K_Treg_IL17/(K_Treg_IL17 + IL17)) * (K_Treg_IL6/(K_Treg_IL6 + IL6));

	dTregdt = a_Treg*pDC1*(act_Treg_IL2 + act_Treg_TGFb) - diff_Th1_Treg - b_Treg*Treg - tr_Treg*vol_alv_ml;
% !SECTION - end

%%%%% Biomarkers %%%%%

%% SECTION - SPD - lung surfactant protein D (pmol)
	dSPDdt = kbasal_SPD + a_SPD * (dAT1 + dAT2 + k_CTL_I_SPD*kill_CTL_I) - tr_SPD ;
% !SECTION - end

%% SECTION - FER - lung Ferritin (pmol)
	dFERdt = kbasal_FER + a_FER * (dAT1 + dAT2 + k_CTL_I_Fer*kill_CTL_I) - tr_FER;
% !SECTION - end

%% SECTION - CRP - Liver CRP Production (pmol)
    Liver_CRP = p.k_livercrp*p.Liver*((IL6_c*p.Liver));    
    prod_CRP_liver = p.kCRPSecretion*Liver_CRP;
    tr_CRP = p.kCRP_BloodtoLiver*Blood_CRP-p.kCRP_LivertoBlood*Liver_CRP/p.Liver;
    prod_CRP_blood = p.kbasal_CRP;
    dCRPExtracellular_dt = prod_CRP_liver + tr_CRP*p.Liver - p.kdeg_CRP*CRPExtracellular;
    dBlood_CRP_dt = (-tr_CRP + prod_CRP_blood - p.kdeg_CRP*Blood_CRP);
% !SECTION - end

%%%%% Cytokine Dynamics %%%%%

%% SECTION - TNF-A (pmol)
	prod_tnf_dat1 = a_tnf_at1*dAT1;
	prod_tnf_i = a_tnf_i*I;
	prod_tnf_dat2 = a_tnf_at2*dAT2;
	prod_tnf_m1 = a_tnf_m1*M1;
	prod_tnf_th1 = a_tnf_th1*Th1;
	prod_tnf_th17 = a_tnf_th17*Th17;
	deg_tnf = b_tnf*TNFa;
   

	dTNFadt = a_tnf*(p.basal_tnfa + prod_tnf_i + prod_tnf_dat2 + prod_tnf_dat1 + prod_tnf_m1 + prod_tnf_th1 +...
	         prod_tnf_th17) - deg_tnf - tr_TNFa;
	     
% !SECTION - end

%% SECTION - IL-6 (pmol)
	prod_il6_dat1 = a_il6_at1*dAT1;
	prod_il6_i = a_il6_i*I;
	prod_il6_dat2 = a_il6_at2*dAT2;
	prod_il6_m1 = a_il6_m1*M1;
	prod_il6_th17 = a_il6_th17*Th17;
	prod_il6_neu = a_il6_neu*N;
	deg_il6 = b_il6*IL6;

	dIL6dt = a_il6*(p.basalil6 + prod_il6_dat1 + prod_il6_i + prod_il6_dat2 + prod_il6_m1 + prod_il6_th17 +...
	         prod_il6_neu) - deg_il6 - tr_IL6;
% !SECTION - end

%% SECTION - IL1-B (pmol)
	prod_il1b_dat1 = a_il1b_at2*dAT1;
	prod_il1b_i = a_il1b_i*I;
	prod_il1b_dat2 = a_il1b_at2*dAT2;
	prod_il1b_m1 = a_il1b_m1*M1;
	prod_il1b_dc = a_il1b_dc*pDC;
	deg_il1b = b_il1b*IL1b;

	dIL1bdt = a_il1b*(p.basalil1 + prod_il1b_dat1 + prod_il1b_i + prod_il1b_dat2 + prod_il1b_m1 + ...
	          prod_il1b_dc) - deg_il1b - tr_IL1b;
% !SECTION - end

%% SECTION - IFN-G (pmol)
	prod_ifng_dc = a_ifng_dc*pDC;
	prod_ifng_th1 = a_ifng_th1*Th1;
	prod_ifng_ctl = a_ifng_ctl*CTL;
	del_ifng = b_ifng*IFNg;

	dIFNgdt = a_ifng*(p.basalifng + prod_ifng_dc + prod_ifng_th1 + prod_ifng_ctl) - del_ifng - tr_IFNg;
% !SECTION - end

%% SECTION - IFN-B (pmol)
	prod_ifnb_i = a_ifnb_i*I;
	prod_ifnb_dc = a_ifnb_dc*pDC;
	del_ifnb = b_ifnb*IFNb;

	dIFNbdt = a_ifnb*(p.basalifnb + prod_ifnb_i + prod_ifnb_dc) - del_ifnb - tr_IFNb;
% !SECTION - end

%% SECTION - IL-2 (pmol)
	prod_il2_dc = a_il2_dc*pDC;
	prod_il2_th1 = a_il2_th1*Th1;
	deg_il2 = b_il2*IL2;

	dIL2dt = a_il2*(p.basalil2 + prod_il2_dc + prod_il2_th1) - deg_il2 - tr_IL2;
% !SECTION - end

%% SECTION - IL-12 (pmol)
	prod_il12_m1 = a_il12_m1*M1;
	prod_il12_dc = a_il12_dc*pDC;
	deg_il12 = b_il12*IL12;

	dIL12dt = a_il12*(p.basalil12 + prod_il12_m1 + prod_il12_dc) - deg_il12 - tr_IL12;
% !SECTION - end

%% SECTION - IL-17 (pmol)
	prod_il17_th17 = a_il17_th17*Th17;
	prod_il17_ctl = a_il17_ctl*CTL;
	deg_il17 = b_il17*IL17;

	dIL17dt = a_il17*(prod_il17_th17 + prod_il17_ctl) - deg_il17 - tr_IL17;
% !SECTION - end

%% SECTION - IL-10 (pmol)
	prod_il10_treg = a_il10_treg*Treg;
	deg_il10 = b_il10*IL10;

	dIL10dt = a_il10*(p.basalil10 + prod_il10_treg) - deg_il10 - tr_IL10;
% !SECTION - end

%% SECTION - TGF-B (pmol)
	prod_tgfb_th17 = a_tgfb_th17*Th17;
	prod_tgfb_treg = a_tgfb_treg*Treg;
	deg_tgfb = b_tgfb*TGFb;

	dTGFbdt = a_tgfb*(p.basaltgfb + prod_tgfb_th17 + prod_tgfb_treg) - deg_tgfb - tr_TGFb;
% !SECTION - end

%% SECTION - GM-CSF (pmol)
	prod_gmcsf_m1 = a_gmcsf_m1*M1;
	prod_gmcsf_th1 = a_gmcsf_th1*Th1;
	prod_gmcsf_th17 = a_gmcsf_th17*Th17;
	deg_gmcsf = b_gmcsf*GMCSF;

	dGMCSFdt = a_gmcsf*(p.basalgmcsf + prod_gmcsf_m1 + prod_gmcsf_th1 + prod_gmcsf_th17) - deg_gmcsf - tr_GMCSF;
% !SECTION - end

%% SECTION - Central compartment cytokines (pmol/mL)
	dTNFa_cdt = tr_TNFa/vol_plasma - b_tnf*TNFa_c;
	dIL6_cdt = tr_IL6/vol_plasma - b_il6*IL6_c;
	dIL1b_cdt = tr_IL1b/vol_plasma - b_il1b*IL1b_c;
	dIFNb_cdt = tr_IFNb/vol_plasma - b_ifnb*IFNb_c;
	dIFNg_cdt = tr_IFNg/vol_plasma - b_ifng*IFNg_c;
	dIL2_cdt = tr_IL2/vol_plasma - b_il2*IL2_c;
	dIL12_cdt = tr_IL12/vol_plasma - b_il12*IL12_c;
	dIL17_cdt = tr_IL17/vol_plasma - b_il17*IL17_c;
	dIL10_cdt = tr_IL10/vol_plasma - b_il10*IL10_c;
	dTGFb_cdt = tr_TGFb/vol_plasma - b_tgfb*TGFb_c;
	dGMCSF_cdt = tr_GMCSF/vol_plasma - b_gmcsf*GMCSF_c;
	dSPD_cdt = tr_SPD/vol_plasma - b_SPD*SPD_c;
	dFER_cdt = p.a_FER_c + tr_FER/vol_plasma - b_FER*FER_c;
% !SECTION - end

%% SECTION - Central compartment cells (# cells/uL)
	dpDC_cdt =  tr_pDC*vol_alv_ml/vol_plasma - s.pDC_c*b_DC;
	dM1_cdt =  tr_M1*vol_alv_ml/vol_plasma - lb_M1*s.M1_c*b_M1;
	dN_cdt = tr_N*vol_alv_ml/vol_plasma -  s.N_c*b_N;
	dTh1_cdt =  tr_Th1*vol_alv_ml/vol_plasma - s.Th1_c*b_Th1;
	dTh17_cdt =  tr_Th17*vol_alv_ml/vol_plasma - s.Th17_c*b_Th17;
	dCTL_cdt =  tr_CTL*vol_alv_ml/vol_plasma - s.CTL_c*b_CTL;
	dTreg_cdt =  tr_Treg*vol_alv_ml/vol_plasma - s.Treg_c*b_Treg;
% !SECTION - end


%% SECTION - Assign Output
	dxdt = [
	dVdt
	dAT1dt
	dAT2dt
	dIdt
	ddAT1dt
	ddAT2dt
	dpDCdt
	dM1dt
	dNdt
	dTh1dt
	dTh17dt
	dCTLdt
	dTregdt
	dTNFadt
	dIL6dt
	dIL1bdt
	dIFNbdt
	dIFNgdt
	dIL2dt
	dIL12dt
	dIL17dt
	dIL10dt
	dTGFbdt
	dGMCSFdt
	dSPDdt
	dFERdt
	dTNFa_cdt
	dIL6_cdt
	dIL1b_cdt
	dIFNb_cdt
	dIFNg_cdt
	dIL2_cdt
	dIL12_cdt
	dIL17_cdt
	dIL10_cdt
	dTGFb_cdt
	dGMCSF_cdt
	dSPD_cdt
	dFER_cdt
	dAbdt
	dpDC_cdt
	dM1_cdt
	dN_cdt
	dTh1_cdt
	dTh17_cdt
	dCTL_cdt
	dTreg_cdt
	dCRPExtracellular_dt
	dBlood_CRP_dt
	dpDC1dt
	dAb_Igxdt
	dAb12dt
	dAb21dt
	dAb22dt];
% !SECTION - end
end
