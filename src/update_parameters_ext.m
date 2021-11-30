function [data_dictionary,virus_innoculation] = update_parameters_ext(data_dictionary,tmp_parameters, sample_index, perturbation_name_vector, lhs_samples) % Updates model dictionary; data_dictionary
	p_org = data_dictionary.parameters;

	for parameter_index = 1:length(perturbation_name_vector)
		tmp_parameters.(perturbation_name_vector(parameter_index)) = p_org.(perturbation_name_vector(parameter_index))*lhs_samples(sample_index,parameter_index);
	end
	data_dictionary.parameters = tmp_parameters;
	virus_innoculation = 1e1; % varying virus innoculation

end