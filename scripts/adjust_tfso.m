% Adjust time from symptom onset
function statep = adjust_tfso(state, tvmax_3s, T_sampleav)

   statep = -1*ones(10,length(tvmax_3s));
    for ii = 1:numel(tvmax_3s)
        % Get time vector indices for trial time points
        ind = [find(T_sampleav/24<=(tvmax_3s(ii)),1,'last');
        find(T_sampleav/24<=(tvmax_3s(ii))+2,1,'last');
        find(T_sampleav/24<=(tvmax_3s(ii))+4,1,'last');
        find(T_sampleav/24<=(tvmax_3s(ii))+6,1,'last');
        find(T_sampleav/24<=(tvmax_3s(ii))+7,1,'last');
        find(T_sampleav/24<=(tvmax_3s(ii))+9,1,'last');
        find(T_sampleav/24<=(tvmax_3s(ii))+10,1,'last');
        find(T_sampleav/24<=(tvmax_3s(ii))+12,1,'last');
        find(T_sampleav/24<=(tvmax_3s(ii))+15,1,'last')
        find(T_sampleav/24<=(tvmax_3s(ii))+20,1,'last')];

        statep(:,ii) = state(ind,ii);
    end

end