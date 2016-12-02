function [jointStateProbs] = StochSim(migrate,coal_rate,t)
% simulate joint state probabilities using Gillespies next time method

% the number of times the simulation must "passed" the end time
reps = 100000;

% initialize the joint state probabilities
jointStateProbs = zeros(2,2,0);

% do the simulations for all different times in t
for ti = 1 : length(t)
    % number of runs that coalesced before t(ti)
    failed=0;
    % number of runs that did not coalesce before t(ti)
    made=0;
    % inititialize the matrix capruting the joint state probs at time=t(ti)
    jointStateProbs_tmp = zeros(2,2);

    disp(ti)
    endtime=t(ti);
    
    for i = 1 : reps
        repeat = true;
        while repeat
            % initial states
            s1 = 1;
            s2 = 0;
            time =  0;
            repeat=false;
            while time <endtime
                
                mig1 = migrate;
                mig2 = migrate;

                if s1==s2
                    coal = coal_rate;
                else
                    coal = 0.0;
                end
                rate = 2*migrate+coal;
                time = time-1/(rate)*log(rand);
                r = rand;
                s1_old = s1;
                s2_old = s2;
                
                if r < migrate/rate % s1 migrates
                    s1 = abs(s1-1); 
                elseif r < (2*migrate)/rate && r >= migrate/rate % s2 migrates
                    s2 = abs(s2-1);
                else % coalescence
                    if time > endtime % if endtime was reached
                        break;
                    else % if the coalescence event happend before t(ti)
                        repeat=true;
                        break;
                    end
                end
            end
            if ~repeat                
                made = made+1;                
                jointStateProbs_tmp(s1_old+1,s2_old+1) = ...
                    jointStateProbs_tmp(s1_old+1,s2_old+1) + 1;
            else
                failed = failed + 1;
            end
        end
    end
    jointStateProbs(:,:,end+1) = jointStateProbs_tmp./(made+failed);
end
