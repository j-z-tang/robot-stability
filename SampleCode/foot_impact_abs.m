    function [impact,terminate,direction] = foot_impact_abs(~,q,p)
    gamma = p.gamma;
        % This function detects a zero crossing, therefore the impact variable
        % needs to be zero when the next foot touches the ramp.
        if q(1) > 0
            impact = (q(1)+q(2)) + 2*gamma;
%                 impact = q(2)  - phaseFinal;
        else
            impact = 1;
        end
        terminate = 1;
        direction = -1;
    end