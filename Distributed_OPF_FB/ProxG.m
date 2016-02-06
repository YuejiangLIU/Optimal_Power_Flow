function p = ProxG(lb, ub)
% modified from cvxgrp@github 
p = @projection_box;
    function v = projection_box(x)
        v = max(lb, min(x, ub));
    end
end