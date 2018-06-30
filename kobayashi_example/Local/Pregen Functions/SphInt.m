% Scalar flux integration error
function [num_error num_sol] = SphInt(quad, refSol, prob)

[num_error num_sol] = calcError(quad, refSol, prob);

end