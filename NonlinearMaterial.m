function [sigma, E_tan] = NonlinearMaterial(epsilon_eff, E0, sigma_max)

epsilon_eff = max(epsilon_eff, 0);

denom = 1 + (E0 * epsilon_eff) / sigma_max;

% Stress
sigma = (E0 * epsilon_eff) / denom;

% Tangent modulus: d(sigma)/d(eps) = E0 / denom^2
E_tan = E0 / denom^2;

end