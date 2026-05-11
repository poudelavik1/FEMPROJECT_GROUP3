function [sigma, E_tan] = NonlinearMaterial(epsilon_eff, E0, sigma_max)
%NONLINEARMATERIAL  Hyperbolic saturation material model
%
%   sigma = E0*epsilon / (1 + E0*epsilon/sigma_max)
% Inputs
%   epsilon_eff : effective strain (scalar or array, >= 0)
%   E0          : initial Young's modulus [Pa]
%   sigma_max   : saturation stress [Pa]
%
% Outputs
%   sigma : equivalent uniaxial saturation stress [Pa]
%   E_tan : tangent modulus d(sigma)/d(epsilon) [Pa]

epsilon_eff = max(epsilon_eff, 0);

denom = 1 + (E0 .* epsilon_eff) ./ sigma_max;

sigma = (E0 .* epsilon_eff) ./ denom;
E_tan = E0 ./ (denom.^2);

end
