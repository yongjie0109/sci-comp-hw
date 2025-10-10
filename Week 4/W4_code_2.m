clc; clear; format long;

g = @(t) arrayfun(@(z) (z==0).*0 + (z~=0).* (4*z*log(z)/(1+25*z^4)), t);

[I,errbnd] = quadgk(g, 0, 1, 'AbsTol',1e-12, 'RelTol',1e-12);

fprintf('I ≈ %.15f\n', I);
fprintf('Estimated error bound ≤ %.3e\n', errbnd);
