clc; clear; format long;

f = @(x) 1 ./ (1 + 25*x.^2);

% 先找 L，使得尾積分誤差 < 1e-10
target_error = 1e-10;
E = @(L) (1/5) * (pi/2 - atan(5*L));

L = fzero(@(L) E(L) - target_error, 1); % 用 fzero 求 L

fprintf('Truncation limit L = %.10f\n', L);
fprintf('Tail error ≈ %.3e\n', E(L));

% 在 [0, L] 上數值積分
I_num = integral(f, 0, L, 'AbsTol',1e-12, 'RelTol',1e-12);

% 加上尾巴誤差 (理論上可忽略)
I_total = I_num + E(L);

fprintf('Numerical approximation = %.12f\n', I_total);
fprintf('Exact value = %.12f\n', pi/10);
fprintf('Absolute error = %.3e\n', abs(I_total - pi/10));
