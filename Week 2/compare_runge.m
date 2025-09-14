clc; clear; close all;
% 比較 Runge 函數 f(x)=1/(1+25x^2) 在 [-1,1] 上的插值近似
% 1) Cubic spline: not-a-knot
% 2) Cubic spline: natural
% 3) Cubic spline: clamped
% 4) Chebyshev-I 
% 5) Chebyshev-II 
f = @(x) 1./(1 + 25*x.^2);
fp = @(x) -50*x ./ (1+25*x.^2).^2;   % 給 clamped 用

tol = 1e-10;
Nmax = 20000;
xx = linspace(-1,1,20001).';
fxx = f(xx);

% ---- Spline: not-a-knot ----
[N_nak, err_nak] = search_minN_spline(@(x,y) spline(x,y), f, fxx, xx, tol, Nmax);

% ---- Spline: natural (Natural) ----
[N_nat, err_nat] = search_minN_spline(@(x,y) csape(x,y,'variational'), f, fxx, xx, tol, Nmax);

% ---- Spline: clamped ----
[N_cla, err_cla] = search_minN_spline(@(x,y) csape(x,[fp(x(1)), y, fp(x(end))],'clamped'), f, fxx, xx, tol, Nmax);

% ---- Chebyshev I ----
[N1, e1] = search_minN_cheb(@nodes_cheb1, f, xx, fxx, tol, Nmax);

% ---- Chebyshev II ----
[N2, e2] = search_minN_cheb(@nodes_cheb2, f, xx, fxx, tol, Nmax);

% ---- Output ----
fprintf('Cubic Spline (Not-a-knot) : minimal N = %d, max err = %.3e\n', N_nak, err_nak);
fprintf('Cubic Spline (Natural) : minimal N = %d, max err = %.3e\n', N_nat, err_nat);
fprintf('Cubic Spline (Clamped) : minimal N = %d, max err = %.3e\n', N_cla, err_cla);
fprintf('Chebyshev-I : minimal N = %d, max_err = %.3e\n', N1, e1);
fprintf('Chebyshev-II : minimal N = %d, max_err = %.3e\n', N2, e2);

%% -------------------------Cubic Spline -----------------------------
function [bestN, bestErr] = search_minN_spline(splinefun, f, fxx, xx, tol, Nmax)
    bestN = NaN; 
    bestErr = Inf;
    for N = 2:Nmax
        x  = linspace(-1,1,N+1);
        y  = f(x);
        pp = splinefun(x,y);
        sxx = ppval(pp, xx);
        err = max(abs(sxx - fxx));
        if err < tol
            bestN = N; 
            bestErr = err;
            break
        end
    end
end

%% ------------------------- Chebyshev  -----------------------------
function [bestN, bestErr] = search_minN_cheb(nodefun, f, xx, fxx, tol, Nmax)
    bestN = NaN; 
    bestErr = Inf;
    for N = 2:Nmax
        [x, w] = nodefun(N);
        y = f(x);
        pxx = bary_eval_point(x, y, w, xx);
        err = max(abs(pxx - fxx));
        if err < tol
            bestN = N; 
            bestErr = err;
            break
        end
    end
end

% Chebyshev-I 
function [x, w] = nodes_cheb1(N)
    k = (0:N).';
    theta = (2*k+1)*pi/(2*(N+1));
    x = cos(theta);
    w = (-1).^k .* sin(theta);
end

% Chebyshev-II 
function [x, w] = nodes_cheb2(N)
    j = (0:N).';
    x = cos(pi*j/N);
    c = ones(N+1,1); c(1)=0.5; c(end)=0.5;
    w = ((-1).^j).*c;
end

% Barycentric evaluation
function p = bary_eval_point(x, y, w, xq)
    x = x(:); y = y(:); w = w(:);
    xq = xq(:);
    n = numel(x);
    m = numel(xq);
    p = zeros(m,1);
    for t = 1:m
        xt = xq(t);
        idx = find(abs(xt - x) < 1e-14, 1);
        if ~isempty(idx)
            p(t) = y(idx); continue
        end
        den = 0; 
        num = 0;
        for j = 1:n
            inv = w(j) / (xt - x(j));
            den = den + inv;
            num = num + inv * y(j);
        end
        p(t) = num / den;
    end
end

