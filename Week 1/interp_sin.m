% Demo for f(x)=sin(x) on [0,1] with 10,20,30 nodes with different algorithms.
clc, clear, close all;
f = @(x) sin(x);
Ns = [10 20 30];
xx = linspace(0,1,2001)';
fx = f(xx);

algos = {'Newton','Lagrange','ModLagrange','Barycentric'};
maxErr = zeros(numel(Ns), numel(algos));

figure('Name','Interpolants of sin(x) on [0,1]','Color','w');

for i = 1:numel(Ns)
    n = Ns(i);
    xk = linspace(0,1,n)'; 
    yk = f(xk);

    aN = newton_divdiff(xk, yk);

    w = barycentric_weights(xk);

    % Evaluate all methods on grid
    p_newton      = newton_eval(xk, aN, xx);
    p_lagrange    = lagrange_eval(xk, yk, xx);
    p_modlagrange = modlagrange_eval_first_form(xk, yk, xx, w);
    p_bary        = barycentric_eval_second_form(xk, yk, xx, w);

    % Errors
    P = [p_newton, p_lagrange, p_modlagrange, p_bary];
    maxErr(i,:) = max(abs(P - fx), [], 1);

    % Plot for this n
    subplot(3,1,i);
    plot(xx, fx, 'k-', 'LineWidth', 1.5); hold on;
    plot(xx, p_newton,      '-', 'LineWidth', 1.0);
    plot(xx, p_lagrange,    '--', 'LineWidth', 1.0);
    plot(xx, p_modlagrange, '-.', 'LineWidth', 1.0);
    plot(xx, p_bary,        ':', 'LineWidth', 1.5);
    plot(xk, yk, 'ko', 'MarkerFaceColor','w', 'MarkerSize',4);
    grid on;
    title(sprintf('n = %d nodes (sin on [0,1])', n));
    xlabel('x'); ylabel('p_n(x)');
    legend({'f(x)','Newton','Lagrange','Modified Lagrange','Barycentric','nodes'}, ...
           'Location','best');
end

% Print max error table
fprintf('\nMax |error| on xx grid (sin on [0,1])\n');
header = sprintf('%10s %12s %12s %16s %14s\n', 'n','Newton','Lagrange','ModLagrange','Barycentric');
fprintf(header);
for i = 1:numel(Ns)
    fprintf('%10d %12.3e %12.3e %16.3e %14.3e\n', Ns(i), maxErr(i,1), maxErr(i,2), maxErr(i,3), maxErr(i,4));
end
fprintf('\nNotes:\n- Lagrange (naive) and Modified/Barycentric should agree numerically; the latter are far more efficient.\n- Barycentric (second form) is typically the most stable for evaluation.\n');


%------------------------- function -----------------------------
%% Newton
%Divided differences matrix
function a = newton_divdiff(x, y)
    x = x(:); y = y(:);
    n = numel(x);
    DD = zeros(n,n);
    DD(:,1) = y;

    for j = 2:n
        for i = 1:(n-j+1)
            DD(i,j) = (DD(i+1,j-1) - DD(i,j-1)) / (x(i+j-1) - x(i));
        end
    end
    a = DD(1,:)';
end

function p = newton_eval(xnodes, a, xeval)
    xnodes = xnodes(:);
    n = numel(xnodes);
    xeval = xeval(:);
    m = numel(xeval);
    p = a(n) * ones(m,1);
    for k = n-1:-1:1
        p = a(k) + (xeval - xnodes(k)).*p;
    end
end

%% Lagrange
function p = lagrange_eval(xk, yk, xeval)
    xk = xk(:); yk = yk(:); xeval = xeval(:);
    n = numel(xk); m = numel(xeval);
    p = zeros(m,1);
    hit_tol = 1e-14;

    for t = 1:m
        x = xeval(t);
        idx = find(abs(x - xk) < hit_tol, 1);
        if ~isempty(idx)
            p(t) = yk(idx);
            continue;
        end

        s = 0.0;
        for j = 1:n
            Lj_num = 1.0; Lj_den = 1.0;
            xj = xk(j);
            for k = 1:n
                if k == j, continue; end
                Lj_num = Lj_num * (x - xk(k));
                Lj_den = Lj_den * (xj - xk(k));
            end
            s = s + yk(j) * (Lj_num / Lj_den);
        end
        p(t) = s;
    end
end

%% Modifiied Lagrange
function p = modlagrange_eval_first_form(xk, yk, xeval, w)
    xk = xk(:); yk = yk(:); xeval = xeval(:); w = w(:);
    m = numel(xeval);
    n = numel(xk);
    p = zeros(m,1);
    hit_tol = 1e-14;

    for t = 1:m
        x = xeval(t);
        idx = find(abs(x - xk) < hit_tol, 1);
        if ~isempty(idx)
            p(t) = yk(idx);
            continue;
        end

        diff = x - xk;                
        Px = prod(diff);              
        s  = sum( (w .* yk) ./ diff );
        p(t) = Px * s;
    end
end



%% Barycentric Lagrange
function w = barycentric_weights(xk)
% Compute barycentric weights
    xk = xk(:);
    n = numel(xk);
    w = ones(n,1);
    for j = 1:n
        diff = xk(j) - xk;
        diff(j) = 1;   
        w(j) = 1 / prod(diff);
    end
end

function p = barycentric_eval_second_form(xk, yk, xeval, w)
    xk = xk(:); yk = yk(:); xeval = xeval(:); w = w(:);
    m = numel(xeval);
    p = zeros(m,1);
    hit_tol = 1e-14;

    for t = 1:m
        x = xeval(t);
        idx = find(abs(x - xk) < hit_tol, 1);
        if ~isempty(idx)
            p(t) = yk(idx);
            continue;
        end
        diff = x - xk;                  
        num = sum( (w .* yk) ./ diff );
        den = sum(  w        ./ diff );
        p(t) = num / den;
    end
end
