function runge_counterexample()
    clc, clear, close all;
    f = @(x) 1 ./ (1 + x.^2);
    xx = linspace(-5, 5, 2000);
    fx = f(xx);

    ns = [10, 12, 16];

    figure('Color','w'); hold on; box on; grid on;
    plot(xx, fx, 'k-', 'LineWidth', 2);
    title('Runge’s Counterexample on [-5,5]');
    xlabel('x'); ylabel('y');
    legendEntries = {'f(x) = 1/(1+x^2)'};

    cmap = lines(numel(ns));
    for i = 1:numel(ns)
        n  = ns(i);
        xk = linspace(-5, 5, n+1);
        yk = f(xk);
        px = lagrange_eval(xk, yk, xx);

        plot(xx, px, 'Color', cmap(i,:), 'LineWidth', 1.2);
        plot(xk, yk, 'o', 'MarkerSize', 4, 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor', cmap(i,:),'HandleVisibility','off');
        legendEntries{end+1} = sprintf('degree %d', n);
    end

    legend(legendEntries, 'Location','best');
    ylim([-0.5, 1.2]);
end

%------------------ functions -------------------------
function p = lagrange_eval(xk, yk, xeval)
    m = numel(xeval);
    n = numel(xk);
    p = zeros(size(xeval));

    hit_tol = 1e-14;

    for t = 1:m
        x = xeval(t);

        idx = find(abs(x - xk) < hit_tol, 1);
        if ~isempty(idx)
            p(t) = yk(idx);
            continue;
        end

        s = 0;
        for j = 1:n
            Lj_num = 1;
            Lj_den = 1;
            xj = xk(j);
            for k = 1:n
                if k == j, continue; end
                Lj_num = Lj_num * (x - xk(k));
                Lj_den = Lj_den * (xj - xk(k));
            end
            Lj = Lj_num / Lj_den;
            s = s + yk(j) * Lj;
        end
        p(t) = s;
    end
end