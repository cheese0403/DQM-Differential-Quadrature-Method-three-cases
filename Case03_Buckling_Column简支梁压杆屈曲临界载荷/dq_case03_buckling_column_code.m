%% Case 03 简介注释版：简支压杆屈曲临界载荷
% 力学模型：两端简支压杆，承受轴向压力 P。
% 无量纲方程：W'''' + Pbar*W'' = 0。
% 广义特征值形式：D4*W = Pbar*(-D2)*W。
% 精确解：Pbar_n = n^2*pi^2，W_n = sin(n*pi*xi)。

clear; clc; close all;

N = 25;
mode_count = 5;
plot_mode_count = 3;
result = zeros(mode_count, 4);

%% 1. DQ法求解
[xi, Pbar, modes] = solve_buckling(N);

%% 2. 精确解及误差对比
for n = 1:mode_count
    Pbar_exact = (n*pi)^2;
    rel_error = abs(Pbar(n) - Pbar_exact)/Pbar_exact;
    result(n, :) = [n, Pbar(n), Pbar_exact, rel_error];
end

xi_fine = linspace(0, 1, 401)';

fid = fopen('dq_case03_buckling_column_code_results.txt', 'w');
fprintf(fid, 'Case 03: simply supported column buckling\n');
fprintf(fid, 'Equation: D4*W = Pbar*(-D2)*W, exact Pbar_n = n^2*pi^2\n\n');
fprintf(fid, '%6s %18s %18s %18s\n', 'mode', 'Pbar DQ', 'Pbar exact', 'rel error');
for n = 1:size(result, 1)
    fprintf(fid, '%6d %18.10f %18.10f %18.6e\n', result(n, 1), result(n, 2), result(n, 3), result(n, 4));
end
fclose(fid);

%% 3. 画图
figure('Color', 'w');
for n = 1:plot_mode_count
    exact_mode = sin(n*pi*xi_fine);
    mode_dq = align_mode_sign(xi, modes(:, n), n);
    plot(xi_fine, exact_mode, 'LineWidth', 1.8); hold on;
    plot(xi, mode_dq, 'o', 'MarkerSize', 5, 'LineWidth', 1.2);
end
grid on; box on;
xlabel('\xi=x/L'); ylabel('normalized buckling mode');
legend('Exact mode 1', 'DQ mode 1', 'Exact mode 2', 'DQ mode 2', 'Exact mode 3', 'DQ mode 3', 'Location', 'best');
title('Simply supported column buckling modes');
saveas(gcf, 'dq_case03_buckling_column_code_buckling_mode_comparison.png');

figure('Color', 'w');
plot(result(:, 1), result(:, 3), 'k-', 'LineWidth', 2.0); hold on;
plot(result(:, 1), result(:, 2), 'go', 'MarkerSize', 6, 'LineWidth', 1.4);
grid on; box on;
xlabel('mode number'); ylabel('critical load parameter Pbar');
legend('Exact Pbar', 'DQ Pbar', 'Location', 'northwest');
title('Buckling critical load comparison');
saveas(gcf, 'dq_case03_buckling_column_code_buckling_load_comparison.png');

%% 函数区
function [xi, Pbar, modes] = solve_buckling(N)
    % DQ 离散：W'' ≈ D2*W，W'''' ≈ D4*W
    xi = (1 - cos(pi*(0:N-1)'/(N-1)))/2;
    D1 = dq_first_derivative_matrix(xi);
    D2 = D1*D1;
    D4 = D2*D2;

    % 广义特征值问题：A*W = Pbar*B*W
    A = D4;
    B = -D2;
    B([1, 2, N-1, N], :) = 0;

    % 简支边界：位移为零，弯矩 W'' 为零
    A(1, :) = 0; A(1, 1) = 1;
    A(2, :) = D2(1, :);
    A(N-1, :) = D2(N, :);
    A(N, :) = 0; A(N, N) = 1;

    [V, E] = eig(A, B);
    lambda = diag(E);
    valid = isfinite(lambda) & abs(imag(lambda)) < 1e-7 & real(lambda) > 1e-8;
    lambda = real(lambda(valid));
    V = real(V(:, valid));
    [lambda, idx] = sort(lambda);
    V = V(:, idx);

    Pbar = lambda;
    modes = normalize_modes(V);
end

function D = dq_first_derivative_matrix(x)
    N = numel(x);
    D = zeros(N, N);
    c = ones(N, 1);

    for i = 1:N
        for j = 1:N
            if i ~= j
                c(i) = c(i)*(x(i) - x(j));
            end
        end
    end

    for i = 1:N
        for j = 1:N
            if i ~= j
                D(i, j) = c(i)/(c(j)*(x(i) - x(j)));
            end
        end
        D(i, i) = -sum(D(i, [1:i-1, i+1:N]));
    end
end

function modes = normalize_modes(modes)
    for k = 1:size(modes, 2)
        modes(:, k) = modes(:, k)/max(abs(modes(:, k)));
    end
end

function mode_dq = align_mode_sign(xi, mode_dq, n)
    [~, idx] = max(abs(mode_dq));
    if mode_dq(idx)*sin(n*pi*xi(idx)) < 0
        mode_dq = -mode_dq;
    end
end
