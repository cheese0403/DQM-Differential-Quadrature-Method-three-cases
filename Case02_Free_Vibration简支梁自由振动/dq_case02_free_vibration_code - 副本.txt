%% Case 02 简介注释版：简支梁自由振动
% 力学模型：Euler-Bernoulli 简支梁自由振动。
% 无量纲特征方程：Phi'''' = beta^4*Phi。
% 边界条件：Phi(0)=0, Phi(1)=0, Phi''(0)=0, Phi''(1)=0。
% 精确解：beta_n = n*pi，Phi_n = sin(n*pi*xi)。

clear; clc; close all;

N = 19;
mode_count = 5;
plot_mode_count = 3;
result = zeros(mode_count, 4);

%% 1. DQ法求解
[xi, beta, modes] = solve_free_vibration(N);

%% 2. 精确解及误差对比
for n = 1:mode_count
    beta_exact = n*pi;
    rel_error = abs(beta(n) - beta_exact)/beta_exact;
    result(n, :) = [n, beta(n), beta_exact, rel_error];
end

xi_fine = linspace(0, 1, 401)';

fid = fopen('dq_case02_free_vibration_code_results.txt', 'w');
fprintf(fid, 'Case 02: simply supported beam free vibration\n');
fprintf(fid, 'Equation: d^4 Phi / d xi^4 = beta^4*Phi, exact beta_n = n*pi\n\n');
fprintf(fid, '%6s %18s %18s %18s\n', 'mode', 'beta DQ', 'beta exact', 'rel error');
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
xlabel('\xi=x/L'); ylabel('normalized mode shape');
legend('Exact mode 1', 'DQ mode 1', 'Exact mode 2', 'DQ mode 2', 'Exact mode 3', 'DQ mode 3', 'Location', 'best');
title('Simply supported beam free vibration modes');
saveas(gcf, 'dq_case02_free_vibration_code_mode_shape_comparison.png');

figure('Color', 'w');
semilogy(result(:, 1), result(:, 4), 'ms-', 'MarkerSize', 6, 'LineWidth', 1.4);
grid on; box on;
xlabel('mode number'); ylabel('relative error of beta');
title('Frequency parameter error');
saveas(gcf, 'dq_case02_free_vibration_code_frequency_error.png');

%% 函数区
function [xi, beta, modes] = solve_free_vibration(N)
    % DQ 离散：Phi'''' ≈ D4*Phi
    xi = (1 - cos(pi*(0:N-1)'/(N-1)))/2;
    D1 = dq_first_derivative_matrix(xi);
    D2 = D1*D1;
    D4 = D2*D2;

    % 广义特征值问题：A*Phi = lambda*B*Phi，其中 lambda = beta^4
    A = D4;
    B = eye(N);
    B([1, 2, N-1, N], :) = 0;

    % 简支边界：振型位移为零，弯矩 Phi'' 为零
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

    beta = lambda.^(1/4);
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
