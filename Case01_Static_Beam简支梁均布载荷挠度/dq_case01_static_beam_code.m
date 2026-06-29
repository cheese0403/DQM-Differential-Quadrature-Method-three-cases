%% Case 01 简介注释版：简支梁均布载荷静挠度
% 力学模型：Euler-Bernoulli 简支梁，承受均布载荷 q。
% 无量纲方程：W'''' = 1。
% 边界条件：W(0)=0, W(1)=0, W''(0)=0, W''(1)=0。
% 精确解：W = (xi - 2*xi^3 + xi^4)/24。

clear; clc; close all;

N_list = [7, 9, 11, 13, 15]; % 用于误差对比的节点数
N_plot = 13;                 % 用于画图的节点数
result = zeros(numel(N_list), 4);
dq_data = cell(numel(N_list), 2);

%% 1. DQ法求解
for k = 1:numel(N_list)
    [xi_k, W_k] = solve_static_beam(N_list(k));
    dq_data{k, 1} = xi_k;
    dq_data{k, 2} = W_k;
end

[xi, W] = solve_static_beam(N_plot);

%% 2. 精确解及误差对比
for k = 1:numel(N_list)
    xi_k = dq_data{k, 1};
    W_k = dq_data{k, 2};
    W_exact_k = exact_solution(xi_k);

    max_error = max(abs(W_k - W_exact_k));
    W_mid = interp1(xi_k, W_k, 0.5, 'linear');
    rel_error = abs(W_mid - 5/384)/(5/384);
    result(k, :) = [N_list(k), max_error, W_mid, rel_error];
end

W_exact = exact_solution(xi);
xi_fine = linspace(0, 1, 401)';
W_fine = exact_solution(xi_fine);

fid = fopen('dq_case01_static_beam_code_results.txt', 'w');
fprintf(fid, 'Case 01: simply supported beam under uniform load\n');
fprintf(fid, 'Equation: d^4 W / d xi^4 = 1, exact Wmax = 5/384\n\n');
fprintf(fid, '%6s %18s %18s %18s\n', 'N', 'max abs error', 'W(0.5) DQ', 'mid rel error');
for k = 1:size(result, 1)
    fprintf(fid, '%6d %18.6e %18.12f %18.6e\n', result(k, 1), result(k, 2), result(k, 3), result(k, 4));
end
fclose(fid);

%% 3. 画图
figure('Color', 'w');
plot(xi_fine, W_fine, 'k-', 'LineWidth', 2.0); hold on;
plot(xi, W, 'ro', 'MarkerSize', 6, 'LineWidth', 1.4);
grid on; box on;
xlabel('\xi=x/L'); ylabel('W=EIw/(qL^4)');
legend('Exact solution', 'DQ nodes', 'Location', 'northwest');
title('Simply supported beam under uniform load');
saveas(gcf, 'dq_case01_static_beam_code_deflection_comparison.png');

figure('Color', 'w');
plot(xi, abs(W - W_exact), 'bs-', 'MarkerSize', 6, 'LineWidth', 1.4);
grid on; box on;
xlabel('\xi=x/L'); ylabel('|W_{DQ}-W_{exact}|');
title('Node error distribution');
saveas(gcf, 'dq_case01_static_beam_code_error_distribution.png');

%% 函数区
function [xi, W] = solve_static_beam(N)
    % DQ 节点和导数矩阵
    xi = (1 - cos(pi*(0:N-1)'/(N-1)))/2;
    D1 = dq_first_derivative_matrix(xi);
    D2 = D1*D1;
    D4 = D2*D2;

    % 矩阵方程 D4*W = 1
    A = D4;
    b = ones(N, 1);

    % 简支边界：位移为零，弯矩 W'' 为零
    A(1, :) = 0; A(1, 1) = 1; b(1) = 0;
    A(2, :) = D2(1, :); b(2) = 0;
    A(N-1, :) = D2(N, :); b(N-1) = 0;
    A(N, :) = 0; A(N, N) = 1; b(N) = 0;

    W = A\b;
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

function W = exact_solution(xi)
    W = (xi - 2*xi.^3 + xi.^4)/24;
end
