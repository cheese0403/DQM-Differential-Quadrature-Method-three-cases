%% Case 02: DQ free vibration of a simply supported Euler-Bernoulli beam
% Chen Yanru - Differential Quadrature Method
% Dimensionless eigenvalue equation: 无量纲自由振动特征方程
%     d^4 Phi / d xi^4 = beta^4 Phi, 0 <= xi <= 1 % Phi 是振型，beta 是频率参数
% Boundary conditions: 边界条件
%     Phi(0)=0, Phi(1)=0, Phi''(0)=0, Phi''(1)=0
% Exact solution: 精确解
%     beta_n = n*pi, Phi_n = sin(n*pi*xi)

clear; clc; close all;

% DQ 节点数量。自由振动是特征值问题，节点数稍多一些结果更稳定。
N = 19;

% 对比前 5 阶频率参数。
mode_count = 5;

% 画图时显示前 3 阶振型。
plot_mode_count = 3;

% result 每一行对应一个振动阶数 n。
% 第1列：阶数 n；第2列：DQ 频率参数 beta；第3列：精确 beta；第4列：相对误差。
result = zeros(mode_count, 4);

%% 1. DQ法求解
% 用 DQ 法求解简支梁自由振动问题，得到节点 xi、频率参数 beta 和振型 modes。
[xi, beta, modes] = solve_free_vibration(N);

%% 2. 精确解及误差对比
% 逻辑：
% ① 对每一阶振型 n，计算精确频率参数 beta_n = n*pi；
% ② 把 DQ 求出的 beta 和精确 beta 比较，计算相对误差；
% ③ 把误差结果写入 results.txt。
for n = 1:mode_count
    % 第 n 阶精确频率参数。
    beta_exact = exact_frequency_parameter(n);

    % 相对误差 = |DQ 解 - 精确解| / 精确解。
    rel_error = abs(beta(n) - beta_exact)/beta_exact;

    % 把当前阶数的结果存进 result 表格。
    result(n, :) = [n, beta(n), beta_exact, rel_error];
end

% 生成更密的点，从 0 到 1 一共取 401 个，用于画平滑的精确振型曲线。
xi_fine = linspace(0, 1, 401)';

% 将误差结果写入 results.txt 文件。
fid = fopen('dq_case02_free_vibration_annotated_results.txt', 'w');

% 往文件中写标题、特征方程和精确解公式。
fprintf(fid, 'Case 02: DQ free vibration of a simply supported Euler-Bernoulli beam\n');
fprintf(fid, 'Eigenvalue equation: D4*Phi = beta^4*Phi\n');
fprintf(fid, 'Exact frequency parameter: beta_n = n*pi\n\n');

% 写表头。
fprintf(fid, '%6s %18s %18s %18s\n', 'mode', 'beta DQ', 'beta exact', 'rel error');

% 循环输出 result 表格的每一行。
for n = 1:size(result, 1)
    % result(n,1)：阶数；result(n,2)：DQ beta；result(n,3)：精确 beta；result(n,4)：相对误差。
    fprintf(fid, '%6d %18.10f %18.10f %18.6e\n', result(n, 1), result(n, 2), result(n, 3), result(n, 4));
end

fclose(fid);

%% 3. 画图
% 图1：DQ 振型和精确振型对比图。曲线是精确解，圆点是 DQ 解。
figure('Color', 'w');
for n = 1:plot_mode_count
    % 第 n 阶精确振型。
    exact_mode = exact_mode_shape(xi_fine, n);

    % 第 n 阶 DQ 振型。
    mode_dq = modes(:, n);

    % 振型可以整体乘以 -1 仍然等价。这里把 DQ 振型方向调到和精确振型一致，方便看图比较。
    mode_dq = align_mode_sign(xi, mode_dq, n);

    plot(xi_fine, exact_mode, 'LineWidth', 1.8); hold on;
    plot(xi, mode_dq, 'o', 'MarkerSize', 5, 'LineWidth', 1.2);
end
grid on; box on;
xlabel('\xi=x/L'); ylabel('normalized mode shape');
legend('Exact mode 1', 'DQ mode 1', 'Exact mode 2', 'DQ mode 2', 'Exact mode 3', 'DQ mode 3', 'Location', 'best');
title('Simply supported beam free vibration modes');
saveas(gcf, 'dq_case02_free_vibration_annotated_mode_shape_comparison.png');

% 图2：频率参数对比图。黑色线是精确频率参数，紫色圆点是 DQ 求出的频率参数。
% 自由振动问题最主要看的不是弯矩剪力，而是每一阶的频率参数 beta 和对应振型。
figure('Color', 'w');
plot(result(:, 1), result(:, 3), 'k-', 'LineWidth', 2.0); hold on;
plot(result(:, 1), result(:, 2), 'mo', 'MarkerSize', 6, 'LineWidth', 1.4);
grid on; box on;
xlabel('mode number'); ylabel('frequency parameter beta');
legend('Exact beta', 'DQ beta', 'Location', 'northwest');
title('Frequency parameter comparison');
saveas(gcf, 'dq_case02_free_vibration_annotated_frequency_comparison.png');

%% 下面是上面 1、2、3 步调用的函数

%% 1. DQ法求解函数：求节点 xi、频率参数 beta 和振型 modes
% DQ法求解简支梁自由振动问题。
% 输入：N - x方向 DQ 节点数量。
% 输出：xi - 切比雪夫洛巴托节点坐标；beta - 频率参数；modes - 振型矩阵。
% 把微分方程 Phi'''' = beta^4*Phi 转化成广义特征值问题 A*Phi = lambda*B*Phi。
% 其中 lambda = beta^4，所以求出 lambda 后再开四次方得到 beta。
function [xi, beta, modes] = solve_free_vibration(N)
    % Chebyshev-Gauss-Lobatto 节点公式。
    % 节点在两端更密、中间较稀，对高阶导数计算更稳定。
    xi = (1 - cos(pi*(0:N-1)'/(N-1)))/2;

    % 构造 DQ 导数矩阵。
    D1 = dq_first_derivative_matrix(xi); % 一阶导数矩阵 Phi' ≈ D1 * Phi
    D2 = D1*D1; % 二阶导数矩阵 Phi'' ≈ D2 * Phi
    D4 = D2*D2; % 四阶导数矩阵 Phi'''' ≈ D4 * Phi

    % 建立广义特征值方程 A*Phi = lambda*B*Phi。
    A = D4;
    B = eye(N); % D4 * Phi = beta^4 * Phi，其中lambda = beta^4

    % 边界条件所在的行不参与右端特征值项，所以 B 对应行置零。
    B([1, 2, N-1, N], :) = 0;

    % 边界 Phi(0)=0。
    A(1, :) = 0; A(1, 1) = 1;
    % 边界 Phi''(0)=0。
    A(2, :) = D2(1, :);
    % 边界 Phi''(1)=0。
    A(N-1, :) = D2(N, :);
    % 边界 Phi(1)=0。
    A(N, :) = 0; A(N, N) = 1;

    % 求广义特征值和特征向量。E 里面是特征值 lambda，V 里面是特征向量，也就是振型
    [V, E] = eig(A, B);
    lambda = diag(E);

    % 只保留有限、实数、正的特征值。
    valid = isfinite(lambda) & abs(imag(lambda)) < 1e-7 & real(lambda) > 1e-8;
    lambda = real(lambda(valid));
    V = real(V(:, valid));

    % 按特征值（频率）从小到大排序，前面对应低阶振型。
    [lambda, idx] = sort(lambda);
    V = V(:, idx);

    % lambda = beta^4，所以 beta = lambda^(1/4)。
    beta = lambda.^(1/4);
    modes = V;

    % 把振型归一化，每一阶振型都除以最大绝对值，使最大幅值为 1，便于画图比较。
    for k = 1:size(modes, 2)
        modes(:, k) = modes(:, k)/max(abs(modes(:, k)));
    end
end

%% 1. DQ法中的权系数矩阵构造函数：根据节点 x 得到一阶导数矩阵 D
function D = dq_first_derivative_matrix(x)
    % 计算节点数。
    N = numel(x);

    % 创建一个 N × N 的全零矩阵。
    D = zeros(N, N);

    % c 是辅助变量，用来计算 Lagrange 插值中的乘积项。
    c = ones(N, 1);

    % 第一组循环计算每个节点对应的 c(i) = Π (x_i - x_j), j ≠ i。
    for i = 1:N
        for j = 1:N
            if i ~= j
                c(i) = c(i)*(x(i) - x(j));
            end
        end
    end

    % 第二组循环填充一阶导数矩阵 D。
    for i = 1:N
        for j = 1:N
            if i ~= j
                % D(i,j) 表示第 i 个节点处一阶导数中，第 j 个节点函数值的权重。
                D(i, j) = c(i)/(c(j)*(x(i) - x(j)));
            end
        end
        % 对角线元素等于这一行其他元素之和的相反数。
        D(i, i) = -sum(D(i, [1:i-1, i+1:N]));
    end
end

%% 2. 精确频率参数函数：用于和 DQ 数值解作对比
function beta_exact = exact_frequency_parameter(n)
    beta_exact = n*pi;
end

%% 2. 精确振型函数：用于画精确振型曲线
function Phi = exact_mode_shape(xi, n)
    Phi = sin(n*pi*xi);
end

%% 3. 振型方向调整函数：让 DQ 振型和精确振型方向一致
function mode_dq = align_mode_sign(xi, mode_dq, n)
    % 找到 DQ 振型绝对值最大的节点。
    [~, idx] = max(abs(mode_dq));

    % 如果该节点处 DQ 振型和精确振型符号相反，就把整个 DQ 振型乘以 -1。
    if mode_dq(idx)*sin(n*pi*xi(idx)) < 0
        mode_dq = -mode_dq;
    end
end
