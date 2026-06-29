%% Case 03: DQ buckling analysis of a simply supported Euler column
% Chen Yanru - Differential Quadrature Method
% Dimensionless governing equation: 无量纲压杆屈曲控制方程
%     W'''' + Pbar*W'' = 0, 0 <= xi <= 1
% Written as a generalized eigenvalue problem: 写成广义特征值问题
%     D4*W = Pbar*(-D2)*W
% Boundary conditions: 边界条件
%     W(0)=0, W(1)=0, W''(0)=0, W''(1)=0
% Exact solution: 精确解
%     Pbar_n = n^2*pi^2, W_n = sin(n*pi*xi)

clear; clc; close all;

% DQ 节点数量。屈曲也是特征值问题，节点数稍多一些结果更稳定。
N = 25;

% 对比前 5 阶无量纲临界载荷。
mode_count = 5;

% 画图时显示前 3 阶屈曲模态。
plot_mode_count = 3;

% result 每一行对应一个屈曲阶数 n。
% 第1列：阶数 n；第2列：DQ 无量纲临界载荷 Pbar；第3列：精确 Pbar；第4列：相对误差。
result = zeros(mode_count, 4);

%% 1. DQ法求解
% 用 DQ 法求解简支压杆屈曲问题，得到节点 xi、无量纲临界载荷 Pbar 和屈曲模态 modes。
[xi, Pbar, modes] = solve_buckling(N);

%% 2. 精确解及误差对比
% 逻辑：
% ① 对每一阶屈曲模态 n，计算精确临界载荷 Pbar_n = n^2*pi^2；
% ② 把 DQ 求出的 Pbar 和精确 Pbar 比较，计算相对误差；
% ③ 把误差结果写入 results.txt。
for n = 1:mode_count
    % 第 n 阶精确无量纲临界载荷。
    Pbar_exact = exact_buckling_load(n);

    % 相对误差 = |DQ 解 - 精确解| / 精确解。
    rel_error = abs(Pbar(n) - Pbar_exact)/Pbar_exact;

    % 把当前阶数的结果存进 result 表格。
    result(n, :) = [n, Pbar(n), Pbar_exact, rel_error];
end

% 生成更密的点，从 0 到 1 一共取 401 个，用于画平滑的精确屈曲模态曲线。
xi_fine = linspace(0, 1, 401)';

% 将误差结果写入 results.txt 文件。
fid = fopen('dq_case03_buckling_column_annotated_results.txt', 'w');

% 往文件中写标题、特征方程和精确解公式。
fprintf(fid, 'Case 03: DQ buckling analysis of a simply supported Euler column\n');
fprintf(fid, 'Eigenvalue equation: D4*W = Pbar*(-D2)*W\n');
fprintf(fid, 'Exact nondimensional critical load: Pbar_n = n^2*pi^2\n\n');

% 写表头。
fprintf(fid, '%6s %18s %18s %18s\n', 'mode', 'Pbar DQ', 'Pbar exact', 'rel error');

% 循环输出 result 表格的每一行。
for n = 1:size(result, 1)
    % result(n,1)：阶数；result(n,2)：DQ Pbar；result(n,3)：精确 Pbar；result(n,4)：相对误差。
    fprintf(fid, '%6d %18.10f %18.10f %18.6e\n', result(n, 1), result(n, 2), result(n, 3), result(n, 4));
end

fclose(fid);

%% 3. 画图
% 图1：DQ 屈曲模态和精确屈曲模态对比图。曲线是精确解，圆点是 DQ 解。
figure('Color', 'w');
for n = 1:plot_mode_count
    % 第 n 阶精确屈曲模态。
    exact_mode = exact_mode_shape(xi_fine, n);

    % 第 n 阶 DQ 屈曲模态。
    mode_dq = modes(:, n);

    % 模态可以整体乘以 -1 仍然等价。这里把 DQ 模态方向调到和精确模态一致，方便看图比较。
    mode_dq = align_mode_sign(xi, mode_dq, n);

    plot(xi_fine, exact_mode, 'LineWidth', 1.8); hold on;
    plot(xi, mode_dq, 'o', 'MarkerSize', 5, 'LineWidth', 1.2);
end
grid on; box on;
xlabel('\xi=x/L'); ylabel('normalized buckling mode');
legend('Exact mode 1', 'DQ mode 1', 'Exact mode 2', 'DQ mode 2', 'Exact mode 3', 'DQ mode 3', 'Location', 'best');
title('Simply supported column buckling modes');
saveas(gcf, 'dq_case03_buckling_column_annotated_buckling_mode_comparison.png');

% 图2：无量纲临界载荷相对误差图。纵坐标用 semilogy，便于显示很小的误差。
figure('Color', 'w');
semilogy(result(:, 1), result(:, 4), 'gs-', 'MarkerSize', 6, 'LineWidth', 1.4);
grid on; box on;
xlabel('mode number'); ylabel('relative error of Pbar');
title('Buckling load parameter error');
saveas(gcf, 'dq_case03_buckling_column_annotated_buckling_load_error.png');

%% 下面是上面 1、2、3 步调用的函数

%% 1. DQ法求解函数：求节点 xi、临界载荷 Pbar 和屈曲模态 modes
% DQ法求解简支压杆屈曲临界载荷问题。
% 输入：N - x方向 DQ 节点数量。
% 输出：xi - 切比雪夫洛巴托节点坐标；Pbar - 无量纲临界载荷；modes - 屈曲模态矩阵。
% 把微分方程 W'''' + Pbar*W'' = 0 转化成广义特征值问题 A*W = Pbar*B*W。
% 其中 A = D4，B = -D2。
function [xi, Pbar, modes] = solve_buckling(N)
    % Chebyshev-Gauss-Lobatto 节点公式。
    % 节点在两端更密、中间较稀，对高阶导数计算更稳定。
    xi = (1 - cos(pi*(0:N-1)'/(N-1)))/2;

    % 构造 DQ 导数矩阵。
    D1 = dq_first_derivative_matrix(xi); % 一阶导数矩阵 W' ≈ D1 * W
    D2 = D1*D1; % 二阶导数矩阵 W'' ≈ D2 * W
    D4 = D2*D2; % 四阶导数矩阵 W'''' ≈ D4 * W

    % 建立广义特征值方程 A*W = Pbar*B*W。
    A = D4;
    B = -D2; % W'''' + Pbar W'' = 0 → W'''' = Pbar(-W'') → D4*W = Pbar*(-D2)*W

    % 边界条件所在的行不参与右端特征值项，所以 B 对应行置零。
    B([1, 2, N-1, N], :) = 0;

    % 边界 W(0)=0。
    A(1, :) = 0; A(1, 1) = 1;
    % 边界 W''(0)=0。
    A(2, :) = D2(1, :);
    % 边界 W''(1)=0。
    A(N-1, :) = D2(N, :);
    % 边界 W(1)=0。
    A(N, :) = 0; A(N, N) = 1;

    % 求广义特征值和特征向量。
    [V, E] = eig(A, B); % lambda 是特征值，V 是特征向量
    lambda = diag(E);

    % 只保留有限、实数、正的特征值。
    valid = isfinite(lambda) & abs(imag(lambda)) < 1e-7 & real(lambda) > 1e-8;
    lambda = real(lambda(valid));
    V = real(V(:, valid));

    % 按特征值（临界载荷）从小到大排序，前面对应低阶屈曲模态。
    [lambda, idx] = sort(lambda);
    V = V(:, idx);

    % 这里的特征值 lambda 就是无量纲临界载荷 Pbar。
    Pbar = lambda;
    modes = V;

    % 模态归一化，每一阶模态都除以最大绝对值，使最大幅值为 1，便于画图比较。
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

%% 2. 精确临界载荷函数：用于和 DQ 数值解作对比
function Pbar_exact = exact_buckling_load(n)
    Pbar_exact = (n*pi)^2;
end

%% 2. 精确屈曲模态函数：用于画精确模态曲线
function W = exact_mode_shape(xi, n)
    W = sin(n*pi*xi);
end

%% 3. 模态方向调整函数：让 DQ 模态和精确模态方向一致
function mode_dq = align_mode_sign(xi, mode_dq, n)
    % 找到 DQ 模态绝对值最大的节点。
    [~, idx] = max(abs(mode_dq));

    % 如果该节点处 DQ 模态和精确模态符号相反，就把整个 DQ 模态乘以 -1。
    if mode_dq(idx)*sin(n*pi*xi(idx)) < 0
        mode_dq = -mode_dq;
    end
end
