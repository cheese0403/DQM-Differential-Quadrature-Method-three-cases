%% Case 01: DQ solution of a simply supported beam under uniform load
% Chen Yanru - Differential Quadrature Method
% Dimensionless governing equation: （无量纲化的微分方程）
%     d^4 W / d xi^4 = 1, 0 <= xi <= 1
% Boundary conditions: （边界条件）
%     W(0)=0, W(1)=0, W''(0)=0, W''(1)=0
% Exact solution: （精确解）
%     W = (xi - 2*xi^3 + xi^4)/24

clear; clc; close all;

% 定义多个节点数，用于观察 DQ 法计算精度随节点数变化的情况。
% 用途：做和精确解的误差对比表。
N_list = [7, 9, 11, 13, 15];

% 画图时使用 N=13 个节点的结果。
N_plot = 13;

% 提前创建结果表。每一行对应一个节点数 N。
% 第1列：N；第2列：最大绝对误差；第3列：跨中挠度 W(0.5)；第4列：跨中相对误差。
result = zeros(numel(N_list), 4);

% 保存不同节点数下的 DQ 计算结果。
% 第1列：xi，节点坐标；第2列：W，DQ 算出来的挠度。
dq_data = cell(numel(N_list), 2);

%% 1. DQ法求解
% 逻辑：对 N_list 里的每一个节点数都算一遍，并把结果存起来；
% 然后单独用 N_plot=13 再算一遍，作为后面画图用的数据。
for k = 1:numel(N_list)
    N = N_list(k);

    % 用 solve_static_beam 函数进行 DQ 求解，返回节点位置 xi 和节点挠度 W。
    [xi, W] = solve_static_beam(N);

    % 将每次计算得到的 xi 和 W 存起来，后面用于误差计算。
    dq_data{k, 1} = xi;
    dq_data{k, 2} = W;
end

% 再计算一次 N_plot 对应的结果，用于画图。
[xi, W] = solve_static_beam(N_plot);

%% 2. 精确解及误差对比
% 逻辑：
% ① 对每个节点数 N，计算精确解；
% ② 把 DQ 解和精确解比较，计算误差；
% ③ 把误差结果写入 results.txt。
for k = 1:numel(N_list)
    % 从前面保存的 dq_data 里取出第 k 次计算的节点坐标。
    xi_k = dq_data{k, 1};

    % 取出第 k 次计算的 DQ 数值挠度。
    W_k = dq_data{k, 2};

    % 在相同节点上计算精确解。
    W_exact_k = exact_solution(xi_k);

    % 计算最大绝对误差：|DQ 数值解 - 精确解| 的最大值。
    max_error = max(abs(W_k - W_exact_k));

    % 求跨中挠度。xi=0.5 不一定正好是一个节点，所以用 interp1 插值估算。
    % linear 表示线性插值。
    W_mid = interp1(xi_k, W_k, 0.5, 'linear');

    % 计算跨中相对误差。简支梁均布载荷无量纲最大挠度精确值为 Wmax = 5/384。
    rel_error = abs(W_mid - 5/384)/(5/384);

    % 把当前 N 的结果存进 result 表格。
    result(k, :) = [N_list(k), max_error, W_mid, rel_error];
end

% W_exact 是 N_plot 节点上的精确解，用于误差分布图。
W_exact = exact_solution(xi);

% 生成更密的点，从 0 到 1 一共取 401 个，用于画一条平滑的精确解曲线。
xi_fine = linspace(0, 1, 401)';

% 在这 401 个密集点上计算精确解。
W_fine = exact_solution(xi_fine);

% 将误差结果写入 results.txt 文件。
fid = fopen('dq_case01_static_beam_annotated_results.txt', 'w');

% 往文件中写标题、控制方程和精确最大挠度。
fprintf(fid, 'Case 01: DQ solution of a simply supported beam under uniform load\n');
fprintf(fid, 'Governing equation: d^4 W / d xi^4 = 1\n');
fprintf(fid, 'Exact maximum deflection: Wmax = 5/384 = %.12f\n\n', 5/384);

% 写表头。
fprintf(fid, '%6s %18s %18s %18s\n', 'N', 'max abs error', 'W(0.5) DQ', 'mid rel error');

% 循环输出 result 表格的每一行。
for k = 1:size(result, 1)
    % result(k,1)：节点数 N；result(k,2)：最大绝对误差；
    % result(k,3)：跨中挠度 W(0.5)；result(k,4)：跨中相对误差。
    fprintf(fid, '%6d %18.6e %18.12f %18.6e\n', result(k, 1), result(k, 2), result(k, 3), result(k, 4));
end

fclose(fid);

%% 3. 画图
% 图1：DQ 解和精确解的挠度对比图。黑色线是精确解，红色点是 DQ 解，用来看两者是否重合。
figure('Color', 'w');
% 精确解 W_exact 随 xi 的变化曲线
plot(xi_fine, W_fine, 'k-', 'LineWidth', 2.0); hold on;
% 画DQ 数值解
plot(xi, W, 'ro', 'MarkerSize', 6, 'LineWidth', 1.4);
grid on; box on;
% 设置坐标轴名称
xlabel('\xi=x/L'); ylabel('W=EIw/(qL^4)');
% 添加图例
legend('Exact solution', 'DQ nodes', 'Location', 'northwest');
% 设置标题
title('Simply supported beam under uniform load');
% 保存成图片文件
saveas(gcf, 'dq_case01_static_beam_annotated_deflection_comparison.png');

% 图2：节点误差分布图。蓝色方块表示每个节点处的误差，用来看 DQ 解和精确解差多少
figure('Color', 'w');
plot(xi, abs(W - W_exact), 'bs-', 'MarkerSize', 6, 'LineWidth', 1.4);
grid on; box on;
xlabel('\xi=x/L'); ylabel('|W_{DQ}-W_{exact}|');
title('Node error distribution');
saveas(gcf, 'dq_case01_static_beam_annotated_error_distribution.png');

%% 下面是上面 1、2、3 步调用的函数

%% 1. DQ法求解函数：求节点 xi 和 DQ 数值挠度 W
% DQ法求解简支梁均布载荷静挠度。
% 输入：N - x方向 DQ 节点数量。
% 输出：xi - 切比雪夫洛巴托节点坐标；W - DQ 数值挠度。
% 给定节点数 N，用 DQ 法求出简支梁在均布载荷下的无量纲挠度 W。
% 即把微分方程 W'''' = 1，加上边界条件 W(0)=0, W(1)=0, W''(0)=0, W''(1)=0，
% 转化成 MATLAB 能解的矩阵方程 A * W = b，然后用 W = A\b 求出来。
% 定义函数solve_static_beam
function [xi, W] = solve_static_beam(N)
    % xi_i = [1 - cos(pi*(i-1)/(N-1))]/2, i = 1,2,...,N。这是 Chebyshev-Gauss-Lobatto 节点公式。
    % 取点不是等间距的，而是两端更密、中间较稀，这样对梁的高阶导数问题更稳定。
    xi = (1 - cos(pi*(0:N-1)'/(N-1)))/2;
		% 构造DQ导数矩阵
    D1 = dq_first_derivative_matrix(xi); % 一阶导数矩阵W' ≈ D1 * W
    D2 = D1*D1; % 二阶导数矩阵 W'' ≈ D1*(D1*W) = D1*D1*W
    D4 = D2*D2; % 四阶导数矩阵 W'''' ≈ D4 * W = D2*D2*W

    % 建立矩阵方程 D4 * W = 1。
    A = D4;
    b = ones(N, 1); % N 行 1 列的全 1 向量

    % 边界 W(0)=0。
    A(1, :) = 0; A(1, 1) = 1; b(1) = 0;
    % 边界 W''(0)=0。
    A(2, :) = D2(1, :); b(2) = 0;
    % 边界 W''(1)=0。
    A(N-1, :) = D2(N, :); b(N-1) = 0;
    % 边界 W(1)=0。
    A(N, :) = 0; A(N, N) = 1; b(N) = 0;

    W = A\b;
end

%% 1. DQ法中的权系数矩阵构造函数：根据节点 x，利用 Lagrange 插值思想，构造 DQ 一阶导数矩阵 D。
% 定义函数dq_first_derivative_matrix，x是节点坐标，D是一阶导数矩阵
function D = dq_first_derivative_matrix(x)
		% 计算节点数
    N = numel(x);
    % 创建一个 N × N 的全零矩阵
    D = zeros(N, N);
    % 创建一个 N × 1 的全 1 向量，c 是辅助变量，用来计算 Lagrange 插值中的乘积项
    c = ones(N, 1);
    
		% 第一组循环计算每个节点对应的 c(i)c_i = Π (x_i - x_j), j ≠ i 把第 i 个节点和其他所有节点的坐标差，全部乘起来
		% eg：x1, x2, x3 c1 = (x1-x2)(x1-x3)；c2 = (x2-x1)(x2-x3)；c3 = (x3-x1)(x3-x2)
    for i = 1:N
        for j = 1:N
            if i ~= j
                c(i) = c(i)*(x(i) - x(j));
            end
        end
    end
		%第二组循环填充导数矩阵D
    for i = 1:N
        for j = 1:N
            if i ~= j
            		% 对应 DQ 一阶导数权系数：a_ij = c_i / [c_j * (x_i - x_j)], i ≠ j
            		% 也就是说第 i 个节点处的一阶导数中，第 j 个节点函数值的权重是多少。
                % 意义是W'(x_i) ≈ D(i,1)*W1 + D(i,2)*W2 + ... + D(i,N)*WN，是D 的每一行对应一个节点处的一阶导数公式。
                D(i, j) = c(i)/(c(j)*(x(i) - x(j)));
            end
        end
        % 计算对角线元素 a_ii = -Σ a_ij, j ≠ i 对角线元素 = 这一行其他元素之和的相反数
        D(i, i) = -sum(D(i, [1:i-1, i+1:N]));
    end
end

%% 2. 精确解函数：用于和 DQ 数值解作对比
% 定义函数exact_solution
function W = exact_solution(xi)
		% 简支梁均布载荷的无量纲解析解
    W = (xi - 2*xi.^3 + xi.^4)/24;
end
