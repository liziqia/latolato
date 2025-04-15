function lato_lato()
    % 参数设置
    L = 1;          % 绳长一半 (m)
    g = 9.81;       % 重力加速度 (m/s²)
    A = 0.1;        % 支点振幅 (m)
    omega = 2*sqrt(g/L); % 驱动频率 (rad/s)
    cr = 0.8;       % 恢复系数（0:完全非弹性, 1:完全弹性）
    theta0 = 0.2;   % 初始摆角 (rad)
    dtheta0 = 0;    % 初始角速度 (rad/s)
    tspan = 0:0.01:100; % 仿真时间 (s)
    
    % 初始状态 [角度，角速度]
    y0 = [theta0; dtheta0];
    % 添加小的随机扰动
    perturbation = 0.01*randn(2,1);
    y0 = y0 + perturbation;
    
    % 求解ODE（包含碰撞事件）
    options = odeset('Events', @(t,y) collision_event(t, y, cr), 'RelTol', 1e-6);
    current_y = y0;
    current_tspan = [tspan(1) tspan(end)];
    t_all = [];
    y_all = [];

    while current_tspan(1) < current_tspan(2)
        [t, y, te, ye, ie] = ode45(@(t,y) lato_ode(t, y, L, g, A, omega), ...
                                    current_tspan, current_y, options);
        
        t_all = [t_all; t];
        y_all = [y_all; y];
        
        if ~isempty(te)
            % 应用恢复系数更新角速度
            current_y = [0; -cr * y(end,2)];
            current_tspan = [te, current_tspan(2)];
        else
            break;
        end
    end
    
    
    [t_unique, idx] = unique(t_all);
    y_interp = interp1(t_unique, y_all(idx,:), tspan);
    
    % 动画可视化
    animate_system(tspan, y_interp, L, A, omega);

    % 相空间图
    figure;
    plot(y_interp(:,1), y_interp(:,2));
    xlabel('\theta (rad)');
    ylabel('d\theta/dt (rad/s)');
    title('相空间轨迹');

    % 计算最大李雅普诺夫指数
    % lambda = lyapunov_exponent(tspan, y_interp);
    % fprintf('Maximum Lyapunov Exponent: %.4f\n', lambda);
end
%% 
% 运动方程（含支点驱动）
function dydt = lato_ode(t, y, L, g, A, omega)
    theta = y(1);
    dtheta = y(2);
    ypp = -A * omega^2 * sin(omega * t); % 支点加速度
    d2theta = -((g + ypp)/L) * sin(theta);
    dydt = [dtheta; d2theta];
end

% 碰撞事件检测（含恢复系数）
function [value, isterminal, direction] = collision_event(~, y, cr)
    value = y(1);       % 检测θ=0
    isterminal = 1;     % 停止积分
    direction = -1;     % 仅检测θ从正方向过零点
end
%%
% 动画生成函数
function animate_system(t, y, L, A, omega)
    figure;
    axis equal;
    xlim([-2*L 2*L]);
    ylim([-2*L 2*L]);
    hold on;
    
    pivot = plot(0, A*sin(omega*t(1)), 'ko', 'MarkerSize', 10);
    balls = plot([-L*sin(y(1,1)), L*sin(y(1,1))], ...
                [A*sin(omega*t(1))-L*cos(y(1,1)), A*sin(omega*t(1))-L*cos(y(1,1))], ...
                'ro', 'MarkerSize', 10);
    
    frame_skip = 10;
    for k = 1:frame_skip:length(t)
        theta = y(k,1);
        y_p = A*sin(omega*t(k));
        x_balls = [-L*sin(theta), L*sin(theta)];
        y_balls = y_p - L*cos(theta);
        
        set(pivot, 'YData', y_p);
        set(balls, 'XData', x_balls, 'YData', [y_balls, y_balls]);
        drawnow limitrate;
        pause(0.01);
    end
end
%%
%{
function lambda = lyapunov_exponent(t, y)
    % 确保时间步长均匀
    dt = mean(diff(t));
    
    % 选择参考轨迹
    ref_traj = y(:,1:2)';  % 转置为2×N矩阵
    
    % 参数设置
    min_norm = 1e-6;  % 最小扰动范数
    max_norm = 1e-2;  % 最大扰动范数
    n = length(t);
    
    % 初始化变量
    lambda_sum = 0;
    total_time = 0;
    ref_idx = 1;
    d1 = ref_traj(:,1) + min_norm * randn(2,1); % 初始扰动
    
    for i = 2:n
        % 计算当前扰动距离
        delta = norm(ref_traj(:,i) - d1);
        
        if delta > max_norm
            % 计算实际时间间隔
            delta_t = (i - ref_idx) * dt;
            
            % 计算增长率并累加
            growth = delta / norm(ref_traj(:,ref_idx) - d1);
            lambda_sum = lambda_sum + log(growth);
            total_time = total_time + delta_t;
            
            % 重新缩放扰动向量
            d1 = ref_traj(:,i) + (d1 - ref_traj(:,i)) * min_norm / delta;
            ref_idx = i;
        end
    end

    if total_time > 0
        lambda = lambda_sum / total_time;
    else
        lambda = NaN;
    end
%}