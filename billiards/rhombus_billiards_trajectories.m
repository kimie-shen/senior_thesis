%% Simulation of billiards on rhombus
% Initial conditions
x_0 = [1 - 0.1, 1 - 0.06, 1 - 0.02, 1 + 0.01, 1 + 0.03, 1 + 0.05];
y_0 = [0, 0, 0, 0.01 * sqrt(3), 0.03 * sqrt(3), 0.05* sqrt(3)];
theta = (2 * pi / 3);
t_total = 0.75;
color_list = [0.635, 0.078, 0.184; 0.850, 0.325, 0.098; 0.929, 0.694, 0.125; 0.466, 0.674, 0.188; 0.000, 0.447, 0.741; 0.494, 0.184, 0.556]; % red, orange, yellow, green, blue, purple
num_lines = 6;

% Length
for i = 1:num_lines
    billiards(x_0(i), y_0(i), theta, t_total, color_list(i, :))
    hold on
end
hold off


function billiards(x_0, y_0, theta, t_total, color_rgb)
    % Set variables
    x = [x_0];
    y = [y_0];
    
    hits = 0;
    t_traveled = 0;
    t_penultimate = 0;
    theta_penultimate = 0;
    
    while (t_traveled <= t_total)
        t_penultimate = t_traveled;
        theta_penultimate = theta;
        t_next = Inf;
        hits = hits + 1;
        
        % Compute time of collision at each wall
        [x_bottom, y_bottom, theta_bottom, t_bottom] = bottom(x(hits), y(hits), theta);
        [x_top, y_top, theta_top, t_top] = top(x(hits), y(hits), theta);
        [x_right, y_right, theta_right, t_right] = right(x(hits), y(hits), theta);
        [x_left, y_left, theta_left, t_left] = left(x(hits), y(hits), theta);
    
        % Determine which collision occurs
        if (t_bottom > 1e-8)
            t_next = t_bottom;
            x(hits + 1) = x_bottom; % Record x corrdinate of collision site
            y(hits + 1) = y_bottom; % Record y corrdinate of collision site
            theta = theta_bottom; % update angle
        end
        if (t_top > 1e-8 && t_top < t_next)
            t_next = t_top;
            x(hits + 1) = x_top;
            y(hits + 1) = y_top;
            theta = theta_top;
        end
        if (t_right > 1e-8 && t_right < t_next)
            t_next = t_right;
            x(hits + 1) = x_right;
            y(hits + 1) = y_right;
            theta = theta_right;
        end
        if (t_left > 1e-8 && t_left < t_next)
            t_next = t_left;
            x(hits + 1) = x_left;
            y(hits + 1) = y_left;
            theta = theta_left;
        end
    
        t_traveled = t_traveled + t_next; % update time traveled
    end
    %% Compute end point
    
    % Remove last collision
    x = x(1:(end-1));
    y = y(1:(end-1));
    
    % Compute time left
    t_remaining = t_total - t_penultimate;
    
    % Compute end location
    x(hits + 1) = x(hits) + t_remaining * cos(theta_penultimate);
    y(hits + 1) = y(hits) + t_remaining * sin(theta_penultimate);
    
    %% Plot
    figure(1)
    plot(x, y, 'Color', color_rgb)
    xlim([0, 3/2])
    ylim([0, sqrt(3)/2])
    daspect([1 1 1])
    hold on
    
    % Plot boundary
    plot([0,1], [0,0], 'LineWidth', 1, 'Color', 'black')
    plot([1, 3/2], [0, sqrt(3) / 2], 'LineWidth', 1, 'Color', 'black')
    plot([0, 1/2], [0, sqrt(3) / 2], 'LineWidth', 1, 'Color', 'black')
    plot([1/2, 3/2], [sqrt(3) / 2, sqrt(3) / 2], 'LineWidth', 1, 'Color', 'black')
end

% Reflection off bottom line
function [x_f, y_f, theta_f, t_hit] = bottom(x_i, y_i, theta_i)
    t_hit = -y_i / sin(theta_i);
    x_f = x_i + t_hit * cos(theta_i);
    y_f = y_i + t_hit * sin(theta_i);
    theta_f = 2 * pi - theta_i;
end

% Reflection off top line
function [x_f, y_f, theta_f, t_hit] = top(x_i, y_i, theta_i)
    t_hit = ((sqrt(3) / 2) - y_i) / sin(theta_i);
    x_f = x_i + t_hit * cos(theta_i);
    y_f = y_i + t_hit * sin(theta_i);
    theta_f = 2 * pi - theta_i;
end

% Reflection off left side
function [x_f, y_f, theta_f, t_hit] = left(x_i, y_i, theta_i)
    t_hit = (sqrt(3) * x_i - y_i) / (sin(theta_i) - sqrt(3) * cos(theta_i));
    x_f = x_i + t_hit * cos(theta_i);
    y_f = y_i + t_hit * sin(theta_i);
    theta_f = (2 * pi / 3) - theta_i;
    %fprintf([num2str(t_hit) ' '  num2str(x_i) ' ' num2str(y_i) ' ' num2str(theta_i) ' \n'])
end

% Reflection off right side
function [x_f, y_f, theta_f, t_hit] = right(x_i, y_i, theta_i)
    t_hit = (sqrt(3) * (x_i - 1) - y_i) / (sin(theta_i) - sqrt(3) * cos(theta_i));
    x_f = x_i + t_hit * cos(theta_i);
    y_f = y_i + t_hit * sin(theta_i);
    theta_f = (2 * pi / 3) - theta_i;
end

