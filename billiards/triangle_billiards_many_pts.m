%% Simulation of billiards on triangle
% Initial conditions
num_points = 1000;
t_total = 1.5;
theta = pi/2 - 0.01;
dtheta = 0;

x_i = zeros(num_points, 1);
y_i = zeros(num_points, 1);

square_width = 0.1;
square_center_x = 0.5;
square_center_y = 0.05;

% Randomly generate points in square 
for i = 1:num_points
    x_i(i) = square_center_x - (square_width / 2) + square_width * rand();
    y_i(i) = square_center_y - (square_width / 2) + square_width * rand();
end

scatter(x_i, y_i)

% Randomly generate theta
theta_i = zeros(num_points, 1);

for i = 1:num_points
    theta_i(i) = theta - dtheta + 2 * dtheta * rand();
end


%% Find endpoint after billiards
x_f = zeros(num_points, 1);
y_f = zeros(num_points, 1);

for i = 1:num_points
    [x_f(i), y_f(i)] = billiards(x_i(i), y_i(i), theta_i(i), t_total);
end

%% Plot
figure(1)
scatter(x_f, y_f, '.')
xlim([0, 1])
ylim([0, sqrt(3)/2])
daspect([1 1 1])
hold on
    
% Plot boundary
plot([0, 1], [0, 0], 'LineWidth', 1, 'Color', 'black')
plot([1, 1/2], [0, sqrt(3) / 2], 'LineWidth', 1, 'Color', 'black')
plot([0, 1/2], [0, sqrt(3) / 2], 'LineWidth', 1, 'Color', 'black')

%% Functions
function [x_end, y_end] = billiards(x_0, y_0, theta, t_total)
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
        [x_right, y_right, theta_right, t_right] = right(x(hits), y(hits), theta);
        [x_left, y_left, theta_left, t_left] = left(x(hits), y(hits), theta);
    
        % Determine which collision occurs
        if (t_bottom > 1e-8)
            t_next = t_bottom;
            x(hits + 1) = x_bottom; % Record x corrdinate of collision site
            y(hits + 1) = y_bottom; % Record y corrdinate of collision site
            theta = theta_bottom; % update angle
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
    x_end = x(hits) + t_remaining * cos(theta_penultimate);
    y_end = y(hits) + t_remaining * sin(theta_penultimate);
    
end

% Reflection off bottom line
function [x_f, y_f, theta_f, t_hit] = bottom(x_i, y_i, theta_i)
    t_hit = -y_i / sin(theta_i);
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
    t_hit = (-sqrt(3) * (x_i - 1) - y_i) / (sin(theta_i) + sqrt(3) * cos(theta_i));
    x_f = x_i + t_hit * cos(theta_i);
    y_f = y_i + t_hit * sin(theta_i);
    theta_f = (4 * pi / 3) - theta_i;
end

