%% X coordinates
% 1
new_x = 2 * (y_coord - 1) + 1 + mod(x_coord - 1, 2);

% 2
new_x = (N + 1) - (2 * (y_coord - 1) + 1) - x_coord + 1;

% 3
new_x = x_coord;

% 4
new_x = x_coord;

% 5
new_x = 2 * (y_coord - 1) + 1 - x_coord + 1;

% 6
new_x = (N + 1) - (2 * (y_coord - 1) + 1) + mod(x_coord - 1, 2);


%% Y coordinates
% 1 
new_y = (N + 1) / 2 - floor((x_coord - 1) / 2); 

% 2
new_y = (N + 1) / 2 - y_coord + 1;

% 3
new_y = (N + 1) / 2 - y_coord + 1 - floor(x_coord / 2);

% 4
new_y = 1 + floor((x_coord - 1) / 2);

% 5
new_y = y_coord + floor(x_coord / 2);

% 6
new_y = y_coord - floor(x_coord / 2);

% 7
new_y = (N + 1) / 2 - floor((x_coord - 1) / 2);

% 8
new_y = (N + 1) / 2 - y_coord + 1 + floor(x_coord / 2);

% 9
new_y = 1 + floor((x_coord - 1) / 2);

% 10
new_y = (N + 1) / 2 + 1 - y_coord;

% 11
new_y = y_coord;