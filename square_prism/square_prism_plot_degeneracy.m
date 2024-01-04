n1 = 2;
n2 = 4;

f1 = @(x, y) cos(pi*(x*n1+y*n2)).*cos(pi*(y*n1-x*n2))...
    - cos(pi*(x*-n1+y*n2)).*cos(pi*(y*-n1-x*n2));

cn = 100;  % number of colors
cm = colormap(parula(cn));

%% Plot tall arm
[x, y] = meshgrid(-2:0.005:3, -1:0.005:0);
color_indices = ceil((2 - f1(x, y)) / 4 * cn);

% Ensure color indices are within the valid range
color_indices(color_indices < 1) = 1;
color_indices(color_indices > cn) = cn;

figure(1)
scatter3(x(:), y(:), reshape(f1(x, y), 1, []), [], cm(color_indices, :), '.')
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
zlabel('$f1(x,y)$', 'interpreter', 'latex');
hold on
view(0,90);
daspect([1 1 0.2])
xlim([-2,3])

%% Plot sides
[x, y] = meshgrid(0:0.005:1, 0:0.005:2);
color_indices = ceil((2 - f1(x, y)) / 4 * cn);

% Ensure color indices are within the valid range
color_indices(color_indices < 1) = 1;
color_indices(color_indices > cn) = cn;

scatter3(x(:), y(:), reshape(f1(x, y), 1, []), [], cm(color_indices, :), '.')
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
zlabel('$f1(x,y)$', 'interpreter', 'latex');
hold on

[x, y] = meshgrid(0:0.005:1, -4:0.005:-1);
color_indices = ceil((2 - f1(x, y)) / 4 * cn);

% Ensure color indices are within the valid range
color_indices(color_indices < 1) = 1;
color_indices(color_indices > cn) = cn;

scatter3(x(:), y(:), reshape(f1(x, y), 1, []), [], cm(color_indices, :), '.')
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
zlabel('$f1(x,y)$', 'interpreter', 'latex');
hold on

view(0, 90)

xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
hold off
