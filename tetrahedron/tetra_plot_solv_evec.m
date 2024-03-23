%% Plot ten tetrahedron eigenvectors
clear;
% Plot tetrahedron eigenvectors
N = 150;
n1 = 1;
n2 = 9;

n3 = 0;
n4 = 2;

index = 50;

X = zeros((N + 1)^2, 1);
Y = zeros((N + 1)^2, 1);

a = 1 / N;

% Compute (x,y) coordinates
for i = 1:(N + 1)^2
    row_num = ceil((N + 1) - sqrt((N + 1)^2 - i));
    rows_above = (N + 1) - row_num;
    %fprintf([num2str(i) ' ' num2str(rows_above) '\n'])
    %fprintf([num2str(i) ' ' num2str(i - ((N + 1)^2 - (rows_above + 1)^2)) '\n'])
    Y(i, 1) = (row_num - 1) * a * sqrt(3);
    X(i, 1) = (row_num - 1) * a + (i - ((N + 1)^2 - (rows_above + 1)^2) - 1) * a;
end

f1 = @(x,y) cos(pi*(x*n1+y*n2 / sqrt(3)));
C3_f1 = @(x,y) cos(pi*(x*((n1 - n2) / 2) + y * ((3 * n1 + n2) / 2) / sqrt(3)));
C3_C3_f1 = @(x,y) cos(pi*(x*((- n1 - n2) / 2) + y * ((3 * n1 - n2) / 2) / sqrt(3)));
f1_all = @(x,y) f1(x,y) + C3_f1(x,y) + C3_C3_f1(x,y);

% Next set (n3, n4)
f2 = @(x,y) cos(pi*(x*n3+y*n4 / sqrt(3)));
C3_f2 = @(x,y) cos(pi*(x*((n3 - n4) / 2) + y * ((3 * n3 + n4) / 2) / sqrt(3)));
C3_C3_f2 = @(x,y) cos(pi*(x*((- n3 - n4) / 2) + y * ((3 * n3 - n4) / 2) / sqrt(3)));
f2_all = @(x,y) f2(x,y) + C3_f2(x,y) + C3_C3_f2(x,y);

% Switch n1 and n2
r_f1 = @(x,y) cos(pi*(x*n2+y*-n1 / sqrt(3)));
r_C3_f1 = @(x,y) cos(pi*(x*((n2 + n1) / 2) + y * ((3 * n2 - n1) / 2) / sqrt(3)));
r_C3_C3_f1 = @(x,y) cos(pi*(x*((- n2 + n1) / 2) + y * ((3 * n2 + n1) / 2) / sqrt(3)));
r_f1_all = @(x,y) x_f1(x,y) + x_C3_f1(x,y) + x_C3_C3_f1(x,y);

% Reflected about x axis
x_f1 = @(x,y) cos(pi*(x*n1+y*-n2 / sqrt(3)));
x_C3_f1 = @(x,y) cos(pi*(x*((n1 + n2) / 2) + y * ((3 * n1 - n2) / 2) / sqrt(3)));
x_C3_C3_f1 = @(x,y) cos(pi*(x*((- n1 + n2) / 2) + y * ((3 * n1 + n2) / 2) / sqrt(3)));
x_f1_all = @(x,y) x_f1(x,y) + x_C3_f1(x,y) + x_C3_C3_f1(x,y);

% Reflected about y axis
y_f1 = @(x,y) cos(pi*(x*-n1+y*n2 / sqrt(3)));
y_C3_f1 = @(x,y) cos(pi*(x*((-n1 - n2) / 2) + y * ((3 * -n1 + n2) / 2) / sqrt(3)));
y_C3_C3_f1 = @(x,y) cos(pi*(x*((n1 - n2) / 2) + y * ((3 * -n1 - n2) / 2) / sqrt(3)));
y_f1_all = @(x,y) y_f1(x,y) + y_C3_f1(x,y) + y_C3_C3_f1(x,y);

% Shifted up by sqrt(3)L
y_f1 = @(x,y) cos(pi*(x*-n1+y*n2 / sqrt(3)));
y_C3_f1 = @(x,y) cos(pi*(x*((-n1 - n2) / 2) + y * ((3 * -n1 + n2) / 2) / sqrt(3)));
y_C3_C3_f1 = @(x,y) cos(pi*(x*((n1 - n2) / 2) + y * ((3 * -n1 - n2) / 2) / sqrt(3)));
y_f1_all = @(x,y) y_f1(x,y) + y_C3_f1(x,y) + y_C3_C3_f1(x,y);

Z = zeros((N + 1)^2, 1);
for i = 1:((N + 1)^2)
    Z(i) = f1_all(X(i), Y(i)) - y_f1_all(X(i), Y(i));
end

    
cn = 100;          % number of colors
cm = colormap(parula(cn)); 
min_of_evec = min(Z);
max_of_evec = max(Z);
color_indices = ceil((max_of_evec - Z) / (max_of_evec - min_of_evec) * cn);
    
% Remove zeros from color indices
for j = 1:((N + 1)^2)
    if (color_indices(j) == 0)
        color_indices(j) = 1;
    end
end
figure(1)
scatter3(X, Y, Z, [], cm(color_indices, :), '.')
xlabel('$x_1$','Interpreter','Latex');
ylabel('$x_2$','Interpreter','Latex');
view(0, 90);
daspect([1 1 0.2])
