% Plot tetrahedron eigenvectors
N = 109;
index = 50;

X = zeros((N + 1)^2, 1);
Y = zeros((N + 1)^2, 1);
Z = eigenvectors();

a = 1 / N;

% Compute (x,y) coordinates
for i = 1:(N + 1)^2
    row_num = ceil((N + 1) - sqrt((N + 1)^2 - i));
    rows_above = (N + 1) - row_num;
    %fprintf([num2str(i) ' ' num2str(rows_above) '\n'])
    fprintf([num2str(i) ' ' num2str(i - ((N + 1)^2 - (rows_above + 1)^2)) '\n'])
    Y(i, 1) = (row_num - 1) * a * sqrt(3);
    X(i, 1) = -1 + (row_num - 1) * a + (i - ((N + 1)^2 - (rows_above + 1)^2) - 1) * a;
end

% Number of colors
cn = 100;          
cm = colormap(parula(cn)); 
min_of_evec = min(eigenvectors(:, index));
max_of_evec = max(eigenvectors(:, index));
color_indices = ceil((max_of_evec - eigenvectors(:, index)) / (max_of_evec - min_of_evec) * cn);
    
% Remove zeros from color indices
for j = 1:((N + 1)^2)
    if (color_indices(j) == 0)
        color_indices(j) = 1;
    end
end
    
scatter3(X, Y, eigenvectors(:, index), [], cm(color_indices, :), '.')
xlabel('x');
ylabel('y');
view(0, 90);
title(['wavefunction index = ' num2str(index)])
daspect([1 1 0.1])

%scatter(X, Y)
%view(0,90);
%daspect([1 1 0.2])
%xlim([-1 - a, 1 + a])
%ylim([-a, sqrt(3) + a])