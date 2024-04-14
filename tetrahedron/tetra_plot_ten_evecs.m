%% Plot ten tetrahedron eigenvectors

% Plot tetrahedron eigenvectors
%N = 31;
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
    %fprintf([num2str(i) ' ' num2str(i - ((N + 1)^2 - (rows_above + 1)^2)) '\n'])
    Y(i, 1) = (row_num - 1) * a * sqrt(3);
    X(i, 1) = -1 + (row_num - 1) * a + (i - ((N + 1)^2 - (rows_above + 1)^2) - 1) * a;
end


% Plot the eigenvector
folderName = 'tetra_sigma_d_evecs_3';
mkdir(folderName);

for k = 0:20
    start_index = k * 10 + 1;
    end_index = start_index + 9;
    
    
    for i=start_index:end_index
        cn = 100;          % number of colors
        cm = colormap(parula(cn)); 
        min_of_evec = min(eigenvectors(:, i));
        max_of_evec = max(eigenvectors(:, i));
        color_indices = ceil((max_of_evec - eigenvectors(:, i)) / (max_of_evec - min_of_evec) * cn);
    
        % Remove zeros from color indices
        for j = 1:((N + 1)^2)
            if (color_indices(j) == 0)
                color_indices(j) = 1;
            end
        end
    
        if (i == start_index)
            figure(1)
            tiledlayout(2, ceil((end_index - start_index + 1)/2), 'TileSpacing', 'tight','Padding','Tight')
            nexttile
            scatter3(X, Y, eigenvectors(:, i), [], cm(color_indices, :), '.')
            xlabel('x');
            ylabel('y');
            view(0,90);
            title(['wavefunction index = ' num2str(i)])
            daspect([1 1 0.2])
        else
            nexttile
            scatter3(X, Y, eigenvectors(:, i), [], cm(color_indices, :), '.')
            xlabel('x');
            ylabel('y');
            view(0,90);
            title(['wavefunction index = ' num2str(i)])
            daspect([1 1 0.2])
        end
    end
    set(figure(1),'position',[0,0,1600,600])
    saveas(gcf, [folderName '/tetra_' num2str(k * 10 + 1) '-' num2str(k * 10 + 10) '.jpeg']);
end