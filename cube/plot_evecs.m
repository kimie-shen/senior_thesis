% Plot the eigenvector

for k = 409:409
    num_sites = 75;
    start_index = k * 10 + 1;
    end_index = start_index + 9;
    
    
    X = zeros(6 * num_sites^2, 1);
    Y = zeros(6 * num_sites^2, 1);
    for i = 1:(6 * num_sites^2)
        X(i) = xcoord(i, num_sites);
        Y(i) = ycoord(i, num_sites);
    end
    
    for i=start_index:end_index
        cn = 100;          % number of colors
        cm = colormap(parula(cn)); 
        min_of_evec = min(eigenvectors(:, i));
        max_of_evec = max(eigenvectors(:, i));
        color_indices = ceil((max_of_evec - eigenvectors(:, i)) / (max_of_evec - min_of_evec) * cn);
    
        % Remove zeros from color indices
        for j = 1:(6 * num_sites^2)
            if (color_indices(j) == 0)
                color_indices(j) = 1;
            end
        end
    
        if (i == start_index)
            figure(1)
            tiledlayout(2, ceil((end_index - start_index + 1)/2))
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
    set(figure(1),'position',[0,0,1600,900])
    saveas(gcf, ['cube_sigma2_evecs/' num2str(k * 10 + 1) '-' num2str(k * 10 + 10) '.jpeg']);
end

% Coordinate functions
function x = xcoord(index, num_sites) 
    a = 1 / num_sites;
    if (index <= 2 * num_sites^2) % index in lower arm
        x = a/2 + mod(index - 1, num_sites) * a;
    elseif (index <= 5 * num_sites^2) % index in middle arms
        x = - num_sites * a + a/2 + mod(index - 1 - 2 * num_sites^2, 3 * num_sites) * a;
    else % index in upper arm
        x = a/2 + mod(index - 1, num_sites) * a;
    end
end

function y = ycoord(index, num_sites) 
    a = 1 / num_sites;
    if (index <= 2 * num_sites^2) % index in lower arm
        y = a/2 + floor((index - 1)/num_sites) * a;
    elseif (index <= 5 * num_sites^2) % index in middle arms
        y = 2 * num_sites * a + a / 2 + floor((index - 1 - 2 * num_sites^2)/(3 * num_sites)) * a;
    else % index in upper arm
        y = 3 * num_sites * a + a / 2 + floor((index - 1 - 5 * num_sites^2)/num_sites) * a;
    end
end