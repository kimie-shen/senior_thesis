% Plot the eigenvector
index = 4093;
N = 75;

X = zeros(6 * N^2, 1);
Y = zeros(6 * N^2, 1);
for i = 1:(6 * N^2)
    X(i) = xcoord(i, N);
    Y(i) = ycoord(i, N);
end
    

cn = 100;          % number of colors
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