% Plot a solvable cube eigenstate
n1 = 3;
n2 = 3;
num_sites = 75;

% Find x-y coordinates
X = zeros(6 * num_sites^2, 1);
Y = zeros(6 * num_sites^2, 1);
for i = 1:(6 * num_sites^2)
    X(i) = xcoord(i, num_sites);
    Y(i) = ycoord(i, num_sites);
end

f1 = @(x,y) cos(pi*(x*n1+y*n2))*cos(pi*(y*n1-x*n2));
eigenvector = zeros(6 * num_sites^2, 1);
for i = 1:(6 * num_sites^2)
    eigenvector(i) = f1(X(i), Y(i));
end

    
cn = 100;          % number of colors
cm = colormap(parula(cn)); 
min_of_evec = min(eigenvector);
max_of_evec = max(eigenvector);
color_indices = ceil((max_of_evec - eigenvector) / (max_of_evec - min_of_evec) * cn);
    
% Remove zeros from color indices
for j = 1:(6 * num_sites^2)
    if (color_indices(j) == 0)
        color_indices(j) = 1;
    end
end
    
scatter3(X, Y, eigenvector, [], cm(color_indices, :), '.')
xlabel('x');
ylabel('y');
view(0,90);
daspect([1 1 0.2])

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