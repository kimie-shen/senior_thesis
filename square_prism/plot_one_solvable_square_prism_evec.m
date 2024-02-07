% Plot a solvable cube eigenstate
n1 = 1;
n2 = 3;
N = 40;

l1 = 2; 
l2 = 1;

n 
l3 = l1;

total_sites = 2 * ((l1 * l2) + (l2 * l3) + (l3 * l1)) * N^2;

% Find x-y coordinates
X = zeros(total_sites, 1);
Y = zeros(total_sites, 1);
for i = 1:total_sites
    X(i) = xcoord(i, l1, l2, l3, N);
    Y(i) = ycoord(i, l1, l2, l3, N);
end

%f1 = @(x,y) cos(pi*(x*n1+y*n2)) + cos(pi*(y*n1-x*n2));
%f1 = @(x,y) cos(pi*(x*n1+y*n2)) + cos(pi*(y*n1-x*n2)) + cos(pi*(x*n2+y*n1)) + cos(pi*(y*n2-x*n1));
f1 = @(x,y) cos(pi*(x*n1+y*n2)) + cos(pi*(y*n1-x*n2)) - (cos(pi*(x*n2+y*n1)) + cos(pi*(y*n2-x*n1)));
eigenvector = zeros(total_sites, 1);
for i = 1:total_sites
    eigenvector(i) = f1(X(i), Y(i));
end

    
cn = 100;          % number of colors
cm = colormap(parula(cn)); 
min_of_evec = min(eigenvector);
max_of_evec = max(eigenvector);
color_indices = ceil((max_of_evec - eigenvector) / (max_of_evec - min_of_evec) * cn);
    
% Remove zeros from color indices
for j = 1:total_sites
    if (color_indices(j) == 0)
        color_indices(j) = 1;
    end
end
    
figure(1)
scatter3(X, Y, eigenvector, [], cm(color_indices, :), '.')
xlabel('x');
ylabel('y');
view(0,90);
daspect([1 1 0.2])
xlim([-l2, l1 + l2])
ylim([0, 2 * (l3 + l2)])

% Coordinate functions
function x = xcoord(index, l1, l2, l3, N) 
    a = 1 / N;
    if (index <= l1 * (l3 + l2) * N^2) % index in lower arm
        x = a/2 + mod(index - 1, l1 * N) * a;
    elseif (index <= l1 * (l3 + l2) * N^2 + l3 * (2 * l2 + l1) * N^2) % index in middle arms
        x = -l2 * N * a + a/2 + mod(index - 1 - l1 * (l3 + l2) * N^2, (2 * l2 + l1) * N) * a;
    else % index in upper arm
        x = a/2 + mod(index - (l1 * (l3 + l2) * N^2 + l3 * (2 * l2 + l1)), l1 * N) * a;
    end
end

function y = ycoord(index, l1, l2, l3, N) 
    a = 1 / N;
    if (index <= l1 * (l3 + l2) * N^2) % index in lower arm
        y = a/2 + floor((index - 1)/(l1 * N)) * a;
    elseif (index <= l1 * (l3 + l2) * N^2 + l3 * (2 * l2 + l1) * N^2) % index in middle arms
        y = (l2 + l3) * N * a + a / 2 + floor((index - 1 - l1 * (l3 + l2) * N^2)/((2 * l2 + l1) * N)) * a;
    else % index in upper arm
        y = (l2 + 2 * l3) * N * a + a / 2 + floor((index - 1 - (l1 * (l3 + l2) * N^2 + l3 * (2 * l2 + l1) * N^2))/ (l1 * N)) * a;
    end
end