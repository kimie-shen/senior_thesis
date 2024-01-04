
% Takes cube wavefunction and discretizes it into a vector

num_sites = 75;
n1 = -1;
n2 = 1;

wavefunction = @(x,y) cos(pi*(x*n1+y*n2))*cos(pi*(y*n1-x*n2));

wavefn_discrete = zeros(6 * num_sites^2, 1);
X = zeros(6 * num_sites^2, 1);
Y = zeros(6 * num_sites^2, 1);
for i = 1:(6 * num_sites^2)
    x = xcoord(i, num_sites);
    y = ycoord(i, num_sites);
    wavefn_discrete(i) = wavefunction(x, y);
    X(i) = x;
    Y(i) = y;
end

% Plot the discretized wavefunction
cn = 100;          % number of colors
cm = colormap(winter(cn)); 
color_indices_test = cm(ceil((wavefn_discrete + 1)*cn));
scatter3(X, Y, wavefn_discrete, [], cm(ceil((wavefn_discrete)*cn)), 'filled')
xlabel('x');
ylabel('y');


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