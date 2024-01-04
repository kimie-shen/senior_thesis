clear
N = 40; % Number of sites per side
L = 1; % Spacing
H = zeros(N^2, N^2); % Hamiltonian matrix

% Set matrix values
for i = 1:N^2
    H(i, i) = -4/L^2;
    H(i, upper(i, N)) = 1/L^2;
    H(i, lower(i, N)) = 1/L^2;
    H(i, right(i, N)) = 1/L^2;
    H(i, left(i, N)) = 1/L^2;
end

% Find eigenvalues
e = eig(H);

% Plot eigenvalues
figure(3)
for i = 1:N^2
    x = [-1, 0];
    y = [-e(i), -e(i)];
    plot(x, y,'-r')
    hold on
end

for i = 0:N
    for j = 0:N
        kx = 2 * pi * i / N;
        ky = 2 * pi * j / N;
        E = -(2/L^2)*(cos(kx*L) + cos(ky*L) - 2);
        x = [0, 1];
        y = [E, E];
        plot(x, y, '-b');
        hold on
    end
end

hold off

% Functions to find indices
function y = upper(x, num_sites)
    if (x > num_sites) % x not in top row
        y = x - num_sites;
    else
        y = num_sites * (num_sites - 1) + x;
    end
end

function y = lower(x, num_sites) 
    if (x <= num_sites * (num_sites - 1)) % x not in bottom row
        y = x + num_sites;
    else
        y = mod(x - 1, num_sites) + 1;
    end
end

function y = right(x, num_sites)
    if (floor(x/num_sites) < x/num_sites) % x not in right column
        y = x + 1;
    else
        y = x - num_sites + 1;
    end
end

function y = left(x, num_sites)
    if (floor((x-1)/num_sites) < (x-1)/num_sites) % x not in left column
        y = x - 1;
    else
        y = x + num_sites - 1;
    end
end
