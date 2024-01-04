% Calculate energy spacings of solvable eigenstates

num_sites = 75;
energies_a1u = zeros(ceil(num_sites^2 / 2), 1);

index_a1u = 1;

for n1 = 1:num_sites
    for n2 = (n1 + 2):2:(num_sites - 1)
        energies_a1u(index_a1u) = 2 * num_sites^2 * (2 - cos(pi * n1 / num_sites) - cos(pi * n2 / num_sites));
        index_a1u = index_a1u + 1;
    end
end

% Remove extra rows of zeros
if (index_a1u > 1)
    energies_a1u = energies_a1u(1:(index_a1u - 1), :);
end 

% Sort energies
energies_a1u = sort(energies_a1u);

%% Find energy level spacings
% Make vector of energy level spacings
spacings_a1u = zeros(size(energies_a1u, 1) - 1, 1);

for i = 1:(size(energies_a1u, 1) - 1)
    spacings_a1u(i) = abs(energies_a1u(i, 1) - energies_a1u(i+1, 1));
end

r = LSR(spacings_a1u);

% Plot histogram
histogram(spacings_a1u(2:end), 100, 'FaceColor','black', 'EdgeColor','none')
title(['A1u: r = ' num2str(r)])
xlabel('Energy spacing')
ylabel('Counts')


function r = LSR(spacings)
    r_n = zeros(size(spacings, 1) - 1, 1);
    
    for i = 1:(size(spacings, 1) - 1)
        r_n(i) = min(spacings(i), spacings(i + 1)) ...
            / max(spacings(i), spacings(i + 1));
    end

    r = mean(r_n);
end