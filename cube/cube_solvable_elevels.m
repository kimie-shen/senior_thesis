% Calculate energy spacings of solvable eigenstates
clear; 

N = 997;
energies_a1g = zeros(ceil(N^2 / 2), 3);
energies_a2g = zeros(ceil(N^2 / 2), 3);
energies_a1u = zeros(ceil(N^2 / 2), 3);
energies_a2u = zeros(ceil(N^2 / 2), 3);

index_a1g = 1;
index_a1u = 1;
index_a2g = 1;
index_a2u = 1;

% Fill in A1g and A1u wavefunctions: (n1, n2) with n1 =/= n2, both non-zero
for n1 = 1:(N - 2)
    for n2 = (n1 + 2):2:N
        if (mod(n1, 2) == 0)
            energies_a1g(index_a1g, 1) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n2 / N));
            energies_a1g(index_a1g, 2) = n1;
            energies_a1g(index_a1g, 3) = n2;
            %fprintf(['n1=' num2str(n1) ' n2=' num2str(n2) '\n'])
            index_a1g = index_a1g + 1;

            energies_a1u(index_a1u, 1) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n2 / N));
            energies_a1u(index_a1u, 2) = n1;
            energies_a1u(index_a1u, 3) = n2;
            index_a1u = index_a1u + 1;
        else
            energies_a2g(index_a2g) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n2 / N));
            energies_a2g(index_a2g, 2) = n1;
            energies_a2g(index_a2g, 3) = n2;
            index_a2g = index_a2g + 1;

            energies_a2u(index_a2u) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n2 / N));
            energies_a2u(index_a2u, 2) = n1;
            energies_a2u(index_a2u, 3) = n2;
            index_a2u = index_a2u + 1;
        end
    end
end

% Fill in A1g wavefunctions: (0, n1) with n1 even
for n1 = 0:2:N
    energies_a1g(index_a1g) = 2 * N^2 * (2 - cos(pi * 0 / N) - cos(pi * n1 / N));
    energies_a1g(index_a1g, 2) = 0;
    energies_a1g(index_a1g, 3) = n1;
    index_a1g = index_a1g + 1;
end

% Fill in A1g wavefunctions: (n1, n1) with n1 even
for n1 = 2:2:N
    energies_a1g(index_a1g) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n1 / N));
    energies_a1g(index_a1g, 2) = n1;
    energies_a1g(index_a1g, 3) = n1;
    index_a1g = index_a1g + 1;
end

% Fill in A2u wavefunctions: (n1, n1) with n1 odd
for n1 = 1:2:N
    energies_a2u(index_a2u) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n1 / N));
    energies_a2u(index_a2u, 2) = n1;
    energies_a2u(index_a2u, 3) = n1;
    index_a2u = index_a2u + 1;
end

% Remove extra rows of zeros
if (index_a1g > 1)
    energies_a1g = energies_a1g(1:(index_a1g - 1), :);
end 
if (index_a1u > 1)
    energies_a1u = energies_a1u(1:(index_a1u - 1), :);
end 
if (index_a2u > 1)
    energies_a2u = energies_a2u(1:(index_a2u - 1), :);
end 
if (index_a2g > 1)
    energies_a2g = energies_a2g(1:(index_a2g - 1), :);
end 

% Sort energies
[scratch_1, id] = sort(energies_a1g(:,1));
energies_a1g = energies_a1g(id, :);

[scratch_2, id] = sort(energies_a1u(:,1));
energies_a1u = energies_a1u(id, :);

[scratch_3, id] = sort(energies_a2g(:,1));
energies_a2g = energies_a2g(id, :);

[scratch_4, id] = sort(energies_a2u(:,1));
energies_a2u = energies_a2u(id, :);

%% Find energy level spacings
% Make vector of energy level spacings
spacings_a1g = zeros(size(energies_a1g, 1) - 1, 1);
spacings_a1u = zeros(size(energies_a1u, 1) - 1, 1);
spacings_a2u = zeros(size(energies_a2u, 1) - 1, 1);
spacings_a2g = zeros(size(energies_a2g, 1) - 1, 1);

for i = 1:(size(energies_a1g, 1) - 1)
    spacings_a1g(i) = abs(energies_a1g(i, 1) - energies_a1g(i + 1, 1));
end

for i = 1:(size(energies_a1u, 1) - 1)
    spacings_a1u(i) = abs(energies_a1u(i, 1) - energies_a1u(i + 1, 1));
end

for i = 1:(size(energies_a2u, 1) - 1)
    spacings_a2u(i) = abs(energies_a2u(i, 1) - energies_a2u(i + 1, 1));
end 

for i = 1:(size(energies_a2g, 1) - 1)
    spacings_a2g(i) = abs(energies_a2g(i, 1) - energies_a2g(i + 1, 1));
end 

% Compute r values
r_a1g = LSR(spacings_a1g(:, 1));
r_a1u = LSR(spacings_a1u(:, 1));
r_a2u = LSR(spacings_a2u(:, 1));
r_a2g = LSR(spacings_a2g(:, 1));

% Plot histogram
tiledlayout(1,4,'TileSpacing', 'tight','Padding','Tight')
nexttile
histogram(spacings_a1g(2:end), 100, 'FaceColor','black', 'EdgeColor','none')
title(['A1g: r = ' num2str(r_a1g)])
xlabel('Energy spacing')
ylabel('Counts')
%text(80, 85, ['r = ' num2str(r_a1g)])

nexttile
histogram(spacings_a1u(2:end), 100, 'FaceColor','black', 'EdgeColor','none')
title(['A1u: r = ' num2str(r_a1u)])
xlabel('Energy spacing')
ylabel('Counts')
%text(80, 80, ['r = ' num2str(r_a1u)])

nexttile
histogram(spacings_a2u(2:end), 100, 'FaceColor','black', 'EdgeColor','none')
title(['A2u: r = ' num2str(r_a2u)])
xlabel('Energy spacing')
ylabel('Counts')
%text(650, 13, ['r = ' num2str(r_a2u)])

nexttile
histogram(spacings_a2g(2:end), 100, 'FaceColor','black', 'EdgeColor','none')
title(['A2g: r = ' num2str(r_a2g)])
xlabel('Energy spacing')
ylabel('Counts')
%text(650, 13, ['r = ' num2str(r_a2g)])

set(figure(1),'position',[0,100,1200,250])

% Calculate mean level spacing ratio of given set of energy spacings

function r = LSR(spacings)
    %Allocate variable for r_n values
    r_n = zeros(size(spacings, 1) - 1, 1);
    
    % Calculate r_n variables
    for i = 1:(size(spacings, 1) - 1)
        r_n(i) = min(spacings(i), spacings(i + 1)) ...
            / max(spacings(i), spacings(i + 1));
    end

    % Calculate r variable
    r = mean(r_n);
end