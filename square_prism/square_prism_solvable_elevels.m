% Calculate energy spacings of solvable eigenstates
N = 497;

% Side length ratios
l1 = 1;
l2 = 2;

% Is l1 even or not?
l1_even = false;
if (mod(l1, 2) == 0)
    l1_even = true;
end

% Is l2 even or not?
l2_even = false;
if (mod(l2, 2) == 0)
    l2_even = true;
end

total_sites = (2 * l1^2 + 4 * l1 * l2) * N^2;


% Allocate variable space
elevels_a1g_solv_all = zeros(total_sites, 3);
elevels_a1u_solv_all = zeros(total_sites, 3);

index_a1g = 1;
index_a1u = 1;

if (l1_even == false && l2_even == true)
    elevels_b2g_solv_all = zeros(total_sites, 3);
    elevels_b2u_solv_all = zeros(total_sites, 3);

    index_b2g = 1;
    index_b2u = 1;
elseif (l1_even == true && l2_even == false)
    elevels_a2g_solv_all = zeros(total_sites, 3);
    elevels_a2u_solv_all = zeros(total_sites, 3);

    index_a2g = 1;
    index_a2u = 1;
elseif (l1_even == false && l2_even == false)
    elevels_b1g_solv_all = zeros(total_sites, 3);
    elevels_b1u_solv_all = zeros(total_sites, 3);

    index_b1g = 1;
    index_b1u = 1;
end


% Fill in A1g and A1u wavefunctions: (n1, n2) with n1 =/= n2, both non-zero
for n1 = 1:(N - 2)
    for n2 = (n1 + 2):2:N
        if (mod(n1, 2) == 0)
            elevels_a1g_solv_all(index_a1g, 1) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n2 / N));
            elevels_a1g_solv_all(index_a1g, 2) = n1;
            elevels_a1g_solv_all(index_a1g, 3) = n2;
            %fprintf(['n1=' num2str(n1) ' n2=' num2str(n2) '\n'])
            index_a1g = index_a1g + 1;

            elevels_a1u_solv_all(index_a1u, 1) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n2 / N));
            elevels_a1u_solv_all(index_a1u, 2) = n1;
            elevels_a1u_solv_all(index_a1u, 3) = n2;
            index_a1u = index_a1u + 1;
        elseif (l1_even == false && l2_even == true)
            elevels_b2g_solv_all(index_b2g) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n2 / N));
            elevels_b2g_solv_all(index_b2g, 2) = n1;
            elevels_b2g_solv_all(index_b2g, 3) = n2;
            index_b2g = index_b2g + 1;

            elevels_b2u_solv_all(index_b2u) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n2 / N));
            elevels_b2u_solv_all(index_b2u, 2) = n1;
            elevels_b2u_solv_all(index_b2u, 3) = n2;
            index_b2u = index_b2u + 1;
            
        elseif (l1_even == true && l2_even == false)
            elevels_a2g_solv_all(index_a2g) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n2 / N));
            elevels_a2g_solv_all(index_a2g, 2) = n1;
            elevels_a2g_solv_all(index_a2g, 3) = n2;
            index_a2g = index_a2g + 1;

            elevels_a2u_solv_all(index_a2u) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n2 / N));
            elevels_a2u_solv_all(index_a2u, 2) = n1;
            elevels_a2u_solv_all(index_a2u, 3) = n2;
            index_a2u = index_a2u + 1;
        elseif (l1_even == false && l2_even == false)
            elevels_b1g_solv_all(index_b1g) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n2 / N));
            elevels_b1g_solv_all(index_b1g, 2) = n1;
            elevels_b1g_solv_all(index_b1g, 3) = n2;
            index_b1g = index_b1g + 1;

            elevels_b1u_solv_all(index_b1u) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n2 / N));
            elevels_b1u_solv_all(index_b1u, 2) = n1;
            elevels_b1u_solv_all(index_b1u, 3) = n2;
            index_b1u = index_b1u + 1;
        end
    end
end

% Fill in A1g wavefunctions: (0, n1) with n1 even
for n1 = 0:2:N
    elevels_a1g_solv_all(index_a1g, 1) = 2 * N^2 * (2 - cos(pi * 0 / N) - cos(pi * n1 / N));
    elevels_a1g_solv_all(index_a1g, 2) = 0;
    elevels_a1g_solv_all(index_a1g, 3) = n1;
    index_a1g = index_a1g + 1;
end

% Fill in A1g wavefunctions: (n1, n1) with n1 even
for n1 = 2:2:N
    elevels_a1g_solv_all(index_a1g, 1) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n1 / N));
    elevels_a1g_solv_all(index_a1g, 2) = n1;
    elevels_a1g_solv_all(index_a1g, 3) = n1;
    index_a1g = index_a1g + 1;
end

% Fill in A2u wavefunctions: (n1, n1) with n1 odd
for n1 = 1:2:N
    if (l1_even == false && l2_even == true)
        elevels_b2g_solv_all(index_b2g, 1) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n1 / N));
        elevels_b2g_solv_all(index_b2g, 2) = n1;
        elevels_b2g_solv_all(index_b2g, 3) = n1;
        index_b2g = index_b2g + 1;
    elseif (l1_even == true && l2_even == false)
        elevels_a2u_solv_all(index_a2u, 1) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n1 / N));
        elevels_a2u_solv_all(index_a2u, 2) = n1;
        elevels_a2u_solv_all(index_a2u, 3) = n1;
        index_a2u = index_a2u + 1;
    elseif (l1_even == false && l2_even == false)
        elevels_b1u_solv_all(index_b1u, 1) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n1 / N));
        elevels_b1u_solv_all(index_b1u, 2) = n1;
        elevels_b1u_solv_all(index_b1u, 3) = n1;
        index_b1u = index_b1u + 1;
    end
end

% Remove extra rows of zeros
if (index_a1g > 1)
    elevels_a1g_solv_all = elevels_a1g_solv_all(1:(index_a1g - 1), :);
end 
if (index_a1u > 1)
    elevels_a1u_solv_all = elevels_a1u_solv_all(1:(index_a1u - 1), :);
end 

if (l1_even == false && l2_even == true)
    if (index_b2u > 1)
        elevels_b2u_solv_all = elevels_b2u_solv_all(1:(index_b2u - 1), :);
    end 
    if (index_b2g > 1)
        elevels_b2g_solv_all = elevels_b2g_solv_all(1:(index_b2g - 1), :);
    end
elseif (l1_even == true && l2_even == false)
    if (index_a2u > 1)
        elevels_a2u_solv_all = elevels_a2u_solv_all(1:(index_a2u - 1), :);
    end 
    if (index_a2g > 1)
        elevels_a2g_solv_all = elevels_a2g_solv_all(1:(index_a2g - 1), :);
    end 
elseif (l1_even == false && l2_even == false)
    if (index_b1u > 1)
        elevels_b1u_solv_all = elevels_b1u_solv_all(1:(index_b1u - 1), :);
    end 
    if (index_b1g > 1)
        elevels_b1g_solv_all = elevels_b1g_solv_all(1:(index_b1g - 1), :);
    end 
end

% Sort energies
elevels_a1g_solv_sorted_all = sortrows(elevels_a1g_solv_all);
elevels_a1u_solv_sorted_all = sortrows(elevels_a1u_solv_all);

if (l1_even == false && l2_even == true)
    elevels_b2g_solv_sorted_all = sortrows(elevels_b2g_solv_all);
    elevels_b2u_solv_sorted_all = sortrows(elevels_b2u_solv_all);
elseif (l1_even == true && l2_even == false)
    elevels_a2g_solv_sorted_all = sortrows(elevels_a2g_solv_all);
    elevels_a2u_solv_sorted_all = sortrows(elevels_a2u_solv_all);
elseif (l1_even == false && l2_even == false)
    elevels_b1g_solv_sorted_all = sortrows(elevels_b1g_solv_all);
    elevels_b1u_solv_sorted_all = sortrows(elevels_b1u_solv_all);
end


%% Find energy level spacings
% Make vector of energy level spacings
spacings_a1g = zeros(size(elevels_a1g_solv_all, 1) - 1, 1);
spacings_a1u = zeros(size(elevels_a1u_solv_all, 1) - 1, 1);

if (l1_even == false && l2_even == true)
    spacings_b2u = zeros(size(elevels_b2u_solv_all, 1) - 1, 1);
    spacings_b2g = zeros(size(elevels_b2g_solv_all, 1) - 1, 1);
elseif (l1_even == true && l2_even == false)
    spacings_a2u = zeros(size(elevels_a2u_solv_all, 1) - 1, 1);
    spacings_a2g = zeros(size(elevels_a2g_solv_all, 1) - 1, 1);
elseif (l1_even == false && l2_even == false)
    spacings_b1u = zeros(size(elevels_b1u_solv_all, 1) - 1, 1);
    spacings_b1g = zeros(size(elevels_b1g_solv_all, 1) - 1, 1);
end


for i = 1:(size(elevels_a1g_solv_all, 1) - 1)
    spacings_a1g(i) = abs(elevels_a1g_solv_sorted_all(i, 1) - elevels_a1g_solv_sorted_all(i + 1, 1));
end

for i = 1:(size(elevels_a1u_solv_all, 1) - 1)
    spacings_a1u(i) = abs(elevels_a1u_solv_sorted_all(i, 1) - elevels_a1u_solv_sorted_all(i + 1, 1));
end


if (l1_even == false && l2_even == true)
    for i = 1:(size(elevels_b2u_solv_all, 1) - 1)
        spacings_b2u(i) = abs(elevels_b2u_solv_sorted_all(i, 1) - elevels_b2u_solv_sorted_all(i + 1, 1));
    end 

    for i = 1:(size(elevels_b2g_solv_all, 1) - 1)
        spacings_b2g(i) = abs(elevels_b2g_solv_sorted_all(i, 1) - elevels_b2g_solv_sorted_all(i + 1, 1));
    end
elseif (l1_even == true && l2_even == false)
    for i = 1:(size(elevels_a2u_solv_all, 1) - 1)
        spacings_a2u(i) = abs(elevels_a2u_solv_sorted_all(i, 1) - elevels_a2u_solv_sorted_all(i + 1, 1));
    end 

    for i = 1:(size(elevels_a2g_solv_all, 1) - 1)
        spacings_a2g(i) = abs(elevels_a2g_solv_sorted_all(i, 1) - elevels_a2g_solv_sorted_all(i + 1, 1));
    end
elseif (l1_even == false && l2_even == false)
    for i = 1:(size(elevels_b1u_solv_all, 1) - 1)
        spacings_b1u(i) = abs(elevels_b1u_solv_sorted_all(i, 1) - elevels_b1u_solv_sorted_all(i + 1, 1));
    end 

    for i = 1:(size(elevels_b1g_solv_all, 1) - 1)
        spacings_b1g(i) = abs(elevels_b1g_solv_sorted_all(i, 1) - elevels_b1g_solv_sorted_all(i + 1, 1));
    end
end

% Compute r values
r_a1g = LSR(spacings_a1g(:, 1));
r_a1u = LSR(spacings_a1u(:, 1));

if (l1_even == false && l2_even == true)
    r_b2u = LSR(spacings_b2u(:, 1));
    r_b2g = LSR(spacings_b2g(:, 1));
elseif (l1_even == true && l2_even == false)
    r_a2u = LSR(spacings_a2u(:, 1));
    r_a2g = LSR(spacings_a2g(:, 1));
elseif (l1_even == false && l2_even == false)
    r_b1u = LSR(spacings_b1u(:, 1));
    r_b1g = LSR(spacings_b1g(:, 1));
end


% Plot histogram
bin_factor = 5;

tiledlayout(1, 4, 'TileSpacing', 'tight','Padding','Tight')
nexttile
histogram(spacings_a1g(2:end), ceil(size(elevels_a1g_solv_sorted_all, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
title('A1g Solvable')
subtitle(['r = ' num2str(r_a1g) '; total = ' num2str(size(elevels_a1g_solv_all, 1))])
xlabel('Energy spacing')
ylabel('Counts')
%text(80, 85, ['r = ' num2str(r_a1g)])

nexttile
histogram(spacings_a1u(2:end), ceil(size(elevels_a1u_solv_sorted_all, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
title('A1u Solvable')
subtitle(['r = ' num2str(r_a1u) '; total = ' num2str(size(elevels_a1u_solv_all, 1))])
xlabel('Energy spacing')
ylabel('Counts')
%text(80, 80, ['r = ' num2str(r_a1u)])

if (l1_even == false && l2_even == true)
    nexttile
    histogram(spacings_b2g(2:end), ceil(size(elevels_b2g_solv_sorted_all, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('B2g Solvable')
    subtitle(['r = ' num2str(r_b2g) '; total = ' num2str(size(elevels_b2g_solv_all, 1))])
    xlabel('Energy spacing')
    ylabel('Counts')

    nexttile
    histogram(spacings_b2u(2:end), ceil(size(elevels_b2u_solv_sorted_all, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('B2u Solvable')
    subtitle(['r = ' num2str(r_b2u) '; total = ' num2str(size(elevels_b2u_solv_all, 1))])
    xlabel('Energy spacing')
    ylabel('Counts')
elseif (l1_even == true && l2_even == false)
    nexttile
    histogram(spacings_a2g(2:end), ceil(size(elevels_a2g_solv_sorted_all, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('A2g Solvable')
    subtitle(['r = ' num2str(r_a2g) '; total = ' num2str(size(elevels_a2g_solv_all, 1))])
    xlabel('Energy spacing')
    ylabel('Counts')

    nexttile
    histogram(spacings_a2u(2:end), ceil(size(elevels_a2u_solv_sorted_all, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('A2u Solvable')
    subtitle(['r = ' num2str(r_a2u) '; total = ' num2str(size(elevels_a2u_solv_all, 1))])
    xlabel('Energy spacing')
    ylabel('Counts')
elseif (l1_even == false && l2_even == false)
    nexttile
    histogram(spacings_b1g(2:end), ceil(size(elevels_b1g_solv_sorted_all, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('B1g Solvable')
    subtitle(['r = ' num2str(r_b1g) '; total = ' num2str(size(elevels_b1g_solv_all, 1))])
    xlabel('Energy spacing')
    ylabel('Counts')

    nexttile
    histogram(spacings_b1u(2:end), ceil(size(elevels_b1u_solv_sorted_all, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('B1u Solvable')
    subtitle(['r = ' num2str(r_b1u) '; total = ' num2str(size(elevels_b1u_solv_all, 1))])
    xlabel('Energy spacing')
    ylabel('Counts')
end

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