%% Calculate r values for various N
clear
close all
tolerance = 1e-8;  % Tolerance for degenerate energy levels
potential = -0;

% Side length ratios
l1 = 1;
l2 = 2;
l3 = 3;

% Find max N
upper_lim = 33000; % Max number of sites allowed in memory
max_N = floor(sqrt(upper_lim / (2 * ((l1 * l2) + (l2 * l3) + (l3 * l1)))));
fprintf(['Max N = ' num2str(max_N) '\n'])

N_nums = primes(max_N);

%% Resume Program
num_irreps = 8;
num_symmetries = 4;
solvable_irreps = 2;

r_array = zeros(size(N_nums, 1), (num_irreps + 1));
size_array = zeros(size(N_nums, 1), (num_irreps + 1));
solvable_prop = zeros(size(N_nums, 1), solvable_irreps);

folderName = ['l1=' num2str(l1) '_l2=' num2str(l2) '_l3=' num2str(l3)];
mkdir(folderName);

index_n = 1;
index_r = 1;
index_size = 1;

%% Start loop
for n = 1:size(N_nums, 2)

    N = N_nums(n);
    fprintf(['\n    N = ' num2str(N) '\n'])
    total_sites = 2 * ((l1 * l2) + (l2 * l3) + (l3 * l1)) * N^2; % Total number of sites in unfolded polyhedron

    %% Diagonalize H and sigma_d matrices
    tic;

    % Generate Hamiltonian matrix
    L = 1/N; % Spacing
    H = zeros(total_sites, total_sites); % Hamiltonian matrix
    
    for i = 1:total_sites
        H(i, i) = 4/L^2;
        %fprintf([num2str(i) ' ' num2str(upper(i, N, l1, l2, l3, total_sites)) ' ' num2str(lower(i, N, l1, l2, l3, total_sites)) ' ' ...
            %num2str(right(i, N, l1, l2, l3, total_sites)) ' ' num2str(left(i, N, l1, l2, l3, total_sites)) '\n'])
        H(i, upper(i, N, l1, l2, l3, total_sites)) = -1/L^2;
        H(i, lower(i, N, l1, l2, l3, total_sites)) = -1/L^2;
        H(i, right(i, N, l1, l2, l3, total_sites)) = -1/L^2;
        H(i, left(i, N, l1, l2, l3, total_sites)) = -1/L^2;       
    end
    
    % Generate matrix for c2x symmetry 
    C2x = zeros(total_sites, total_sites);
    for i = 1:(total_sites)
        C2x(i, c2x_index(i, N, l1, l2, l3, total_sites)) = 1;
        %fprintf([num2str(i) ' ' num2str(c2x_index(i, N, l1, l2, l3, total_sites)) '\n'])
    end
    
    %% Calculate eigenvectors and eigenvalues
    [eigenvectors, eigenvalues] = eig(H + ((1 + sqrt(5))/2) * C2x);
    fprintf('evals and evecs done \n')
    toc
    
    %% Generate Symmetry Matrices
    
    % Generate matrix for C2y
    C2y = zeros(total_sites, total_sites);
    for i = 1:(total_sites)
        C2y(i, c2y_index(i, N, l1, l2, l3, total_sites)) = 1;
        %fprintf([num2str(i) ' ' num2str(c2y_index(i, N, l1, l2, l3, total_sites)) '\n'])
    end

    % Generate matrix for C2z
    C2z = zeros(total_sites, total_sites);
    for i = 1:(total_sites)
        C2z(i, c2z_index(i, N, l1, l2, l3, total_sites)) = 1;
        %fprintf([num2str(i) ' ' num2str(c2z_index(i, N, l1, l2, l3, total_sites)) '\n'])
    end

    % Generate matrix for i
    inv = zeros(total_sites, total_sites);
    for i = 1:(total_sites)
        inv(i, inv_index(i, N, l1, l2, l3, total_sites)) = 1;
        %fprintf([num2str(i) ' ' num2str(inv_index(i, N, l1, l2, l3, total_sites)) '\n'])
    end

    
    %% Calculate characters of symmetries
    tic;
    % Preallocate data array
    c2x_evals = zeros(total_sites, 1);
    c2y_evals = zeros(total_sites, 1);
    c2z_evals = zeros(total_sites, 1);
    inv_evals = zeros(total_sites, 1);
    
    % Find energy and c2x eigenvluaes
    for i = 1:(total_sites)
        c2x_evals(i) = (eigenvectors(:, i).' * C2x * eigenvectors(:, i));
        eigenvalues(i, i) = eigenvalues(i, i) - ((1 + sqrt(5))/2) * c2x_evals(i);
    end
    fprintf('c2x done \n')
    
    % Find c2y eigenvalues
    for i = 1:(total_sites)
        c2y_evals(i) = (eigenvectors(:, i).' * C2y * eigenvectors(:, i));
    end
    fprintf('c2y done \n')

    % Find c2z eigenvalues
    for i = 1:(total_sites)
        c2z_evals(i) = (eigenvectors(:, i).' * C2z * eigenvectors(:, i));
    end
    fprintf('c2z done \n')

    % Find inv eigenvalues
    for i = 1:(total_sites)
        inv_evals(i) = (eigenvectors(:, i).' * inv * eigenvectors(:, i));
    end
    fprintf('inv done \n')
    
    
    %% Analyze Degeneracies
    % Reorder energies and evecs
    energy_levels = zeros(total_sites, num_symmetries + 3);
    [energies, id] = sort(diag(eigenvalues));
    eigenvectors = eigenvectors(:, id);
    c2x_evals_sorted = c2x_evals(id);
    c2y_evals_sorted = c2y_evals(id);
    c2z_evals_sorted = c2z_evals(id);
    inv_evals_sorted = inv_evals(id);

    % Analyze for degeneracies
    energy = energies(1);
    degeneracy = 1;
    index = 1;
    trace_c2x = c2x_evals_sorted(1);
    trace_c2y = c2y_evals_sorted(1);
    trace_c2z = c2z_evals_sorted(1);
    trace_inv = inv_evals_sorted(1);

    
    for i = 2:(total_sites)
        if (abs(energies(i) - energy) < tolerance) % Next energy is degenerate
            degeneracy = degeneracy + 1;
            trace_c2x = trace_c2x + c2x_evals_sorted(i);
            trace_c2y = trace_c2y + c2y_evals_sorted(i);
            trace_c2z = trace_c2z + c2z_evals_sorted(i);
            trace_inv = trace_inv + inv_evals_sorted(i);
        else % Next energy is new
            % Record stats for previous energy level
            energy_levels(index, 1) = energy;
            energy_levels(index, 2) = degeneracy;
            energy_levels(index, 3) = trace_c2x;
            energy_levels(index, 4) = trace_c2y;
            energy_levels(index, 5) = trace_c2z;
            energy_levels(index, 6) = trace_inv;
            energy_levels(index, 7) = i - 1;
    
            % Reset variables for new energy level
            energy = energies(i);
            degeneracy = 1;
            index = index + 1;
            trace_c2x = c2x_evals_sorted(i);
            trace_c2y = c2y_evals_sorted(i);
            trace_c2z = c2z_evals_sorted(i);
            trace_inv = inv_evals_sorted(i);            
        end
    
        % Record energy level if we reach the end
        if (i == total_sites)
            energy_levels(index, 1) = energy;
            energy_levels(index, 2) = degeneracy;
            energy_levels(index, 3) = trace_c2x;
            energy_levels(index, 4) = trace_c2y;
            energy_levels(index, 5) = trace_c2z;
            energy_levels(index, 6) = trace_inv;
            energy_levels(index, 7) = i;
            
        end
    end
    
    energy_levels = energy_levels(1:index, :);
    
    %% Separate the representation classes
    % Allocate space for each irrep
    elevels_ag = zeros(total_sites, num_symmetries + 3);
    elevels_b1g = zeros(total_sites, num_symmetries + 3);
    elevels_b2g = zeros(total_sites, num_symmetries + 3);
    elevels_b3g = zeros(total_sites, num_symmetries + 3);
    elevels_au = zeros(total_sites, num_symmetries + 3);
    elevels_b1u = zeros(total_sites, num_symmetries + 3);
    elevels_b2u = zeros(total_sites, num_symmetries + 3);
    elevels_b3u = zeros(total_sites, num_symmetries + 3);
    
    index_ag = 1;
    index_b1g = 1;
    index_b2g = 1;
    index_b3g = 1;
    index_au = 1;
    index_b1u = 1;
    index_b2u = 1;
    index_b3u = 1;

    % Counter for accidental degeneracies
    ag_au = 0;
    b2g_b2u = 0;
    b2g_b3g = 0;
    b2u_b3u = 0;
    b1g_b1u = 0;
    ag_b3g = 0;
    b1g_b2g = 0;
    au_b3u = 0;
    b1u_b2u = 0;
    b3g_b3u = 0;
    
    % Round the energy_levels characters
    energy_levels_rounded = zeros(size(energy_levels, 1), num_symmetries + 3);
    energy_levels_rounded(:, 2:(num_symmetries + 2)) = round(energy_levels(:, 2:(num_symmetries + 2)));
    energy_levels_rounded(:, 1) = energy_levels(:, 1);
    energy_levels_rounded(:, num_symmetries + 3) = energy_levels(:, num_symmetries + 3);
    
    % Fill in elevel info in each irrep
    for i = 1:size(energy_levels_rounded, 1)
        trace_e = energy_levels_rounded(i, 2);
        trace_c2x = energy_levels_rounded(i, 3);
        trace_c2y = energy_levels_rounded(i, 4);
        trace_c2z = energy_levels_rounded(i, 5);
        trace_inv = energy_levels_rounded(i, 6);
    
        traces = [trace_e, trace_c2x, trace_c2y, trace_c2z, trace_inv];
    
        if (isequal(traces, [1, 1, 1, 1, 1])) % Ag
            elevels_ag(index_ag, :) = energy_levels_rounded(i, :);
            index_ag = index_ag + 1;
        elseif (isequal(traces, [1, 1, -1, -1, 1])) % B1g
            elevels_b1g(index_b1g, :) = energy_levels_rounded(i, :);
            index_b1g = index_b1g + 1;
        elseif (isequal(traces, [1, -1, 1, -1, 1])) % B2g
            elevels_b2g(index_b2g, :) = energy_levels_rounded(i, :);
            index_b2g = index_b2g + 1;
        elseif (isequal(traces, [1, -1, -1, 1, 1])) % B3g
            elevels_b3g(index_b3g, :) = energy_levels_rounded(i, :);
            index_b3g = index_b3g + 1;
        elseif (isequal(traces, [1, 1, 1, 1, -1])) % Au
            elevels_au(index_au, :) = energy_levels_rounded(i, :);
            index_au = index_au + 1;
        elseif (isequal(traces, [1, 1, -1, -1, -1])) % B1u
            elevels_b1u(index_b1u, :) = energy_levels_rounded(i, :);
            index_b1u = index_b1u + 1;
        elseif (isequal(traces, [1, -1, 1, -1, -1])) % B2u
            elevels_b2u(index_b2u, :) = energy_levels_rounded(i, :);
            index_b2u = index_b2u + 1;
        elseif (isequal(traces, [1, -1, -1, +1, -1])) % B3u
            elevels_b3u(index_b3u, :) = energy_levels_rounded(i, :);
            index_b3u = index_b3u + 1;
        elseif (isequal(traces, [2, 2, 2, 2, 0])) % Accidental degeneracy Ag + Au
            elevels_ag(index_ag, :) = energy_levels_rounded(i, :);
            index_ag = index_ag + 1;
            elevels_au(index_au, :) = energy_levels_rounded(i, :); 
            index_au = index_au + 1;
            ag_au = ag_au + 1;
        elseif (isequal(traces, [2, 2, -2, -2, 0])) % Accidental degeneracy B1g + B1u
            elevels_b1g(index_b1g, :) = energy_levels_rounded(i, :);
            index_b1g = index_b1g + 1;
            elevels_b1u(index_b1u, :) = energy_levels_rounded(i, :); 
            index_b1u = index_b1u + 1;
            b1g_b1u = b1g_b1u + 1;
        elseif (isequal(traces, [2, -2, 2, -2, 0])) % Accidental degeneracy B2g + B2u
            elevels_b2g(index_b2g, :) = energy_levels_rounded(i, :);
            index_b2g = index_b2g + 1;
            elevels_b2u(index_b2u, :) = energy_levels_rounded(i, :); 
            index_b2u = index_b2u + 1;
            b2g_b2u = b2g_b2u + 1;
        elseif (isequal(traces, [2, -2, -2, 2, 0])) % Accidental degeneracy B3g + B3u
            elevels_b3g(index_b3g, :) = energy_levels_rounded(i, :);
            index_b3g = index_b3g + 1;
            elevels_b3u(index_b3u, :) = energy_levels_rounded(i, :); 
            index_b3u = index_b3u + 1;
            b3g_b3u = b3g_b3u + 1;
        else 
            fprintf([num2str(energy_levels_rounded(i, 1)) ' ' num2str(energy_levels_rounded(i, 2)) ' ' ...
                num2str(energy_levels_rounded(i, 3)) ' ' num2str(energy_levels_rounded(i, 4)) ' ' ...
                num2str(energy_levels_rounded(i, 5)) ' ' num2str(energy_levels_rounded(i, 6)) ' ' num2str(energy_levels_rounded(i, 7)) '\n'])
        end
    end
    
    % Remove extra rows of zeros
    if (index_ag > 1)
        elevels_ag = elevels_ag(1:(index_ag - 1), :);
    end 
    
    if (index_b1g > 1)
        elevels_b1g = elevels_b1g(1:(index_b1g - 1), :);
    end
    
    if (index_b2g > 1)
        elevels_b2g = elevels_b2g(1:(index_b2g - 1), :);
    end
    
    if (index_b3g > 1)
        elevels_b3g = elevels_b3g(1:(index_b3g - 1), :);
    end
    
    if (index_au > 1)
        elevels_au = elevels_au(1:(index_au - 1), :);
    end 
        
    if (index_b1u > 1)
        elevels_b1u = elevels_b1u(1:(index_b1u - 1), :);
    end
    
    if (index_b2u > 1)
        elevels_b2u = elevels_b2u(1:(index_b2u - 1), :);
    end
    
    if (index_b3u > 1)
        elevels_b3u = elevels_b3u(1:(index_b3u - 1), :);
    end
    
    %% Find energy level spacings
    % Make vector of energy level spacings
    spacings_ag = zeros(size(elevels_ag, 1) - 1, 1);
    spacings_b1g = zeros(size(elevels_b1g, 1) - 1, 1);
    spacings_b2g = zeros(size(elevels_b2g, 1) - 1, 1);
    spacings_b3g = zeros(size(elevels_b3g, 1) - 1, 1);
    spacings_au = zeros(size(elevels_au, 1) - 1, 1);
    spacings_b1u = zeros(size(elevels_b1u, 1) - 1, 1);
    spacings_b2u = zeros(size(elevels_b2u, 1) - 1, 1);
    spacings_b3u = zeros(size(elevels_b3u, 1) - 1, 1);
    
    for i = 1:(size(elevels_ag, 1) - 1)
        spacings_ag(i) = abs(elevels_ag(i, 1) - elevels_ag(i+1, 1));
    end
    
    for i = 1:(size(elevels_b1g, 1) - 1)
        spacings_b1g(i) = abs(elevels_b1g(i, 1) - elevels_b1g(i+1, 1));
    end
    
    for i = 1:(size(elevels_b2g, 1) - 1)
        spacings_b2g(i) = abs(elevels_b2g(i, 1) - elevels_b2g(i+1, 1));
    end
    
    for i = 1:(size(elevels_b3g, 1) - 1)
        spacings_b3g(i) = abs(elevels_b3g(i, 1) - elevels_b3g(i+1, 1));
    end
    
    for i = 1:(size(elevels_au, 1) - 1)
        spacings_au(i) = abs(elevels_au(i, 1) - elevels_au(i+1, 1));
    end
    
    for i = 1:(size(elevels_b1u, 1) - 1)
        spacings_b1u(i) = abs(elevels_b1u(i, 1) - elevels_b1u(i+1, 1));
    end
    
    for i = 1:(size(elevels_b2u, 1) - 1)
        spacings_b2u(i) = abs(elevels_b2u(i, 1) - elevels_b2u(i+1, 1));
    end
    
    for i = 1:(size(elevels_b3u, 1) - 1)
        spacings_b3u(i) = abs(elevels_b3u(i, 1) - elevels_b3u(i+1, 1));
    end
    
    % Compute level spacing ratios
    r_ag = LSR(spacings_ag);
    r_b1g = LSR(spacings_b1g);
    r_b2g = LSR(spacings_b2g);
    r_b3g = LSR(spacings_b3g);
    r_au = LSR(spacings_au);
    r_b1u = LSR(spacings_b1u);
    r_b2u = LSR(spacings_b2u);
    r_b3u = LSR(spacings_b3u);

    % Fill in r array values
    r_array(index_r, :) = [N, r_ag, r_b1g, r_b2g, r_b3g, r_au, r_b1u, r_b2u, r_b3u];
   
    % Fill in size of irreps
    size_array(index_size, :) = [N, size(elevels_ag, 1), size(elevels_b1g, 1), ...
        size(elevels_b2g, 1), size(elevels_b3g, 1), size(elevels_au, 1), ...
        size(elevels_b1u, 1), size(elevels_b2u, 1), size(elevels_b3u, 1)];

    % Find proportion of solvable states
    solvable_prop(index_size, :) = [size(elevels_ag, 1), size(elevels_au, 1)] / (total_sites);

    %% Plot the level spacings histogram and show r values
    % Plot histograms
    figure(index_n)
    tiledlayout(2, 4,'TileSpacing', 'tight','Padding','Tight')
    bin_factor = 5;
    
    nexttile
    histogram(spacings_ag, ceil(size(elevels_ag, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('Ag')
    subtitle(['r = ' num2str(r_ag) '; total = ' num2str(size_array(index_size, 2))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_b1g, ceil(size(elevels_b1g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('B1g')
    subtitle(['r = ' num2str(r_b1g) '; total = ' num2str(size_array(index_size, 3))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_b2g, ceil(size(elevels_b2g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('B2g')
    subtitle(['r = ' num2str(r_b2g) '; total = ' num2str(size_array(index_size, 4))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_b3g, ceil(size(elevels_b3g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('B3g')
    subtitle(['r = ' num2str(r_b3g) '; total = ' num2str(size_array(index_size, 5))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_au, ceil(size(elevels_au, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('Au')
    subtitle(['r = ' num2str(r_au) '; total = ' num2str(size_array(index_size, 6))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_b1u, ceil(size(elevels_b1u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('B1u')
    subtitle(['r = ' num2str(r_b1u) '; total = ' num2str(size_array(index_size, 7))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_b2u, ceil(size(elevels_b2u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('B2u')
    subtitle(['r = ' num2str(r_b2u) '; total = ' num2str(size_array(index_size, 8))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_b3u, ceil(size(elevels_b3u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('B3u')
    subtitle(['r = ' num2str(r_b3u) '; total = ' num2str(size_array(index_size, 9))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    set(figure(index_n),'position',[0,100,25000,400])
    saveas(gcf, [folderName '/N=' num2str(N) '_elevel_hist.jpeg']);

    % Increment indices
    index_n = index_n + 1;
    index_size = index_size + 1;
    index_r = index_r + 1;
end

%% Plot r as a function of size in each irrep
figure(index_n)
tiledlayout(2,4,'TileSpacing', 'tight','Padding','Tight')

nexttile
plot(N_nums.', r_array(:, 2).', '-o', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Ag')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([0.3, 0.8])

nexttile
plot(N_nums, r_array(:, 3).', '-o', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('B1g')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([0.3, 0.8])

nexttile
plot(N_nums, r_array(:, 4).', '-o', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('B2g')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([0.3, 0.8])

nexttile
plot(N_nums, r_array(:, 5).', '-o', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('B3g')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([0.3, 0.8])

nexttile
plot(N_nums, r_array(:, 6).', '-o', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Au')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([0.3, 0.8])

nexttile
plot(N_nums, r_array(:, 7).', '-o', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('B1u')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([0.3, 0.8])

nexttile
plot(N_nums, r_array(:, 8).', '-o', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('B2u')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([0.3, 0.8])

nexttile
plot(N_nums, r_array(:, 9).', '-o', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('B3u')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([0.3, 0.8])

set(figure(index_n),'position',[0,100,3000,400])
saveas(gcf, [folderName '/LSR_plot.jpeg']);

%% Plot proportions of solvable states
figure(index_n + 1)
plot(N_nums, sum(solvable_prop, 2), '-o')
title('Proportion of states which are solvable')
xlabel('N')
yline(1/4, '--', '1/4')
ylim([0.248, 0.253])
ylabel('Proportion')
saveas(gcf, [folderName '/solvable_proportion.jpeg']);



%% Functions

% Gives x coordinate of index within each face
function x = xfacecoord(index, N, l1, l2, l3, total_sites)
    if (index <= l1 * (l2 + l3) * N^2) % In face 1 or 2
        x = mod(index - 1, l1 * N) + 1;
    elseif (index > total_sites - l1 * l2 * N^2) % In face 6
        x = mod(index - (total_sites - l1 * l2 * N^2) - 1, l1 * N) + 1;
    elseif (mod(index - l1 * (l2 + l3) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= l2 * N) % In face 3
        x = mod(index - l1 * (l2 + l3) * N^2 - 1, (2 * l2 + l1) * N) + 1;
    elseif (mod(index - l1 * (l2 + l3) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= (l2 + l1) * N) % In face 4
        x = mod(index - l1 * (l2 + l3) * N^2 - 1, (2 * l2 + l1) * N) + 1 - l2 * N;
    else % In face 5
        x = mod(index - l1 * (l2 + l3) * N^2 - 1, (2 * l2 + l1) * N) + 1 - (l2 + l1) * N;
    end 
end

% Gives y coordinate of index within each face
function y = yfacecoord(index, N, l1, l2, l3, total_sites) 
    if (index <= l1 * l3 * N^2) % in face 1
        y = ceil(index / (l1 * N));
    elseif (index <= l1 * (l2 + l3) * N^2) % in face 2 
        y = ceil((index - l1 * l3 * N^2) / (l1 * N));
    elseif (index <= total_sites - l1 * l2 * N^2) % in faces 3, 4, or 5
        y = ceil((index - l1 * (l2 + l3) * N^2) / ((2 * l2 + l1) * N));
    else % in face 6
        y = ceil((index - l1 * (l2 + l3) * N^2 - l3 * (2 * l2 + l1) * N^2) / (l1 * N));
    end
end

% From face coordinate and face number give index
function i = index_from_coord(x, y, face, N, l1, l2, l3, total_sites)
    if (face == 1)
        i = x + (l1 * N) * (y - 1);
    elseif (face == 2)
        i = (l1 * l3 * N^2) + x + (l1 * N) * (y - 1);
    elseif (face == 3)
        i = l1 * (l3 + l2) * N^2 + x + (2 * l2 + l1) * N * (y - 1);
    elseif (face == 4)
        i = l1 * (l3 + l2) * N^2 + l2 * N + x + (2 * l2 + l1) * N * (y - 1);
    elseif (face == 5)
        i = l1 * (l3 + l2) * N^2 + (l2 + l1) * N + x + (2 * l2 + l1) * N * (y - 1);
    else 
        i = total_sites - l1 * l2 * N^2 + x + l1 * N * (y - 1);
    end
end

% Given index on discretized cube returns new index after C2x symmetry
function i = c2x_index(index, N, l1, l2, l3, total_sites)
    if (index <= l1 * l3 * N^2) % in face 1
        new_face = 4;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l3 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (index <= l1 * (l3 + l2) * N^2) % in face 2 
        new_face = 2;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l2 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (index > total_sites - l1 * l2 * N^2) % in face 6
        new_face = 6;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l2 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (mod(index - l1 * (l3 + l2) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= l2 * N) % in face 3
        new_face = 5;
        new_x = xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (mod(index - l1 * (l3 + l2) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= (l2 + l1) * N) % in face 4
        new_face = 1;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l3 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    else  % in face 5
        new_face = 3;
        new_x = xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    end
end

% Given index on discretized cube returns new index after inv symmetry
function i = inv_index(index, N, l1, l2, l3, total_sites)
    if (index <= l1 * l3 * N^2) % in face 1
        new_face = 4;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (index <= l1 * (l3 + l2) * N^2) % in face 2  
        new_face = 6;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (index > total_sites - l1 * l2 * N^2) % in face 6
        new_face = 2;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (mod(index - l1 * (l3 + l2) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= l2 * N) % in face 3
        new_face = 5;
        new_x = xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l3 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (mod(index - l1 * (l3 + l2) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= (l2 + l1) * N) % in face 4
        new_face = 1;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    else  % in face 5
        new_face = 3;
        new_x = xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l3 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    end
end

% Given index on discretized cube returns new index after C2y symmetry
function i = c2y_index(index, N, l1, l2, l3, total_sites)
    if (index <= l1 * l3 * N^2) % in face 1
        new_face = 1;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l3 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (index <= l1 * (l3 + l2) * N^2) % in face 2 
        new_face = 6;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l2 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (index > total_sites - l1 * l2 * N^2) % in face 6
        new_face = 2;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l2 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (mod(index - l1 * (l3 + l2) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= l2 * N) % in face 3
        new_face = 5;
        new_x = l2 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l3 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (mod(index - l1 * (l3 + l2) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= (l2 + l1) * N) % in face 4
        new_face = 4;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l3 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    else  % in face 5
        new_face = 3;
        new_x = l2 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l3 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    end
end

% Given index on discretized cube returns new index after C2z symmetry
function i = c2z_index(index, N, l1, l2, l3, total_sites)
    if (index <= l1 * l3 * N^2) % in face 1
        new_face = 4;
        new_x = xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (index <= l1 * (l3 + l2) * N^2) % in face 2 
        new_face = 6;
        new_x = xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (index > total_sites - l1 * l2 * N^2) % in face 6
        new_face = 2;
        new_x = xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = yfacecoord(index, N, l1, l2, l3, total_sites);
        %fprintf(['new_x = ' num2str(new_x) ' new_y = ' num2str(new_y) '\n'])
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (mod(index - l1 * (l3 + l2) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= l2 * N) % in face 3
        new_face = 3;
        new_x = l2 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l3 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (mod(index - l1 * (l3 + l2) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= (l2 + l1) * N) % in face 4
        new_face = 1;
        new_x = xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    else  % in face 5
        new_face = 5;
        new_x = l2 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l3 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    end
end

% Functions to find indices of H matrix
function y = upper(x, N, l1, l2, l3, total_sites) 
    if (x > l1 * (l2 + l3) * N^2 - l1 * N && x <= l1 * (l2 + l3) * N^2) % x in top row of lower arm 
        y = l1 * (l2 + l3) * N^2 + l2 * N + mod(x - 1, l1 * N) + 1;
    elseif (x > l1 * (l2 + l3) * N^2 && x <= (total_sites - l2 * l1 * N^2 - (2 * l2 + l1) * N)) % x in middle arm but not top row
        y = x + (2 * l2 + l1) * N;
    elseif (x > total_sites - l1 * l2 * N^2 - (2 * l2 + l1) * N && x <= total_sites - l1 * l2 * N^2 - (l2 + l1) * N) % x top row of middle arms, left side
        dist_from_vertex = x - (total_sites - l1 * l2 * N^2 - (2 * l2 + l1) * N); % distance from leftmost vertex
        y = total_sites - dist_from_vertex * l1 * N + 1; 
    elseif (x > total_sites - l1 * l2 * N^2 - (2 * l2 + l1) * N && x <= total_sites - l1 * l2 * N^2 - l2 * N) % x in top row of middle arms, center side
        y = x + (l1 + l2) * N;
    elseif (x > total_sites - l1 * l2 * N^2 - (2 * l2 + l1) * N && x <= total_sites - l1 * l2 * N^2) % x in top row of middle arms, right side
        dist_from_fold = x - (total_sites - l1 * l2 * N^2 - l2 * N);
        y = total_sites - l1 * l2 * N^2 + dist_from_fold * l1 * N;
    elseif (x > total_sites - l1 * N) % x in top row of upper arm
        y = mod((x - total_sites - l1 * N) - 1, l1 * N) + 1;
    else % x in non-top row upper arm or non-top row lower arm
        y = x + l1 * N;
    end
end

function y = lower(x, N, l1, l2, l3, total_sites) 
    if (x > 0 && x <= l1 * N) % x in lower arm bottom row
        y = total_sites - l1 * N + mod(x - 1, l1 * N) + 1;
    elseif (x > l1 * (l3 + l2) * N^2 && x <= l1 * (l3 + l2) * N^2 + l2 * N) % x in middle arms, bottom row, left side
        dist_from_vertex = x - l1 * (l3 + l2) * N^2;
        y = l1 * l3 * N^2 + (dist_from_vertex - 1) * l1 * N + 1;
    elseif (x > l1 * (l3 + l2) * N^2 && x <= l1 * (l3 + l2) * N^2 + (l2 + l1) * N) % x in middle arms, bottom row, center
        y = x - (l2 + l1) * N;
    elseif (x > l1 * (l3 + l2) * N^2 && x <= l1 * (l3 + l2) * N^2 + (2 * l2 + l1) * N) % x in middle arms, bottom row, right side
        dist_from_fold = x - (l1 * (l3 + l2) * N^2 + (l2 + l1) * N);
        y = l1 * (l3 + l2) * N^2 - (dist_from_fold - 1) * l1 * N ;
    elseif (x > l1 * (l3 + l2) * N^2 + (2 * l2 + l1) * N && x <= total_sites - l1 * l2 * N^2) % x in middle arms, above bottom row
        y = x - (2 * l2 + l1) * N;
    elseif (x > total_sites - l1 * l2 * N^2 && x <= total_sites - l1 * l2 * N^2 + l1 * N) % x in upper arm, bottom row
        y = x - (l1 + l2) * N;
    else % x in upper arm, above bottom row or in lower arm above bottom row
        y = x - l1 * N;
    end
end

function y = right(x, N, l1, l2, l3, total_sites)
    if (x <= l3 * l1 * N^2 && floor(x / (l1 * N)) == x / (l1 * N)) % x in lower half of lower arm, rightmost column
        row_from_bottom = x / (l1 * N);
        y = total_sites - l2 * l1 * N^2 - (row_from_bottom - 1) * (2 * l2 + l1) * N;
    elseif (x > l1 * l3 * N^2 && x <= l1 * (l2 + l3) * N^2 && floor(x / (l1 * N)) == x / (l1 * N)) % x in upper half of lower arm, rightmost column
        row_from_top = (l1 * (l2 + l3) * N^2 - x) / (l1 * N) + 1;
        y = (l3 + l2) * l1 * N^2 + (l2 + l1) * N + row_from_top;
    elseif (x > (l3 + l2) * l1 * N^2 && x <= total_sites - l2 * l1 * N^2 && floor((x - (l2 + l3) * l1 * N^2)/((2 * l2 + l1) * N)) == (x - (l2 + l3) * l1 * N^2)/((2 * l2 + l1) * N)) % x in middle arm, rightmost column 
        row_from_top_2 = (total_sites - l2 * l1 * N^2 - x)/ ((2 * l2 + l1) * N) + 1;
        y = row_from_top_2 * l1 * N;
    elseif (x > total_sites - l2 * l1 * N^2 && floor((x - (total_sites - l2 * l1 * N^2)) / (l1 * N)) == (x - (total_sites - l2 * l1 * N^2)) / (l1 * N))  % x in upper arm, rightmost column 
        row_from_bottom_2 = (x - (total_sites - l2 * l1 * N^2)) / (l1 * N);
        y = total_sites - l2 * l1 * N^2 - l2 * N + row_from_bottom_2;
    else 
        y = x + 1;
    end
end

function y = left(x, N, l1, l2, l3, total_sites)
    if (x <= l3 * l1 * N^2 && (floor((x - 1) / (l1 * N)) == (x - 1) / (l1 * N))) % x in lower half of lower arm, leftmost column
        row_from_bottom = (x - 1) / (l1 * N) + 1;
        y = total_sites - l2 * l1 * N^2 - (row_from_bottom) * (2 * l2 + l1) * N + 1;
    elseif (x > l3 * l1 * N^2 && x <= (l3 + l2) * l1 * N^2 && floor((x - 1) / (l1 * N)) == (x - 1) / (l1 * N)) % x in upper half of lower arm, leftmost column
        row_from_bottom_2 = ceil(x / (l1 * N)) - l3 * N;
        y = (l3 + l2) * l1 * N^2 + row_from_bottom_2;
    elseif (x > (l3 + l2) * l1 * N^2 && x <= total_sites - l2 * l1 * N^2 && floor((x - 1 - (l3 + l2) * l1 * N^2)/((2 * l2 + l1) * N)) == (x - 1 - (l2 + l3) * l1 * N^2)/((2 * l2 + l1) * N)) % x in middle arm, leftmost column 
        row_from_top = ceil((total_sites - l2 * l1 * N^2 - x)/((2 * l2 + l1) * N));
        y = l1 * N * (row_from_top - 1) + 1;
    elseif (x > total_sites - l2 * l1 * N^2 && floor((x - (total_sites - l2 * l1 * N^2) - 1) / (l1 * N)) == (x - (total_sites - l2 * l1 * N^2) - 1) / (l1 * N))  % x in upper arm, leftmost column 
        row_from_top_2 = ceil((total_sites - x) / (l1 * N));
        y = total_sites - l2 * l1 * N^2 - (2 * l2 + l1) * N + row_from_top_2;
    else  
        y = x - 1;
    end
end

% Function to calculate mean level spacing ratio
function r = LSR(spacings)
    r_n = zeros(size(spacings, 1) - 1, 1);
    
    for i = 1:(size(spacings, 1) - 1)
        r_n(i) = min(spacings(i), spacings(i + 1)) ...
            / max(spacings(i), spacings(i + 1));
    end

    r = mean(r_n);
end