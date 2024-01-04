% Calculate r values for various N
tolerance = 1e-8;  % Tolerance for degenerate energy levels
potential = -0;
%N_nums = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53];
N_nums = [17];

r_array = zeros(size(N_nums, 1), 11);
size_array = zeros(size(N_nums, 1), 11);
solvable_prop = zeros(size(N_nums, 1), 4);

index_n = 1;
index_r = 1;
index_size = 1;

for n = 1:size(N_nums, 2)

    N = N_nums(n);
    fprintf(['\n    N = ' num2str(N) '\n'])
    total_sites = 10 * N^2; % Total number of sites in unfolded polyhedron

    %% Diagonalize H and sigma_d matrices
    tic;

    % Generate Hamiltonian matrix
    L = 1/N; % Spacing
    H = zeros(total_sites, total_sites); % Hamiltonian matrix
    
    for i = 1:total_sites
        H(i, i) = 4/L^2;
        H(i, upper(i, N)) = -1/L^2;
        H(i, lower(i, N)) = -1/L^2;
        H(i, right(i, N)) = -1/L^2;
        H(i, left(i, N)) = -1/L^2;
    end

    % Add corner potentials
    % Add potentials at all vertices
    for face = 1:6
        % Bottom left corner
        x = 1; 
        y = 1;
        index = index_from_coord(x, y, face, N);
        H(index, index) = H(index, index) + potential;
    
        % Top left corner
        x = 1; 
        y = N;
        index = index_from_coord(x, y, face, N);
        H(index, index) = H(index, index) + potential;
    
        % Bottom right corner
        x = N; 
        y = 1;
        index = index_from_coord(x, y, face, N);
        H(index, index) = H(index, index) + potential;
    
        % Top right corner
        x = N; 
        y= N;
        index = index_from_coord(x, y, face, N);
        H(index, index) = H(index, index) + potential;
    end
    
    % Generate matrix for sigma_d symmetry 
    sigma_d = zeros(total_sites, total_sites);
    for i = 1:(total_sites)
        sigma_d(i, sigma_d_index(i, N)) = 1;
    end
    
    %% Calculate eigenvectors and eigenvalues
    [eigenvectors, eigenvalues] = eig(H + ((1 + sqrt(5))/2) * sigma_d);
    fprintf('evals and evecs done \n')
    toc
    
    %% Generate Symmetry Matrices
    
    % Generate matrix for C2x
    C2x = zeros(total_sites, total_sites);
    for i = 1:(total_sites)
        C2x(i, c2x_index(i, N)) = 1;
    end

    % Generate matrix for i
    inv = zeros(total_sites, total_sites);
    for i = 1:(total_sites)
        inv(i, inv_index(i, N)) = 1;
    end

    % Generate matrix for C2y
    C2y = zeros(total_sites, total_sites);
    for i = 1:(total_sites)
        C2y(i, c2y_index(i, N)) = 1;
    end
    
    %% Calculate characters of symmetries
    tic;
    % Preallocate data array
    sigma_d_evals = zeros(total_sites, 1);
    c2x_evals = zeros(total_sites, 1);
    inv_evals = zeros(total_sites, 1);
    c2y_evals = zeros(total_sites, 1);
    
    % Find inv eigenvalues
    for i = 1:(total_sites)
        inv_evals(i) = (eigenvectors(:, i).' * inv * eigenvectors(:, i));
    end
    fprintf('inv done \n')
    
     % Find energy and sigma_d eigenvluaes
    for i = 1:(total_sites)
        sigma_d_evals(i) = (eigenvectors(:, i).' * sigma_d * eigenvectors(:, i));
        eigenvalues(i, i) = eigenvalues(i, i) - ((1 + sqrt(5))/2) * sigma_d_evals(i);
    end
    fprintf('sigma_d done \n')
    
    % Find c2x eigenvalues
    for i = 1:(total_sites)
        c2x_evals(i) = (eigenvectors(:, i).' * C2x * eigenvectors(:, i));
    end
    fprintf('c2x done \n')
    
    % Find c2y eigenvalues
    for i = 1:(total_sites)
        c2y_evals(i) = (eigenvectors(:, i).' * C2y * eigenvectors(:, i));
    end
    fprintf('c2y done \n')
    toc;
    
    %% Analyze Degeneracies
    % Reorder energies and evecs
    energy_levels = zeros(total_sites, 7);
    [energies, id] = sort(diag(eigenvalues));
    eigenvectors = eigenvectors(:, id);
    sigma_d_evals_sorted = sigma_d_evals(id);
    c2x_evals_sorted = c2x_evals(id);
    inv_evals_sorted = inv_evals(id);
    c2y_evals_sorted = c2y_evals(id);
    
    % Analyze for degeneracies
    energy = energies(1);
    degeneracy = 1;
    index = 1;
    trace_sigma_d = sigma_d_evals_sorted(1);
    trace_c2x = c2x_evals_sorted(1);
    trace_inv = inv_evals_sorted(1);
    trace_c2y = c2y_evals_sorted(1);
    
    for i = 2:(total_sites)
        if (abs(energies(i) - energy) < tolerance) % Next energy is degenerate
            degeneracy = degeneracy + 1;
            trace_sigma_d = trace_sigma_d + sigma_d_evals_sorted(i);
            trace_c2x = trace_c2x + c2x_evals_sorted(i);
            trace_inv = trace_inv + inv_evals_sorted(i);
            trace_c2y = trace_c2y + c2y_evals_sorted(i);
        else % Next energy is new
            % Record stats for previous energy level
            energy_levels(index, 1) = energy;
            energy_levels(index, 2) = degeneracy;
            energy_levels(index, 3) = trace_c2x;
            energy_levels(index, 4) = trace_c2y;
            energy_levels(index, 5) = trace_inv;
            energy_levels(index, 6) = trace_sigma_d;
            energy_levels(index, 7) = i - 1;
    
            % Reset variables for new energy level
            energy = energies(i);
            degeneracy = 1;
            index = index + 1;
            trace_c2x = c2x_evals_sorted(i);
            trace_c2y = c2y_evals_sorted(i);
            trace_inv = inv_evals_sorted(i);
            trace_sigma_d = sigma_d_evals_sorted(i);
            
        end
    
        % Record energy level if we reach the end
        if (i == total_sites)
            energy_levels(index, 1) = energy;
            energy_levels(index, 2) = degeneracy;
            energy_levels(index, 3) = trace_c2x;
            energy_levels(index, 4) = trace_c2y;
            energy_levels(index, 5) = trace_inv;
            energy_levels(index, 6) = trace_sigma_d;
            energy_levels(index, 7) = i;
            
        end
    end
    
    energy_levels = energy_levels(1:index, :);
    
    %% Separate the representation classes
    % Allocate space for each irrep
    elevels_a1g = zeros(total_sites, 7);
    elevels_a2g = zeros(total_sites, 7);
    elevels_b1g = zeros(total_sites, 7);
    elevels_b2g = zeros(total_sites, 7);
    elevels_eg = zeros(total_sites, 7);
    elevels_a1u = zeros(total_sites, 7);
    elevels_a2u = zeros(total_sites, 7);
    elevels_b1u = zeros(total_sites, 7);
    elevels_b2u = zeros(total_sites, 7);
    elevels_eu = zeros(total_sites, 7);
    
    index_a1g = 1;
    index_a2g = 1;
    index_b1g = 1;
    index_b2g = 1;
    index_eg = 1;
    index_a1u = 1;
    index_a2u = 1;
    index_b1u = 1;
    index_b2u = 1;
    index_eu = 1;
    
    % Round the energy_levels characters
    energy_levels_rounded = zeros(size(energy_levels, 1), 7);
    energy_levels_rounded(:, 2:6) = round(energy_levels(:, 2:6));
    energy_levels_rounded(:, 1) = energy_levels(:, 1);
    energy_levels_rounded(:, 7) = energy_levels(:, 7);
    
    % Fill in elevel info in each irrep
    for i = 1:size(energy_levels_rounded, 1)
        trace_e = energy_levels_rounded(i, 2);
        trace_c2x = energy_levels_rounded(i, 3);
        trace_c2y = energy_levels_rounded(i, 4);
        trace_inv = energy_levels_rounded(i, 5);
        trace_sigma_d = energy_levels_rounded(i, 6);
    
        traces = [trace_e, trace_c2x, trace_c2y, trace_inv, trace_sigma_d];
    
        if (isequal(traces, [1, 1, 1, 1, 1])) % A1g
            elevels_a1g(index_a1g, :) = energy_levels_rounded(i, :);
            index_a1g = index_a1g + 1;
        elseif (isequal(traces, [2, 2, 2, 0, 0])) % Accidental degeneracy A1g + A1u
            elevels_a1g(index_a1g, :) = energy_levels_rounded(i, :);
            index_a1g = index_a1g + 1;
            elevels_a1u(index_a1u, :) = energy_levels_rounded(i, :); 
            index_a1u = index_a1u + 1;
        elseif (isequal(traces, [2, -2, -2, 0, 0])) % Accidental degeneracy A2g + A2u
            elevels_a2g(index_a2g, :) = energy_levels_rounded(i, :);
            index_a2g = index_a2g + 1;
            elevels_a2u(index_a2u, :) = energy_levels_rounded(i, :);
            index_a2u = index_a2u + 1;
        elseif (isequal(traces, [2, 2, -2, 0, 0])) % Accidental degeneracy B1g + B1u
            elevels_b1g(index_b1g, :) = energy_levels_rounded(i, :);
            index_b1g = index_b1g + 1;
            elevels_b1u(index_b1u, :) = energy_levels_rounded(i, :); 
            index_b1u = index_b1u + 1;
        elseif (isequal(traces, [2, -2, 2, 0, 0])) % Accidental degeneracy B2g + B2u
            elevels_b2g(index_b2g, :) = energy_levels_rounded(i, :);
            index_b2g = index_b2g + 1;
            elevels_b2u(index_b2u, :) = energy_levels_rounded(i, :);
            index_b2u = index_b2u + 1;
        elseif (isequal(traces, [1, -1, -1, 1, -1])) % A2g
            elevels_a2g(index_a2g, :) = energy_levels_rounded(i, :);
            index_a2g = index_a2g + 1;
        elseif (isequal(traces, [1, 1, -1, 1, -1])) % B1g
            elevels_b1g(index_b1g, :) = energy_levels_rounded(i, :);
            index_b1g = index_b1g + 1;
        elseif (isequal(traces, [1, -1, 1, 1, 1])) % B2g
            elevels_b2g(index_b2g, :) = energy_levels_rounded(i, :);
            index_b2g = index_b2g + 1;
        elseif (isequal(traces, [2, 0, 0, 2, 0])) % Eg
            elevels_eg(index_eg, :) = energy_levels_rounded(i, :);
            index_eg = index_eg + 1;
        elseif (isequal(traces, [1, 1, 1, -1, -1])) % A1u
            elevels_a1u(index_a1u, :) = energy_levels_rounded(i, :);
            index_a1u = index_a1u + 1;
        elseif (isequal(traces, [1, -1, -1, -1, 1])) % A2u
            elevels_a2u(index_a2u, :) = energy_levels_rounded(i, :);
            index_a2u = index_a2u + 1;
        elseif (isequal(traces, [1, 1, -1, -1, 1])) % B1u
            elevels_b1u(index_b1u, :) = energy_levels_rounded(i, :);
            index_b1u = index_b1u + 1;
        elseif (isequal(traces, [1, -1, 1, -1, -1])) % B2u
            elevels_b2u(index_b2u, :) = energy_levels_rounded(i, :);
            index_b2u = index_b2u + 1;
        elseif (isequal(traces, [2, 0, 0, -2, 0])) % Eu
            elevels_eu(index_eu, :) = energy_levels_rounded(i, :);
            index_eu = index_eu + 1;
        else 
            fprintf([num2str(energy_levels_rounded(i, 1)) ' ' num2str(energy_levels_rounded(i, 2)) ' ' ...
                num2str(energy_levels_rounded(i, 3)) ' ' num2str(energy_levels_rounded(i, 4)) ' ' ...
                num2str(energy_levels_rounded(i, 5)) ' ' num2str(energy_levels_rounded(i, 6)) ' ' ...
                num2str(energy_levels_rounded(i, 7)) '\n'])
        end
    end
    
    % Remove extra rows of zeros
    if (index_a1g > 1)
        elevels_a1g = elevels_a1g(1:(index_a1g - 1), :);
    end 
    
    if (index_a2g > 1)
        elevels_a2g = elevels_a2g(1:(index_a2g - 1), :);
    end
    
    if (index_b1g > 1)
        elevels_b1g = elevels_b1g(1:(index_b1g - 1), :);
    end
    
    if (index_b2g > 1)
        elevels_b2g = elevels_b2g(1:(index_b2g - 1), :);
    end
    
    if (index_eg > 1)
        elevels_eg = elevels_eg(1:(index_eg - 1), :);
    end
    
    if (index_a1u > 1)
        elevels_a1u = elevels_a1u(1:(index_a1u - 1), :);
    end 
    
    if (index_a2u > 1)
        elevels_a2u = elevels_a2u(1:(index_a2u - 1), :);
    end
    
    if (index_b1u > 1)
        elevels_b1u = elevels_b1u(1:(index_b1u - 1), :);
    end
    
    if (index_b2u > 1)
        elevels_b2u = elevels_b2u(1:(index_b2u - 1), :);
    end
    
    if (index_eu > 1)
        elevels_eu = elevels_eu(1:(index_eu - 1), :);
    end
    
    %% Find energy level spacings
    % Make vector of energy level spacings
    spacings_a1g = zeros(size(elevels_a1g, 1) - 1, 1);
    spacings_a2g = zeros(size(elevels_a2g, 1) - 1, 1);
    spacings_b1g = zeros(size(elevels_b1g, 1) - 1, 1);
    spacings_b2g = zeros(size(elevels_b2g, 1) - 1, 1);
    spacings_eg = zeros(size(elevels_eg, 1) - 1, 1);
    spacings_a1u = zeros(size(elevels_a1u, 1) - 1, 1);
    spacings_a2u = zeros(size(elevels_a2u, 1) - 1, 1);
    spacings_b1u = zeros(size(elevels_b1u, 1) - 1, 1);
    spacings_b2u = zeros(size(elevels_b2u, 1) - 1, 1);
    spacings_eu = zeros(size(elevels_eu, 1) - 1, 1);
    
    for i = 1:(size(elevels_a1g, 1) - 1)
        spacings_a1g(i) = abs(elevels_a1g(i, 1) - elevels_a1g(i+1, 1));
    end
    
    for i = 1:(size(elevels_a2g, 1) - 1)
        spacings_a2g(i) = abs(elevels_a2g(i, 1) - elevels_a2g(i+1, 1));
    end
    
    for i = 1:(size(elevels_b1g, 1) - 1)
        spacings_b1g(i) = abs(elevels_b1g(i, 1) - elevels_b1g(i+1, 1));
    end
    
    for i = 1:(size(elevels_b2g, 1) - 1)
        spacings_b2g(i) = abs(elevels_b2g(i, 1) - elevels_b2g(i+1, 1));
    end
    
    for i = 1:(size(elevels_eg, 1) - 1)
        spacings_eg(i) = abs(elevels_eg(i, 1) - elevels_eg(i+1, 1));
    end
    
    for i = 1:(size(elevels_a1u, 1) - 1)
        spacings_a1u(i) = abs(elevels_a1u(i, 1) - elevels_a1u(i+1, 1));
    end
    
    for i = 1:(size(elevels_a2u, 1) - 1)
        spacings_a2u(i) = abs(elevels_a2u(i, 1) - elevels_a2u(i+1, 1));
    end
    
    for i = 1:(size(elevels_b1u, 1) - 1)
        spacings_b1u(i) = abs(elevels_b1u(i, 1) - elevels_b1u(i+1, 1));
    end
    
    for i = 1:(size(elevels_b2u, 1) - 1)
        spacings_b2u(i) = abs(elevels_b2u(i, 1) - elevels_b2u(i+1, 1));
    end
    
    for i = 1:(size(elevels_eu, 1) - 1)
        spacings_eu(i) = abs(elevels_eu(i, 1) - elevels_eu(i+1, 1));
    end
    
    % Compute level spacing ratios
    r_a1g = LSR(spacings_a1g);
    r_a2g = LSR(spacings_a2g);
    r_b1g = LSR(spacings_b1g);
    r_b2g = LSR(spacings_b2g);
    r_eg = LSR(spacings_eg);
    r_a1u = LSR(spacings_a1u);
    r_a2u = LSR(spacings_a2u);
    r_b1u = LSR(spacings_b1u);
    r_b2u = LSR(spacings_b2u);
    r_eu = LSR(spacings_eu);

    % Fill in r array values
    r_array(index_r, :) = [N, r_a1g, r_a2g, r_b1g, r_b2g, r_eg, r_a1u, r_a2u, r_b1u, r_b2u, r_eu];
   
    % Fill in size of irreps
    size_array(index_size, :) = [N, size(elevels_a1g, 1), size(elevels_a2g, 1), size(elevels_b1g, 1), ...
        size(elevels_b2g, 1), size(elevels_eg, 1), size(elevels_a1u, 1), size(elevels_a2u, 1), ...
        size(elevels_b1u, 1), size(elevels_b2u, 1), size(elevels_eu, 1)];

    % Find proportion of solvable states
    solvable_prop(index_size, :) = [size(elevels_a1g, 1), size(elevels_a2g, 1), size(elevels_a1u, 1), size(elevels_a2u, 1)] / (total_sites);

    %% Plot the level spacings histogram and show r values
    % Plot histograms
    figure(index_n)
    tiledlayout(2,5)
    bin_factor = 5;
    
    nexttile
    histogram(spacings_a1g, ceil(size(elevels_a1g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('A1g')
    subtitle(['r = ' num2str(r_a1g) '; total = ' num2str(size_array(index_size, 2))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_a2g, ceil(size(elevels_a2g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('A2g')
    subtitle(['r = ' num2str(r_a2g) '; total = ' num2str(size_array(index_size, 3))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_b1g, ceil(size(elevels_b1g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('B1g')
    subtitle(['r = ' num2str(r_b1g) '; total = ' num2str(size_array(index_size, 4))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_b2g, ceil(size(elevels_b2g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('B2g')
    subtitle(['r = ' num2str(r_b2g) '; total = ' num2str(size_array(index_size, 5))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_eg, ceil(size(elevels_eg, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('Eg')
    subtitle(['r = ' num2str(r_eg) '; total = ' num2str(size_array(index_size, 6))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_a1u, ceil(size(elevels_a1u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('A1u')
    subtitle(['r = ' num2str(r_a1u) '; total = ' num2str(size_array(index_size, 7))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_a2u, ceil(size(elevels_a2u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('A2u')
    subtitle(['r = ' num2str(r_a2u) '; total = ' num2str(size_array(index_size, 8))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_b1u, ceil(size(elevels_b1u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('B1u')
    subtitle(['r = ' num2str(r_b1u) '; total = ' num2str(size_array(index_size, 9))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_b2u, ceil(size(elevels_b2u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('B2u')
    subtitle(['r = ' num2str(r_b2u) '; total = ' num2str(size_array(index_size, 10))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_eu, ceil(size(elevels_eu, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('Eu')
    subtitle(['r = ' num2str(r_eu) '; total = ' num2str(size_array(index_size, 11))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    set(figure(index_n),'position',[0,100,3000,400])
    
    % Incremenet indices
    index_n = index_n + 1;
    index_size = index_size + 1;
    index_r = index_r + 1;
end

% Plot r as a function of size in each irrep
figure(index_n)
tiledlayout(2,5)

nexttile
plot(N_nums.', r_array(:, 2).', '-o', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('A1g')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([0.3, 1])

nexttile
plot(N_nums, r_array(:, 3).', '-o', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('A2g')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([0.3, 1])

nexttile
plot(N_nums, r_array(:, 4).', '-o', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('B1g')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([0.3, 1])

nexttile
plot(N_nums, r_array(:, 5).', '-o', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('B2g')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([0.3, 1])

nexttile
plot(N_nums, r_array(:, 6).', '-o', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Eg')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([0.3, 1])

nexttile
plot(N_nums, r_array(:, 7).', '-o', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('A1u')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([0.3, 1])

nexttile
plot(N_nums, r_array(:, 8).', '-o', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('A2u')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([0.3, 1])

nexttile
plot(N_nums, r_array(:, 9).', '-o', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('B1u')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([0.3, 1])

nexttile
plot(N_nums, r_array(:, 10).', '-o', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('B2u')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([0.3, 1])

nexttile
plot(N_nums, r_array(:, 11).', '-o', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Eu')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([0.3, 1])

set(figure(index_n),'position',[0,100,3000,400])

% Plot proportions of solvable states
figure(index_n + 1)
plot(N_nums, sum(solvable_prop, 2), '-o')
title('Proportion of states which are solvable')
xlabel('N')
yline(1/12, '--', '1/12')
ylim([0, 0.12])
ylabel('Proportion')


%% Functions

% Gives x coordinate of index within each face
function x = xfacecoord(index, N)
    if (index <= 3 * N^2) % In face 1 or 2
        x = mod(index - 1, N) + 1;
    elseif (index > 8 * N^2) % In face 6
        x = mod(index - 1, N) + 1;
    elseif (mod(index - 3 * N^2 - 1, 5 * N) + 1 <= 2 * N) % In face 3
        x = mod(index - 3 * N^2 - 1, 5 * N) + 1;
    elseif (mod(index - 3 * N^2 - 1, 5 * N) + 1 <= 3 * N) % In face 4
        x = mod(index - 3 * N^2 - 1, 5 * N) + 1 - 2 * N;
    else % In face 5
        x = mod(index - 3 * N^2 - 1, 5 * N) + 1 - 3 * N;
    end 
end

% Gives y coordinate of index within each face
function y = yfacecoord(index, N) 
    if (index <= N^2) % in face 1
        y = ceil(index / N);
    elseif (index <= 3 * N^2) % in face 2 
        y = ceil((index - N^2) / N);
    elseif (index <= 8 * N^2) % in faces 3, 4, or 5
        y = ceil((index - 3 * N^2) / (5 * N));
    else % in face 6
        y = ceil((index - 8 * N^2) / N);
    end
end

% From face coordinate and face number give index
function i = index_from_coord(x, y, face, N)
    if (face == 1)
        i = x + N * (y - 1);
    elseif (face == 2)
        i = N^2 + x + N * (y - 1);
    elseif (face == 3)
        i = 3 * N^2 + x + 5 * N * (y - 1);
    elseif (face == 4)
        i = 3 * N^2 + 2 * N + x + 5 * N * (y - 1);
    elseif (face == 5)
        i = 3 * N^2 + 3 * N + x + 5 * N * (y - 1);
    else 
        i = 8 * N^2 + x + N * (y - 1);
    end
end

% Given index on discretized cube returns new index after sigma_d symmetry
function i = sigma_h_index(index, N)
    if (index <= N^2) % in face 1
        new_face = 4;
        new_x = xfacecoord(index, N);
        new_y = yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (index <= 3 * N^2) % in face 2 
        new_face = 2;
        new_x = xfacecoord(index, N);
        new_y = 2 * N + 1 - yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (index > 8 * N^2) % in face 6
        new_face = 6;
        new_x = xfacecoord(index, N);
        new_y = 2 * N + 1 - yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (mod(index - 3 * N^2 - 1, 5 * N) + 1 <= 2 * N) % in face 3
        new_face = 3;
        new_x = 2 * N + 1 - xfacecoord(index, N);
        new_y = yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (mod(index - 3 * N^2 - 1, 5 * N) + 1 <= 3 * N) % in face 4
        new_face = 1;
        new_x = xfacecoord(index, N);
        new_y = yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    else  % in face 5
        new_face = 5;
        new_x = 2 * N + 1 - facecoord(index, N);
        new_y = yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    end
end

% Given index on discretized cube returns new index after sigma_d symmetry
function i = sigma_d_index(index, N)
    if (index <= N^2) % in face 1
        new_face = 1;
        new_x = N + 1 - yfacecoord(index, N);
        new_y = N + 1 - xfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (index <= 3 * N^2) % in face 2 
        new_face = 3;
        new_x = yfacecoord(index, N);
        new_y = xfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (index > 8 * N^2) % in face 6
        new_face = 5;
        new_x = yfacecoord(index, N);
        new_y = xfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (mod(index - 3 * N^2 - 1, 5 * N) + 1 <= 2 * N) % in face 3
        new_face = 2;
        new_x = yfacecoord(index, N);
        new_y = xfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (mod(index - 3 * N^2 - 1, 5 * N) + 1 <= 3 * N) % in face 4
        new_face = 4;
        new_x = yfacecoord(index, N);
        new_y = xfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    else  % in face 5
        new_face = 6;
        new_x = yfacecoord(index, N);
        new_y = xfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    end
end

% Given index on discretized cube returns new index after C2x symmetry
function i = c2x_index(index, N)
    if (index <= N^2) % in face 1
        new_face = 4;
        new_x = N + 1 - xfacecoord(index, N);
        new_y = N + 1 - yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (index <= 3 * N^2) % in face 2 
        new_face = 2;
        new_x = N + 1 - xfacecoord(index, N);
        new_y = 2 * N + 1 - yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (index > 8 * N^2) % in face 6
        new_face = 6;
        new_x = N + 1 - xfacecoord(index, N);
        new_y = 2 * N + 1 - yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (mod(index - 3 * N^2 - 1, 5 * N) + 1 <= 2 * N) % in face 3
        new_face = 5;
        new_x = xfacecoord(index, N);
        new_y = yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (mod(index - 3 * N^2 - 1, 5 * N) + 1 <= 3 * N) % in face 4
        new_face = 1;
        new_x = N + 1 - xfacecoord(index, N);
        new_y = N + 1 - yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    else  % in face 5
        new_face = 3;
        new_x = xfacecoord(index, N);
        new_y = yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    end
end

% Given index on discretized cube returns new index after sigma_v symmetry
function i = inv_index(index, N)
    if (index <= N^2) % in face 1
        new_face = 4;
        new_x = N + 1 - xfacecoord(index, N);
        new_y = yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (index <= 3 * N^2) % in face 2 
        new_face = 6;
        new_x = N + 1 - xfacecoord(index, N);
        new_y = yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (index > 8 * N^2) % in face 6
        new_face = 2;
        new_x = N + 1 - xfacecoord(index, N);
        new_y = yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (mod(index - 3 * N^2 - 1, 5 * N) + 1 <= 2 * N) % in face 3
        new_face = 5;
        new_x = xfacecoord(index, N);
        new_y = N + 1 - yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (mod(index - 3 * N^2 - 1, 5 * N) + 1 <= 3 * N) % in face 4
        new_face = 1;
        new_x = N + 1 - xfacecoord(index, N);
        new_y = yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    else  % in face 5
        new_face = 3;
        new_x = xfacecoord(index, N);
        new_y = N + 1 - yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    end
end

% Given index on discretized cube returns new index after C2y symmetry
function i = c2y_index(index, N)
    if (index <= N^2) % in face 1
        new_face = 4;
        new_x = yfacecoord(index, N);
        new_y = N + 1 - xfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (index <= 3 * N^2) % in face 2 
        new_face = 5;
        new_x = yfacecoord(index, N);
        new_y = N + 1 - xfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (index > 8 * N^2) % in face 6
        new_face = 3;
        new_x = yfacecoord(index, N);
        new_y = N + 1 - xfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (mod(index - 3 * N^2 - 1, 5 * N) + 1 <= 2 * N) % in face 3
        new_face = 6;
        new_x = N + 1 - yfacecoord(index, N);
        new_y = xfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (mod(index - 3 * N^2 - 1, 5 * N) + 1 <= 3 * N) % in face 4
        new_face = 1;
        new_x = N + 1 - yfacecoord(index, N);
        new_y = xfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    else  % in face 5
        new_face = 2;
        new_x = N + 1 - yfacecoord(index, N);
        new_y = xfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    end
end
% Functions to find indices of H matrix
function y = upper(x, N) 
    if (x > 3 * N^2 - N && x <= 3 * N^2) % x in top row of lower arm 
        y = 3 * N^2 + 2 * N + mod(x - 1, N) + 1;
    elseif (x > 3 * N^2 && x <= (8 * N^2 - 5 * N)) % x in middle arm but not top row
        y = x + 5 * N;
    elseif (x > 8 * N^2 - 5 * N && x <= 8 * N^2 - 3 * N) % x top row of middle arms, left side
        dist_from_vertex = x - (8 * N^2 - 5 * N); % distance from leftmost vertex
        y = 10 * N^2 - dist_from_vertex * N + 1; 
    elseif (x > 8 * N^2 - 3 * N && x <= 8 * N^2 - 2 * N) % x in top row of middle arms, center side
        y = 8 * N^2 + mod(x - 1, N) + 1;
    elseif (x > 8 * N^2 - 2 * N && x <= 8 * N^2) % x in top row of middle arms, right side
        dist_from_fold = x - (8 * N^2 - 2 * N);
        y = 8 * N^2 + dist_from_fold * N;
    elseif (x > 8 * N^2 + N * (2 * N - 1)) % x in top row of upper arm
        y = mod(x - 1, N) + 1;
    else % x in non-top row upper arm or non-top row lower arm
        y = x + N;
    end
end

function y = lower(x, N) 
    if (x > 0 && x <= N) % x in lower arm bottom row
        y = 10 * N^2 - N + mod(x - 1, N) + 1;
    elseif (x > 3 * N^2 && x <= 3 * N^2 + 2 * N) % x in middle arms, bottom row, left side
        dist_from_vertex = x - 3 * N^2;
        y = N^2 + (dist_from_vertex - 1) * N + 1;
    elseif (x > 3 * N^2 + 2 * N && x <= 3 * N^2 + 3 * N) % x in middle arms, bottom row, center
        y = x - 3 * N;
    elseif (x > 3 * N^2 + 3 * N && x <= 3 * N^2 + 5 * N) % x in middle arms, bottom row, right side
        dist_from_fold = x - (3 * N^2 + 3 * N);
        y = 3 * N^2 - (dist_from_fold - 1) * N ;
    elseif (x > 3 * N^2 + 5 * N && x <= 8 * N^2) % x in middle arms, above bottom row
        y = x - 5 * N;
    elseif (x > 8 * N^2 && x <= 8 * N^2 + N) % x in upper arm, bottom row
        y = x - 3 * N;
    else % x in upper arm, above bottom row or in lower arm above bottom row
        y = x - N;
    end
end

function y = right(x, N)
    if (x <= N^2 && floor(x / N) == x / N) % x in lower half of lower arm, rightmost column
        row_from_bottom = x / N;
        y = 8 * N^2 - (row_from_bottom - 1) * 5 * N;
    elseif (x > N^2 && x <= 3 * N^2 && floor(x / N) == x / N) % x in upper half of lower arm, rightmost column
        row_from_top = (3 * N^2 - x) / N + 1;
        y = 3 * N^2 + 3 * N + row_from_top;
    elseif (x > 3 * N^2 && x <= 8 * N^2 && floor((x - 3 * N^2)/(5 * N)) == (x - 3 * N^2)/(5 * N)) % x in middle arm, rightmost column 
        row_from_top_2 = (8 * N^2 - x)/ (5 * N) + 1;
        y = row_from_top_2 * N;
    elseif (x > 8 * N^2 && floor(x / N) == x / N)  % x in upper arm, rightmost column 
        row_from_bottom_2 = (x - 8 * N^2) / N;
        y = 8 * N^2 - 2 * N + row_from_bottom_2;
    else  
        y = x + 1;
    end
end

function y = left(x, N)
    if (x <= N^2 && floor((x - 1) / N) == (x - 1) / N) % x in lower half of lower arm, leftmost column
        row_from_bottom = (x - 1) / N + 1;
        y = 8 * N^2 - (row_from_bottom) * 5 * N + 1;
    elseif (x > N^2 && x <= 3 * N^2 && floor((x - 1) / N) == (x - 1) / N) % x in upper half of lower arm, leftmost column
        row_from_bottom_2 = ceil(x / N) - N;
        y = 3 * N^2 + row_from_bottom_2;
    elseif (x > 3 * N^2 && x <= 8 * N^2 && floor((x - 1 - 3 * N^2)/(5 * N)) == (x - 1 - 3 * N^2)/(5 * N)) % x in middle arm, leftmost column 
        row_from_top = ceil((8 * N^2 - x)/(5 * N));
        y = N * (row_from_top - 1) + 1;
    elseif (x > 8 * N^2 && floor((x - 1) / N) == (x - 1) / N)  % x in upper arm, leftmost column 
        row_from_top_2 = ceil((10 * N^2 - x) / N);
        y = 8 * N^2 - 5 * N + row_from_top_2;
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