% Calculate r values for various N
clear;
close all
tolerance = 1e-8;  % Tolerance for degenerate energy levels
corner_potential = 1;

% Find max N
upper_lim = 30000; % Max number of sites allowed in laptop memory
max_N = floor(sqrt(upper_lim / 2) - 1);
fprintf(['Max N = ' num2str(max_N) '\n'])
N_nums = primes(max_N);
N_nums = N_nums(5:end);
N_nums = [7];

% Make directory
folderName = 'octahedron_plots';
mkdir(folderName);

r_array = zeros(size(N_nums, 1), 11);
size_array = zeros(size(N_nums, 1), 11);
solvable_prop = zeros(size(N_nums, 1), 4);

index_n = 1;
index_r = 1;
index_size = 1;


for n = 1:size(N_nums, 2)
    N = N_nums(n);

    %% Diagonalize H and sigma_d matrices
    tic;
    fprintf(['\n    N = ' num2str(N) '\n'])

    % Generate Hamiltonian matrix
    L = 1/N; % Spacing
    H = zeros(2 * (N + 1)^2, 2 * (N + 1)^2); % Hamiltonian matrix
    
    for i = 1:2 * (N + 1)^2
        %fprintf([num2str(i) ' '])
        H(i, i) = 4/L^2;
        %fprintf([num2str(right(i, N)) ' '])
        H(i, right(i, N)) = -(1 + 2/sqrt(3))/L^2;
        %fprintf([num2str(left(i, N)) ' '])
        H(i, left(i, N)) = -(1 + 2/sqrt(3))/L^2;
        if (has_upper(i, N))
            %fprintf('upper ')
        end
        if (has_upper(i, N))
            %fprintf([num2str(upper(i, N)) ' \n'])
            H(i, upper(i, N)) = -(1 + 2/sqrt(3))/L^2;
        else
            fprintf([num2str(i) ' ' num2str(lower(i, N)) ' \n'])
            H(i, lower(i, N)) = -(1 + 2/sqrt(3))/L^2;
        end
    end

    % Add corner potentials
    faces = [3, 5, 7, 8];
    for i = 1:4
        face = faces(i);

        lower_left_corner = index_from_coord(1, 1, face, N); 
        lower_right_corner = index_from_coord(N, 1, face, N); 
        top_corner = index_from_coord(1, (N + 1) / 2, face, N);

        %fprintf([num2str(lower_left_corner) ' ' num2str(lower_right_corner) ' ' num2str(top_corner) '\n'])

        H(lower_left_corner, lower_left_corner) = corner_potential; 
        H(lower_right_corner, lower_right_corner) = corner_potential;
        H(top_corner, top_corner) = corner_potential;
    end

    % Face 2
    faces = [1, 2, 4, 6];
    for i = 1:4
        face = faces(i);

        bottom_corner = index_from_coord(1, 1, face, N);
        top_right_corner = index_from_coord(N, (N + 1) / 2, face, N);
        top_left_corner = index_from_coord(1, (N + 1) / 2, face, N);

        H(bottom_corner, bottom_corner) = corner_potential; 
        H(top_right_corner, top_right_corner) = corner_potential;
        H(top_left_corner, top_left_corner) = corner_potential;
    end
    
    % Generate matrix for sigma_d symmetry 
    sigma_d = zeros(2 * (N + 1)^2, 2 * (N + 1)^2);
    for i = 1: (2 * (N + 1)^2)
        sigma_d(i, sigma_d_index(i, N)) = 1;
    end
    
    % Calculate eigenvectors and eigenvalues
    [eigenvectors, eigenvalues] = eig(H + ((1 + sqrt(5))/2)*sigma_d);
    fprintf('evals and evecs done \n')
    toc
    
    %% Generate Symmetry Matrices
    
    % Generate matrix for C3
    C3 = zeros(2 * (N + 1)^2, 2 * (N + 1)^2);
    
    for i = 1: (2 * (N + 1)^2)
        C3(i, c3_index(i, N)) = 1;
    end
    
    % Generate matrix for i
    inv = zeros(2 * (N + 1)^2, 2 * (N + 1)^2);
    for i = 1: (2 * (N + 1)^2)
        inv(i, inv_index(i, N)) = 1;
    end
    
    % Generate matrix for C2
    C2 = zeros(2 * (N + 1)^2, 2 * (N + 1)^2);
    for i = 1: (2 * (N + 1)^2)
        C2(i, c2_index(i, N)) = 1;
    end
    
    %% Calculate characters of symmetries
    tic;
    % Preallocate data array
    sigma2_evals = zeros(2 * (N + 1)^2, 1);
    c3_evals = zeros(2 * (N + 1)^2, 1);
    inv_evals = zeros(2 * (N + 1)^2, 1);
    c2_evals = zeros(2 * (N + 1)^2, 1);
    
    % Subtract off sigma2 eigenvalues from eigenvalues 
    for i = 1:(2 * (N + 1)^2)
        inv_evals(i) = (eigenvectors(:, i).' * inv * eigenvectors(:, i));
    end
    fprintf('inv done \n')
    
    for i = 1:(2 * (N + 1)^2)
        sigma2_evals(i) = (eigenvectors(:, i).' * sigma_d * eigenvectors(:, i));
        eigenvalues(i, i) = eigenvalues(i, i) - ((1 + sqrt(5))/2) * sigma2_evals(i);
    end
    fprintf('sigma_d done \n')
    
    for i = 1:(2 * (N + 1)^2)
        c3_evals(i) = (eigenvectors(:, i).' * C3 * eigenvectors(:, i));
    end
    fprintf('c3 done \n')
    
    for i = 1:(2 * (N + 1)^2)
        c2_evals(i) = (eigenvectors(:, i).' * C2 * eigenvectors(:, i));
    end
    fprintf('c2 done \n')
    toc;
    
    %% Analyze Degeneracies
    % Reorder energies and evecs
    energy_levels = zeros(2 * (N + 1)^2, 7);
    [energies, id] = sort(diag(eigenvalues));
    eigenvectors = eigenvectors(:, id);
    sigma2_evals_sorted = sigma2_evals(id);
    c3_evals_sorted = c3_evals(id);
    inv_evals_sorted = inv_evals(id);
    c2_evals_sorted = c2_evals(id);
    
    % Analyze for degeneracies
    energy = energies(1);
    degeneracy = 1;
    index = 1;
    trace_sigma2 = sigma2_evals_sorted(1);
    trace_c3 = c3_evals_sorted(1);
    trace_inv = inv_evals_sorted(1);
    trace_c2 = c2_evals_sorted(1);
    
    for i = 2:(2 * (N + 1)^2)
        if (abs(energies(i) - energy) < tolerance) % Next energy is degenerate
            degeneracy = degeneracy + 1;
            trace_sigma2 = trace_sigma2 + sigma2_evals_sorted(i);
            trace_c3 = trace_c3 + c3_evals_sorted(i);
            trace_inv = trace_inv + inv_evals_sorted(i);
            trace_c2 = trace_c2 + c2_evals_sorted(i);
        else % Next energy is new
            % Record stats for previous energy level
            energy_levels(index, 1) = energy;
            energy_levels(index, 2) = degeneracy;
            energy_levels(index, 3) = trace_sigma2;
            energy_levels(index, 4) = trace_c3;
            energy_levels(index, 5) = trace_inv;
            energy_levels(index, 6) = trace_c2;
            energy_levels(index, 7) = i - 1;
    
            % Reset variables for new energy level
            energy = energies(i);
            degeneracy = 1;
            index = index + 1;
            trace_sigma2 = sigma2_evals_sorted(i);
            trace_c3 = c3_evals_sorted(i);
            trace_inv = inv_evals_sorted(i);
            trace_c2 = c2_evals_sorted(i);
        end
    
        % Record energy level if we reach the end
        if (i == 2 * (N + 1)^2)
            energy_levels(index, 1) = energy;
            energy_levels(index, 2) = degeneracy;
            energy_levels(index, 3) = trace_sigma2;
            energy_levels(index, 4) = trace_c3;
            energy_levels(index, 5) = trace_inv;
            energy_levels(index, 6) = trace_c2;
            energy_levels(index, 7) = i;
        end
    end
    
    energy_levels = energy_levels(1:index, :);
    
    %% Separate the representation classes
    % Allocate space for each irrep
    elevels_a1g = zeros(2 * (N + 1)^2, 7);
    elevels_a2g = zeros(2 * (N + 1)^2, 7);
    elevels_eg = zeros(2 * (N + 1)^2, 7);
    elevels_t1g = zeros(2 * (N + 1)^2, 7);
    elevels_t2g = zeros(2 * (N + 1)^2, 7);
    elevels_a1u = zeros(2 * (N + 1)^2, 7);
    elevels_a2u = zeros(2 * (N + 1)^2, 7);
    elevels_eu = zeros(2 * (N + 1)^2, 7);
    elevels_t1u = zeros(2 * (N + 1)^2, 7);
    elevels_t2u = zeros(2 * (N + 1)^2, 7);
    
    index_a1g = 1;
    index_a2g = 1;
    index_eg = 1;
    index_t1g = 1;
    index_t2g = 1;
    index_a1u = 1;
    index_a2u = 1;
    index_eu = 1;
    index_t1u = 1;
    index_t2u = 1;
    
    % Round the energy_levels characters
    energy_levels_rounded = zeros(size(energy_levels, 1), 7);
    energy_levels_rounded(:, 2:6) = round(energy_levels(:, 2:6));
    energy_levels_rounded(:, 1) = energy_levels(:, 1);
    energy_levels_rounded(:, 7) = energy_levels(:, 7);
    
    % Fill in elevel info in each irrep
    for i = 1:size(energy_levels_rounded, 1)
        trace_e = energy_levels_rounded(i, 2);
        trace_sigma_d = energy_levels_rounded(i, 3);
        trace_c3 = energy_levels_rounded(i, 4);
        trace_inv = energy_levels_rounded(i, 5);
        trace_c2 = energy_levels_rounded(i, 6);
    
        traces = [trace_e, trace_sigma_d, trace_c3, trace_inv, trace_c2];
    
        if (isequal(traces, [1, 1, 1, 1, 1])) % A1g
            elevels_a1g(index_a1g, :) = energy_levels_rounded(i, :);
            index_a1g = index_a1g + 1;
        elseif (isequal(traces, [1, -1, 1, 1, -1])) % A2g
            elevels_a2g(index_a2g, :) = energy_levels_rounded(i, :);
            index_a2g = index_a2g + 1;
        elseif (isequal(traces, [2, 0, -1, 2, 0])) % Eg
            elevels_eg(index_eg, :) = energy_levels_rounded(i, :);
            index_eg = index_eg + 1;
        elseif (isequal(traces, [3, -1, 0, 3, -1])) % T1g
            elevels_t1g(index_t1g, :) = energy_levels_rounded(i, :);
            index_t1g = index_t1g + 1;
        elseif (isequal(traces, [3, 1, 0, 3, 1])) %T2g
            elevels_t2g(index_t2g, :) = energy_levels_rounded(i, :);
            index_t2g = index_t2g + 1;
        elseif (isequal(traces, [1, -1, 1, -1, 1])) %A1u
            elevels_a1u(index_a1u, :) = energy_levels_rounded(i, :);
            index_a1u = index_a1u + 1;
        elseif (isequal(traces, [1, 1, 1, -1, -1])) %A2u
            elevels_a2u(index_a2u, :) = energy_levels_rounded(i, :);
            index_a2u = index_a2u + 1;
        elseif (isequal(traces, [2, 0, -1, -2, 0])) %Eu
            elevels_eu(index_eu, :) = energy_levels_rounded(i, :);
            index_eu = index_eu + 1;
        elseif (isequal(traces, [3, 1, 0, -3, -1])) %T1u
            elevels_t1u(index_t1u, :) = energy_levels_rounded(i, :);
            index_t1u = index_t1u + 1;
        elseif (isequal(traces, [3, -1, 0, -3, +1])) %T2u
            elevels_t2u(index_t2u, :) = energy_levels_rounded(i, :);
            index_t2u = index_t2u + 1;
        elseif (isequal(traces, [2, 0, 2, 0, 2])) % Accidental degeneracy A1g + A1u
            elevels_a1g(index_a1g, :) = energy_levels_rounded(i, :);
            index_a1g = index_a1g + 1;
            elevels_a1u(index_a1u, :) = energy_levels_rounded(i, :); 
            index_a1u = index_a1u + 1;
        elseif (isequal(traces, [2, 0, 2, 0, -2])) % Accidental degeneracy A2g + A2u
            elevels_a2g(index_a2g, :) = energy_levels_rounded(i, :);
            index_a2g = index_a2g + 1;
            elevels_a2u(index_a2u, :) = energy_levels_rounded(i, :);
            index_a2u = index_a2u + 1;
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
    
    if (index_eg > 1)
        elevels_eg = elevels_eg(1:(index_eg - 1), :);
    end
    
    if (index_t1g > 1)
        elevels_t1g = elevels_t1g(1:(index_t1g - 1), :);
    end
    
    if (index_t2g > 1)
        elevels_t2g = elevels_t2g(1:(index_t2g - 1), :);
    end
    
    if (index_a1u > 1)
        elevels_a1u = elevels_a1u(1:(index_a1u - 1), :);
    end 
    
    if (index_a2u > 1)
        elevels_a2u = elevels_a2u(1:(index_a2u - 1), :);
    end
    
    if (index_eu > 1)
        elevels_eu = elevels_eu(1:(index_eu - 1), :);
    end
    
    if (index_t1u > 1)
        elevels_t1u = elevels_t1u(1:(index_t1u - 1), :);
    end
    
    if (index_t2u > 1)
        elevels_t2u = elevels_t2u(1:(index_t2u - 1), :);
    end
    
    %% Find energy level spacings
    % Make vector of energy level spacings
    spacings_a1g = zeros(size(elevels_a1g, 1) - 1, 1);
    spacings_a2g = zeros(size(elevels_a2g, 1) - 1, 1);
    spacings_eg = zeros(size(elevels_eg, 1) - 1, 1);
    spacings_t1g = zeros(size(elevels_t1g, 1) - 1, 1);
    spacings_t2g = zeros(size(elevels_t2g, 1) - 1, 1);
    spacings_a1u = zeros(size(elevels_a1u, 1) - 1, 1);
    spacings_a2u = zeros(size(elevels_a2u, 1) - 1, 1);
    spacings_eu = zeros(size(elevels_eu, 1) - 1, 1);
    spacings_t1u = zeros(size(elevels_t1u, 1) - 1, 1);
    spacings_t2u = zeros(size(elevels_t2u, 1) - 1, 1);
    
    for i = 1:(size(elevels_a1g, 1) - 1)
        spacings_a1g(i) = abs(elevels_a1g(i, 1) - elevels_a1g(i+1, 1));
    end
    
    for i = 1:(size(elevels_a2g, 1) - 1)
        spacings_a2g(i) = abs(elevels_a2g(i, 1) - elevels_a2g(i+1, 1));
    end
    
    for i = 1:(size(elevels_eg, 1) - 1)
        spacings_eg(i) = abs(elevels_eg(i, 1) - elevels_eg(i+1, 1));
    end
    
    for i = 1:(size(elevels_t1g, 1) - 1)
        spacings_t1g(i) = abs(elevels_t1g(i, 1) - elevels_t1g(i+1, 1));
    end
    
    for i = 1:(size(elevels_t2g, 1) - 1)
        spacings_t2g(i) = abs(elevels_t2g(i, 1) - elevels_t2g(i+1, 1));
    end
    
    for i = 1:(size(elevels_a1u, 1) - 1)
        spacings_a1u(i) = abs(elevels_a1u(i, 1) - elevels_a1u(i+1, 1));
    end
    
    for i = 1:(size(elevels_a2u, 1) - 1)
        spacings_a2u(i) = abs(elevels_a2u(i, 1) - elevels_a2u(i+1, 1));
    end
    
    for i = 1:(size(elevels_eu, 1) - 1)
        spacings_eu(i) = abs(elevels_eu(i, 1) - elevels_eu(i+1, 1));
    end
    
    for i = 1:(size(elevels_t1u, 1) - 1)
        spacings_t1u(i) = abs(elevels_t1u(i, 1) - elevels_t1u(i+1, 1));
    end
    
    for i = 1:(size(elevels_t2u, 1) - 1)
        spacings_t2u(i) = abs(elevels_t2u(i, 1) - elevels_t2u(i+1, 1));
    end
    
    % Compute level spacing ratios
    r_a1g = LSR(spacings_a1g);
    r_a1u = LSR(spacings_a1u);
    r_a2g = LSR(spacings_a2g);
    r_a2u = LSR(spacings_a2u);
    r_eg = LSR(spacings_eg);
    r_eu = LSR(spacings_eu);
    r_t1g = LSR(spacings_t1g);
    r_t1u = LSR(spacings_t1u);
    r_t2g = LSR(spacings_t2g);
    r_t2u = LSR(spacings_t2u);

    % Fill in r array values
    r_array(index_r, :) = [N, r_a1g, r_a2g, r_eg, r_t1g, r_t2g, r_a1u, r_a2u, r_eu, r_t1u, r_t2u];
   
    % Fill in size of irreps
    size_array(index_size, :) = [N, size(elevels_a1g, 1), size(elevels_a2g, 1), size(elevels_eg, 1), ...
        size(elevels_t1g, 1), size(elevels_t2g, 1), size(elevels_a1u, 1), size(elevels_a2u, 1), ...
        size(elevels_eu, 1), size(elevels_t1u, 1), size(elevels_t2u, 1)];

    % Find proportion of solvable states
    solvable_prop(index_size, :) = [size(elevels_a1g, 1), size(elevels_a2g, 1), size(elevels_a1u, 1), size(elevels_a2u, 1)] / (2 * (N + 1)^2);

    %% Plot the level spacings histogram and show r values
    % Plot histograms
    figure(index_n)
    tiledlayout(2, 5,'TileSpacing', 'tight','Padding','Tight')
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
    histogram(spacings_eg, ceil(size(elevels_eg, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('Eg')
    subtitle(['r = ' num2str(r_eg) '; total = ' num2str(size_array(index_size, 4))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_t1g, ceil(size(elevels_t1g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('T1g')
    subtitle(['r = ' num2str(r_t1g) '; total = ' num2str(size_array(index_size, 5))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_t2g, ceil(size(elevels_t2g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('T2g')
    subtitle(['r = ' num2str(r_t2g) '; total = ' num2str(size_array(index_size, 6))])
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
    histogram(spacings_eu, ceil(size(elevels_eu, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('Eu')
    subtitle(['r = ' num2str(r_eu) '; total = ' num2str(size_array(index_size, 9))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_t1u, ceil(size(elevels_t1u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('T1u')
    subtitle(['r = ' num2str(r_t1u) '; total = ' num2str(size_array(index_size, 10))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_t2u, ceil(size(elevels_t2u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('T2u')
    subtitle(['r = ' num2str(r_t2u) '; total = ' num2str(size_array(index_size, 11))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    set(figure(index_n),'position',[0,100,3000,400])
    saveas(gcf, [folderName '/N=' num2str(N) '_elevel_hist.jpeg']);
    
    % Increment indices
    index_n = index_n + 1;
    index_size = index_size + 1;
    index_r = index_r + 1;
end

%% Plot r as a function of size in each irrep
figure(index_n)
tiledlayout(2, 5, 'TileSpacing', 'tight', 'Padding', 'Tight')
ylow = 0.3;
yhigh = 0.7;

nexttile
plot(N_nums.', r_array(:, 2).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('A1g')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_array(:, 3).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('A2g')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_array(:, 4).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Eg')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_array(:, 5).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('T1g')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_array(:, 6).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('T2g')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_array(:, 7).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('A1u')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_array(:, 8).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('A2u')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_array(:, 9).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Eu')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_array(:, 10).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('T1u')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_array(:, 11).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('T2u')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow, yhigh])

set(figure(index_n),'position',[0,100,3000,400])
saveas(gcf, [folderName '/LSR_plot.jpeg']);

%% Plot proportions of solvable states
figure(index_n + 1)
plot(N_nums, sum(solvable_prop, 2), 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Proportion of states which are solvable')
xlabel('N')
yline(1/12, '--', '1/12')
ylim([0.0, 0.1])
ylabel('Proportion')
saveas(gcf, [folderName '/solv_prop.jpeg']);

%% Save variables
save([folderName '/r_array.mat'],'r_array','-v7.3')
save([folderName '/size_array.mat'],'size_array','-v7.3')
save([folderName '/elevels_a1g.mat'],'elevels_a1g','-v7.3')
save([folderName '/elevels_a1u.mat'],'elevels_a1u','-v7.3')
save([folderName '/elevels_a2g.mat'],'elevels_a2g','-v7.3')
save([folderName '/elevels_a2u.mat'],'elevels_a2u','-v7.3')
save([folderName '/elevels_t1g.mat'],'elevels_t1g','-v7.3')
save([folderName '/elevels_t1u.mat'],'elevels_t1u','-v7.3')
save([folderName '/elevels_t2g.mat'],'elevels_t2g','-v7.3')
save([folderName '/elevels_t2u.mat'],'elevels_t2u','-v7.3')
save([folderName '/elevels_eg.mat'],'elevels_eg','-v7.3')
save([folderName '/elevels_eu.mat'],'elevels_eu','-v7.3')

%% Functions

% Gives width of row within face CHECKED
function w = row_width(index, N) 
    face = face_from_index(index, N);
    y = yfacecoord(index, N);
    if ((face == 1) || (face == 2) || (face == 4) || (face == 6))
        w = 1 + 2 * (y - 1);
    else
        w = N - 2 * (y - 1);
    end
end

% Gives y coordinate of index within each face CHECKED
function y = yfacecoord(index, N) 
    sites_per_face = ((N + 1) / 2)^2;
    if (index <= sites_per_face)
        y = ceil(sqrt(index));
    elseif (index > 7 * sites_per_face)
        y = (N + 1) / 2 + 1 - ceil(sqrt(2 * (N + 1)^2 + 1 - index));
    else
        % remove bottom triangle
        index = index - sites_per_face;
        row_length = 3 * (N + 1);

        y = floor((index - 1) / row_length) + 1;
    end
end

% Gives x coordinate of index within each face CHECKED
function x = xfacecoord(index, N) 
    face = face_from_index(index, N);
    y_coord = yfacecoord(index, N);
    sites_per_face = ((N + 1) / 2)^2;
    total_sites = 2 * (N + 1)^2;
    row_w = row_width(index, N);

    if (face == 1)
        x = index - (y_coord - 1)^2; 
    elseif (face == 8)
        x = row_w - (total_sites + 1 - index - ((N + 1) / 2 - y_coord)^2) + 1;
    elseif (face == 2) 
        x = index - sites_per_face - (y_coord - 1) * 3 * (N + 1);
    elseif (face == 4)
        x = index - sites_per_face - (y_coord - 1) * 3 * (N + 1) - (N + 1);
    elseif (face == 6)
        x = index - sites_per_face - (y_coord - 1) * 3 * (N + 1) - 2 * (N + 1);
    elseif (face == 3)
        x = index - sites_per_face - (y_coord - 1) * 3 * (N + 1) - (N + 1 - row_w);
    elseif (face == 5)
        x = index - sites_per_face - (y_coord - 1) * 3 * (N + 1) - (N + 1) - (N + 1 - row_w);
    elseif (face == 7)
        x = index - sites_per_face - (y_coord - 1) * 3 * (N + 1) - 2 * (N + 1) - (N + 1 - row_w);
    end
end

% Gives face containing given index CHECKED
function f = face_from_index(index, N) 
    sites_per_face = ((N + 1) / 2)^2;
    if (index <= sites_per_face)
        f = 1;
    elseif (index > 7 * sites_per_face)
        f = 8;
    else
        % remove bottom triangle
        index = index - sites_per_face;
        row_length = 3 * (N + 1);

        face_pair = mod(ceil(index / (N + 1)) - 1, 3) + 1;
        row_num = ceil(index / row_length);

        short_width = 1 + 2 * (row_num - 1);
        long_width = (N + 1) - short_width;

        index = mod(index - 1, short_width + long_width) + 1;
        if (index <= short_width)
            f = 2 * face_pair; 
        else
            f = 2 * face_pair + 1;
        end
    end
end

% From face coordinate and face number give index CHECKED
function i = index_from_coord(x, y, face, N)
    sites_per_face = ((N + 1) / 2)^2;
    total_sites = 2 * (N + 1)^2;

    if (face == 1)
        i = (y - 1)^2 + x;
    elseif (face == 2)
        i = sites_per_face + (y - 1) * (3 * (N + 1)) + x;
    elseif (face == 3)
        i = sites_per_face + (y - 1) * (3 * (N + 1)) + (1 + 2 * (y - 1)) + x;
    elseif (face == 4)
        i = sites_per_face + (y - 1) * (3 * (N + 1)) + (N + 1) + x;
    elseif (face == 5)
        i = sites_per_face + (y - 1) * (3 * (N + 1)) + (N + 1) + (1 + 2 * (y - 1)) + x;
    elseif (face == 6)
        i = sites_per_face + (y - 1) * (3 * (N + 1)) + 2 * (N + 1) + x;
    elseif (face == 7)
        i = sites_per_face + (y - 1) * (3 * (N + 1)) + 2 * (N + 1) + (1 + 2 * (y - 1)) + x;
    elseif (face == 8)
        i = total_sites - ((N + 1) / 2 + 1 - y)^2 + x;
    end
end

% Given index on discretized octahedron returns new index after sigma_d symmetry
% CHECKED
function i = sigma_d_index(index, N)
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    if (face == 1) 
        new_face = 3;
        new_x = 2 * (y_coord - 1) + 1 - x_coord + 1;
        new_y = ceil(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 2) 
        new_face = 7;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) + mod(x_coord - 1, 2);
        new_y = y_coord - floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 3) 
        new_face = 1;
        new_x = 2 * (y_coord - 1) + 1 + mod(x_coord - 1, 2);
        new_y = y_coord + floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 4) 
        new_face = 5;
        new_x = 2 * (y_coord - 1) + 1 - (x_coord - 1);
        new_y = ceil(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 5) 
        new_face = 4;
        new_x = 2 * (y_coord - 1) + 1 + mod(x_coord - 1, 2);
        new_y = y_coord + floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 6) 
        new_face = 8;
        new_x = 2 * (y_coord - 1) + 1 - (x_coord - 1);
        new_y = ceil(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 7) 
        new_face = 2;
        new_x = N + 1 - (2 * (y_coord - 1) + 1) - x_coord + 1;
        new_y = (N + 1) / 2 - ceil(x_coord / 2) + 1;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 8) 
        new_face = 6;
        new_x = 2 * (y_coord - 1) + 1 + mod(x_coord - 1, 2);
        new_y = y_coord + floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    end
end

% Given index on discretized octahedron returns new index after C3 symmetry
% CHECKED
function i = c3_index(index, N)
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    if (face == 1) 
        new_face = 6;
        new_x = x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 2) 
        new_face = 1;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) + mod(x_coord - 1, 2);
        new_y = (N + 1) / 2 + 1 - y_coord + floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 3) 
        new_face = 5;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) - x_coord + 1;
        new_y = ceil(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 4) 
        new_face = 4;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) + mod(x_coord - 1, 2);
        new_y = (N + 1) / 2 + 1 - y_coord + floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 5) 
        new_face = 8;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) - x_coord + 1;
        new_y = ceil(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 6) 
        new_face = 2;
        new_x = 1 + 2 * (y_coord - 1) - x_coord + 1;
        new_y = (N + 1) / 2 + 1 - ceil(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 7) 
        new_face = 7;
        new_x = 2 * (y_coord - 1) + 1 + mod(x_coord - 1, 2);
        new_y = (N + 1) / 2 + 1 - y_coord - floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 8) 
        new_face = 3;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) - x_coord + 1;
        new_y = ceil(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    end
end

% Given index on discretized octahedron returns new index after C2 symmetry
% CHECKED
function i = c2_index(index, N)
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    face = face_from_index(index, N);
    row_w = row_width(index, N);
    if (face == 1) % in face 1
        new_face = 8;
        new_x = row_w + 1 - x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 2) % in face 2 
        new_face = 7;
        new_x = row_w + 1 - x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        %fprintf(['right_col = ' num2str(1 + 2 * (new_y - 1)) ' index = ' num2str(index) ' new_x = ' num2str(new_x) ' new_y = ' num2str(new_y) '\n'])
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 3)
        new_face = 6;
        new_x = row_w + 1 - x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 4) % in face 4
        new_face = 5;
        row_w = row_width(index, N);
        new_x = row_w + 1 - x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        %fprintf(['left_col = ' num2str(left_col(x_coord, y_coord, N)) ' index = ' num2str(index) ' new_x = ' num2str(new_x) ' new_y = ' num2str(new_y) '\n'])
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 5) % in face 5
        new_face = 4;
        new_x = row_w + 1 - x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 6) % in face 6 
        new_face = 3;
        new_x = row_w + 1 - x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        %fprintf(['right_col = ' num2str(1 + 2 * (new_y - 1)) ' index = ' num2str(index) ' new_x = ' num2str(new_x) ' new_y = ' num2str(new_y) '\n'])
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 7)
        new_face = 2;
        new_x = row_w + 1 - x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 8) % in face 8
        new_face = 1;
        new_x = row_w + 1 - x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        %fprintf(['left_col = ' num2str(left_col(x_coord, y_coord, N)) ' index = ' num2str(index) ' new_x = ' num2str(new_x) ' new_y = ' num2str(new_y) '\n'])
        i = index_from_coord(new_x, new_y, new_face, N);
    end
end

% Given index on discretized octahedron returns new index after inv symmetry
% CHECKED
function i = inv_index(index, N)
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    if (face == 1) 
        new_face = 7;
        new_x = 2 * (y_coord - 1) + 1 - x_coord + 1;
        new_y = ceil(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 2) 
        new_face = 5;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) + mod(x_coord - 1, 2);
        new_y = y_coord - floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 3) 
        new_face = 1;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) - x_coord + 1;
        new_y = (N + 1) / 2 - ceil(x_coord / 2) + 1;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 4) 
        new_face = 3;
        new_x = 2 * (y_coord - 1) + 1 - x_coord + 1;
        new_y = ceil(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 5) 
        new_face = 2;
        new_x = 2 * (y_coord - 1) + 1 + mod(x_coord - 1, 2);
        new_y = y_coord + floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 6) 
        new_face = 8;
        new_x = x_coord;
        new_y = (N + 1) / 2 - y_coord + 1;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 7) 
        new_face = 6;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) - x_coord + 1;
        new_y = (N + 1) / 2 - ceil(x_coord / 2) + 1;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 8) 
        new_face = 4;
        new_x = 2 * (y_coord - 1) + 1 + mod(x_coord - 1, 2);
        new_y = y_coord + floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    end
end

% Functions to find indices of H matrix
function y = upper(index, N) % CHECKED
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    row_w = row_width(index, N);

    if (face == 1 && y_coord < (N + 1) / 2) % in face 1 but not top row
        y = index_from_coord(x_coord + 1, y_coord + 1, 1, N);
    elseif (face == 1 && y_coord == (N + 1) / 2) % in face 1 top row
        y = index_from_coord(x_coord, 1, 5, N);
    elseif (face == 8)
        y = index_from_coord(x_coord - 1, y_coord + 1, 8, N);
    elseif (face == 2 && y_coord == (N + 1) / 2) % face 2 top row
        new_y = (N + 1) / 2 - (x_coord - 1) / 2;
        y = index_from_coord(1, new_y, 8, N);
    elseif (face == 6 && y_coord == (N + 1) / 2) % face 6 top row
        new_y = (x_coord - 1) / 2 + 1;
        y = index_from_coord(row_w + 1 - x_coord, new_y, 8, N);
    elseif (face == 4 && y_coord == (N + 1) / 2) % face 4 top row
        y = index + 2 * (N + 1);
    elseif (mod(face, 2) == 0) % face = 2, 4, or 6
        y = index_from_coord(x_coord + 1, y_coord + 1, face, N);
    else % face = 3, 5, or 7
        y = index_from_coord(x_coord - 1, y_coord + 1, face, N);
    end
end

function y = lower(index, N) % CHECKED
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);

    if (face == 1)
        y = index_from_coord(x_coord - 1, y_coord - 1, 1, N);
    elseif (face == 3 && y_coord == 1) % face 3 bottom row
        new_y = (x_coord - 1) / 2 + 1;
        y = index_from_coord(1, new_y, 1, N);
    elseif (face == 5 && y_coord == 1) % face 5 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 1, N);
    elseif (face == 7 && y_coord == 1) % face 7 bottom row
        new_y = (N + 1) / 2 + 1 - ((x_coord - 1) / 2 + 1); 
        y = index_from_coord(N + 1 - x_coord, new_y, 1, N);
    elseif (face == 8 && y_coord == 1) % face 8 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 4, N);
    elseif (face == 8)
        y = index_from_coord(x_coord, y_coord - 1, 8, N) + 1;
    else
        y = index - 3 * (N + 1) - 1;
    end
end
function y = right(index, N) % CHECKED
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    row_w = row_width(index, N);

    if (face == 1 && x_coord == row_w) % face 1 right side
        new_x = N + 1 - (1 + 2 * (y_coord - 1));
        y = index_from_coord(new_x, 1, 7, N);
    elseif (face == 8 && x_coord == row_w) % face 8 right side
        new_x = 1 + 2 * (y_coord - 1);
        y = index_from_coord(new_x, (N + 1) / 2, 6, N);
    elseif (face == 7 && x_coord == row_w) % face 7 right side
        y = index_from_coord(1, y_coord, 2, N);
    else
        y = index + 1; 
    end
end

function y = left(index, N) % CHECKED
    face = face_from_index(index, N);
    row_w = row_width(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);

    if (face == 1 && x_coord == 1) % face 1 left side
        new_x = 1 + 2 * (y_coord - 1);
        y = index_from_coord(new_x, 1, 3, N);
    elseif (face == 2 && x_coord == 1) % face 2 left side
        y = index_from_coord(N + 1 - row_w, y_coord, 7, N);
    elseif (face == 8 && x_coord == 1)
        new_x = N + 1 - (1 + 2 * (y_coord - 1));
        y = index_from_coord(new_x, (N + 1) / 2, 2, N);
    else
        y = index - 1;
    end  
end

% Indicate whether nearest vertical neighbor is above or below CHECKED
function answer = has_upper(index, N)
    y = yfacecoord(index, N);
    face = face_from_index(index, N);

    if (face == 1 || face == 8)
        if (mod(y, 2) - mod(index, 2) == 0)
            answer = true;
        else 
            answer = false;
        end
    elseif (mod(index, 2) == 0)
        answer = true;
    else
        answer = false;
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