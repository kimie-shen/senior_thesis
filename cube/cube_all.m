% Calculate r values for various N
tolerance = 1e-8;  % Tolerance for degenerate energy levels
N_start = 5;
N_spacing = 5;
N_end = 60;

r_array = zeros(ceil((N_end - N_start)/N_spacing), 11);
size_array = zeros(ceil((N_end - N_start)/N_spacing), 11);

index_n = 1;
index_r = 1;
index_size = 1;

for N = N_start:N_spacing:N_end
    %% Diagonalize H and sigma_d matrices
    tic;
    fprintf(['\n    N = ' num2str(N) '\n'])
    % Generate Hamiltonian matrix
    L = 1/N; % Spacing
    H = zeros(6 * N^2, 6 * N^2); % Hamiltonian matrix
    
    for i = 1:6*N^2
        H(i, i) = 4/L^2;
        H(i, upper(i, N)) = -1/L^2;
        H(i, lower(i, N)) = -1/L^2;
        H(i, right(i, N)) = -1/L^2;
        H(i, left(i, N)) = -1/L^2;
    end
    
    % Generate matrix for sigma_d symmetry 
    sigma_d = zeros(6 * N^2, 6 * N^2);
    for i = 1: (6 * N^2)
        sigma_d(i, sigma_d_index(i, N)) = 1;
    end
    
    % Calculate eigenvectors and eigenvalues
    [eigenvectors, eigenvalues] = eig(H + ((1 + sqrt(5))/2)*sigma_d);
    fprintf('evals and evecs done \n')
    toc
    
    %% Generate Symmetry Matrices
    
    % Generate matrix for C3
    C3 = zeros(6 * N^2, 6 * N^2);
    
    for i = 1: (6 * N^2)
        C3(i, c3_index(i, N)) = 1;
    end
    
    % Generate matrix for i
    inv = zeros(6 * N^2, 6 * N^2);
    for i = 1: (6 * N^2)
        inv(i, inv_index(i, N)) = 1;
    end
    
    % Generate matrix for C2
    C2 = zeros(6 * N^2, 6 * N^2);
    for i = 1: (6 * N^2)
        C2(i, c2_index(i, N)) = 1;
    end
    
    %% Calculate characters of symmetries
    tic;
    % Preallocate data array
    sigma2_evals = zeros(6 * N^2, 1);
    c3_evals = zeros(6 * N^2, 1);
    inv_evals = zeros(6 * N^2, 1);
    c2_evals = zeros(6 * N^2, 1);
    
    % Subtract off sigma2 eigenvalues from eigenvalues 
    for i = 1:(6 * N^2)
        inv_evals(i) = (eigenvectors(:, i).' * inv * eigenvectors(:, i));
    end
    fprintf('inv done \n')
    
    for i = 1:(6 * N^2)
        sigma2_evals(i) = (eigenvectors(:, i).' * sigma_d * eigenvectors(:, i));
        eigenvalues(i, i) = eigenvalues(i, i) - ((1 + sqrt(5))/2) * sigma2_evals(i);
    end
    fprintf('sigma_d done \n')
    
    for i = 1:(6 * N^2)
        c3_evals(i) = (eigenvectors(:, i).' * C3 * eigenvectors(:, i));
    end
    fprintf('c3 done \n')
    
    for i = 1:(6 * N^2)
        c2_evals(i) = (eigenvectors(:, i).' * C2 * eigenvectors(:, i));
    end
    fprintf('c2 done \n')
    toc;
    
    %% Analyze Degeneracies
    % Reorder energies and evecs
    energy_levels = zeros(6 * N^2, 7);
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
    
    for i = 2:(6 * N^2)
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
        if (i == 6 * N^2)
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
    elevels_a1g = zeros(6 * N^2, 7);
    elevels_a2g = zeros(6 * N^2, 7);
    elevels_eg = zeros(6 * N^2, 7);
    elevels_t1g = zeros(6 * N^2, 7);
    elevels_t2g = zeros(6 * N^2, 7);
    elevels_a1u = zeros(6 * N^2, 7);
    elevels_a2u = zeros(6 * N^2, 7);
    elevels_eu = zeros(6 * N^2, 7);
    elevels_t1u = zeros(6 * N^2, 7);
    elevels_t2u = zeros(6 * N^2, 7);
    
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
    
    %% Plot the level spacings histogram and show r values
    % Plot histograms
    figure(index_n)
    tiledlayout(2,5)
    bin_factor = 5;
    
    nexttile
    histogram(spacings_a1g, ceil(size(elevels_a1g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title(['A1g: r = ' num2str(r_a1g)])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_a2g, ceil(size(elevels_a2g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title(['A2g: r = ' num2str(r_a2g)])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_eg, ceil(size(elevels_eg, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title(['Eg: r = ' num2str(r_eg)])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_t1g, ceil(size(elevels_t1g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title(['T1g: r = ' num2str(r_t1g)])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_t2g, ceil(size(elevels_t2g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title(['T2g: r = ' num2str(r_t2g)])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_a1u, ceil(size(elevels_a1u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title(['A1u: r = ' num2str(r_a1u)])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_a2u, ceil(size(elevels_a2u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title(['A2u: r = ' num2str(r_a2u)])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_eu, ceil(size(elevels_eu, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title(['Eu: r = ' num2str(r_eu)])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_t1u, ceil(size(elevels_t1u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title(['T1u: r = ' num2str(r_t1u)])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_t2u, ceil(size(elevels_t2u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title(['T2u: r = ' num2str(r_t2u)])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    set(figure(index_n),'position',[0,100,1500,400])
    index_n = index_n + 1;

    % Fill in r array values
    r_array(index_r, :) = [N, r_a1g, r_a2g, r_eg, r_t1g, r_t2g, r_a1u, r_a2u, r_eu, r_t1u, r_t2u];
    index_r = index_r + 1;

    % Fill in size of irreps
    size_array(index_size, :) = [N, size(elevels_a1g, 1), size(elevels_a2g, 1), size(elevels_eg, 1), ...
        size(elevels_t1g, 1), size(elevels_t2g, 1), size(elevels_a1u, 1), size(elevels_a2u, 1), ...
        size(elevels_eu, 1), size(elevels_t1u, 1), size(elevels_t2u, 1)];
    index_size = index_size + 1;
end

%% Plot r as a function of size in each irrep
figure(index_n)
tiledlayout(2,5)

nexttile
plot(N_start:N_spacing:N_end, r_array(:, 2), '-o')
title('A1g')
xlabel('N')
ylabel('r')
hold on 
yline(0.39)
ylim([0.3, 0.8])

nexttile
plot(N_start:N_spacing:N_end, r_array(:, 3), '-o')
title('A2g')
xlabel('N')
ylabel('r')
hold on 
yline(0.39)
ylim([0.3, 0.8])

nexttile
plot(N_start:N_spacing:N_end, r_array(:, 4), '-o')
title('Eg')
xlabel('N')
ylabel('r')
hold on 
yline(0.53)
ylim([0.3, 0.8])

nexttile
plot(N_start:N_spacing:N_end, r_array(:, 5), '-o')
title('T1g')
xlabel('N')
ylabel('r')
hold on 
yline(0.53)
ylim([0.3, 0.8])

nexttile
plot(N_start:N_spacing:N_end, r_array(:, 6), '-o')
title('T2g')
xlabel('N')
ylabel('r')
hold on 
yline(0.53)
ylim([0.3, 0.8])

nexttile
plot(N_start:N_spacing:N_end, r_array(:, 7), '-o')
title('A1u')
xlabel('N')
ylabel('r')
hold on 
yline(0.39)
ylim([0.3, 0.8])

nexttile
plot(N_start:N_spacing:N_end, r_array(:, 8), '-o')
title('A2u')
xlabel('N')
ylabel('r')
hold on 
yline(0.39)
ylim([0.3, 0.8])

nexttile
plot(N_start:N_spacing:N_end, r_array(:, 9), '-o')
title('Eu')
xlabel('N')
ylabel('r')
hold on 
yline(0.53)
ylim([0.3, 0.8])

nexttile
plot(N_start:N_spacing:N_end, r_array(:, 10), '-o')
title('T1u')
xlabel('N')
ylabel('r')
hold on 
yline(0.53)
ylim([0.3, 0.8])

nexttile
plot(N_start:N_spacing:N_end, r_array(:, 11), '-o')
title('T2u')
xlabel('N')
ylabel('r')
hold on 
yline(0.53)
ylim([0.3, 0.8])

set(figure(index_n),'position',[0,100,1500,400])

%% Plot proportions of elevels in each state
n_runs = (N_end - N_start) / N_spacing + 1;
solvable_sizes = zeros(n_runs, 4);

for i = 1:n_runs
    solvable_sizes(i, 1:2) = size_array(i, 2:3);
    solvable_sizes(i, 3:4) = size_array(i, 6:7);
end

figure(index_n + 1)
bar(solvable_sizes, 'stacked')
title('Proportions')
xlabel('N')
ylabel('Proportion')



%% Functions

% Gives x coordinate of index within each face
function x = xfacecoord(index, N)
    x = mod(index - 1, N) + 1;
end

% Gives y coordinate of index within each face
function y = yfacecoord(index, N) 
    if (index <= N^2) % in face 1
        y = ceil(index / N);
    elseif (index <= 2 * N^2) % in face 2 
        y = ceil((index - N^2) / N);
    elseif (index <= 5 * N^2) % in faces 3, 4, or 5
        y = ceil((index - 2 * N^2) / (3 * N));
    else % in face 6
        y = ceil((index - 5 * N^2) / N);
    end
end

% From face coordinate and face number give index
function i = index_from_coord(x, y, face, N)
    if (face == 1)
        i = x + N * (y - 1);
    elseif (face == 2)
        i = N^2 + x + N * (y - 1);
    elseif (face == 3)
        i = 2 * N^2 + x + 3 * N * (y - 1);
    elseif (face == 4)
        i = 2 * N^2 + N + x + 3 * N * (y - 1);
    elseif (face == 5)
        i = 2 * N^2 + 2 * N + x + 3 * N * (y - 1);
    else 
        i = 5 * N^2 + x + N * (y - 1);
    end
end

% Given index on discretized cube returns new index after sigma_d symmetry
function i = sigma_d_index(index, N)
    if (index <= N^2) % in face 1
        new_face = 1;
        new_x = N + 1 - yfacecoord(index, N);
        new_y = N + 1 - xfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (index <= 2 * N^2) % in face 2 
        new_face = 3;
        new_x = yfacecoord(index, N);
        new_y = xfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (index > 5 * N^2) % in face 6
        new_face = 5;
        new_x = yfacecoord(index, N);
        new_y = xfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (mod(index - 2 * N^2 - 1, 3 * N) + 1 <= N) % in face 3
        new_face = 2;
        new_x = yfacecoord(index, N);
        new_y = xfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (mod(index - 2 * N^2 - 1, 3 * N) + 1 <= 2 * N) % in face 4
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

% Given index on discretized cube returns new index after C3 symmetry
function i = c3_index(index, N)
    if (index <= N^2) % in face 1
        new_face = 5;
        new_x = N + 1 - yfacecoord(index, N);
        new_y = xfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (index <= 2 * N^2) % in face 2 
        new_face = 4;
        new_x = N + 1 - yfacecoord(index, N);
        new_y = xfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (index > 5 * N^2) % in face 6
        new_face = 1;
        new_x = yfacecoord(index, N);
        new_y = N + 1 - xfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (mod(index - 2 * N^2 - 1, 3 * N) + 1 <= N) % in face 3
        new_face = 2;
        new_x = N + 1 - xfacecoord(index, N);
        new_y = N + 1 - yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (mod(index - 2 * N^2 - 1, 3 * N) + 1 <= 2 * N) % in face 4
        new_face = 3;
        new_x = N + 1 - yfacecoord(index, N);
        new_y = xfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    else  % in face 5
        new_face = 6;
        new_x = xfacecoord(index, N);
        new_y = yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    end
end

% Given index on discretized cube returns new index after inv symmetry
function i = inv_index(index, N)
    if (index <= N^2) % in face 1
        new_face = 4;
        new_x = N + 1 - xfacecoord(index, N);
        new_y = yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (index <= 2 * N^2) % in face 2 
        new_face = 6;
        new_x = N + 1 - xfacecoord(index, N);
        new_y = yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (index > 5 * N^2) % in face 6
        new_face = 2;
        new_x = N + 1 - xfacecoord(index, N);
        new_y = yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (mod(index - 2 * N^2 - 1, 3 * N) + 1 <= N) % in face 3
        new_face = 5;
        new_x = xfacecoord(index, N);
        new_y = N + 1 - yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (mod(index - 2 * N^2 - 1, 3 * N) + 1 <= 2 * N) % in face 4
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

% Given index on discretized cube returns new index after C2 symmetry
function i = c2_index(index, N)
    if (index <= N^2) % in face 1
        new_face = 6;
        new_x = N + 1 - xfacecoord(index, N);
        new_y = N + 1 - yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (index <= 2 * N^2) % in face 2 
        new_face = 4;
        new_x = N + 1 - xfacecoord(index, N);
        new_y = N + 1 - yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (index > 5 * N^2) % in face 6
        new_face = 1;
        new_x = N + 1 - xfacecoord(index, N);
        new_y = N + 1 - yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (mod(index - 2 * N^2 - 1, 3 * N) + 1 <= N) % in face 3
        new_face = 5;
        new_x = yfacecoord(index, N);
        new_y = N + 1 - xfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (mod(index - 2 * N^2 - 1, 3 * N) + 1 <= 2 * N) % in face 4
        new_face = 2;
        new_x = N + 1 - xfacecoord(index, N);
        new_y = N + 1 -yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    else  % in face 5
        new_face = 3;
        new_x = N + 1 - yfacecoord(index, N);
        new_y = xfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    end
end

% Functions to find indices of H matrix
function y = upper(x, N) 
    if (x > 2 * N^2 - N && x <= 2 * N^2) % x in top row of lower arm 
        y = 2 * N^2 + N + mod(x - 1, N) + 1;
    elseif (x > 2 * N^2 && x <= 5 * N^2 - 3 * N) % x in middle arm but not top row
        y = x + 3 * N;
    elseif (x > 5 * N^2 - 3 * N && x <= 5 * N^2 - 2 * N) % x top row of middle arms, left side
        dist_from_vertex = x - (5 * N^2 - 3 * N); % distance from leftmost vertex
        y = 6 * N^2 - dist_from_vertex * N + 1; 
    elseif (x > 5 * N^2 - 2 * N && x <= 5 * N^2 - 1 * N) % x in top row of middle arms, center side
        y = 5 * N^2 + mod(x - 1, N) + 1;
    elseif (x > 5 * N^2 - 1 * N && x <= 5 * N^2) % x in top row of middle arms, right side
        dist_from_fold = x - (5 * N^2 - 1 * N);
        y = 5 * N^2 + dist_from_fold * N;
    elseif (x > 5 * N^2 + N * (N - 1)) % x in top row of upper arm
        y = mod(x - 1, N) + 1;
    else % x in non-top row upper arm or non-top row lower arm
        y = x + N;
    end
end

function y = lower(x, N) 
    if (x > 0 && x <= N) % x in lower arm bottom row
        y = 6 * N^2 - N + mod(x - 1, N) + 1;
    elseif (x > 2 * N^2 && x <= 2 * N^2 + N) % x in middle arms, bottom row, left side
        dist_from_vertex = x - 2 * N^2;
        y = N^2 + (dist_from_vertex - 1) * N + 1;
    elseif (x > 2 * N^2 + N && x <= 2 * N^2 + 2 * N) % x in middle arms, bottom row, center
        y = x - 2 * N;
    elseif (x > 2 * N^2 + 2 * N && x <= 2 * N^2 + 3 * N) % x in middle arms, bottom row, right side
        dist_from_fold = x - (2 * N^2 + 2 * N);
        y = 2 * N^2 - (dist_from_fold - 1) * N ;
    elseif (x > 2 * N^2 + 3 * N && x <= 5 * N^2) % x in middle arms, above bottom row
        y = x - 3 * N;
    elseif (x > 5 * N^2 && x <= 5 * N^2 + N) % x in upper arm, bottom row
        y = x - 2 * N;
    else % x in upper arm, above bottom row or in lower arm above bottom row
        y = x - N;
    end
end

function y = right(x, N)
    if (x <= N^2 && floor(x / N) == x / N) % x in lower half of lower arm, rightmost column
        row_from_bottom = x/N;
        y = 5 * N^2 - (row_from_bottom - 1) * 3 * N;
    elseif (x > N^2 && x <= 2 * N^2 && floor(x / N) == x / N) % x in upper half of lower arm, rightmost column
        row_from_top = (2 * N^2 - x)/N + 1;
        y = 2 * N^2 + 2 * N + row_from_top;
    elseif (x > 2 * N^2 && x <= 5 * N^2 && floor((x - 2 * N^2)/(3 * N)) == (x - 2 * N^2)/(3 * N)) % x in middle arm, rightmost column 
        row_from_top_2 = (5 * N^2 - x)/ (3 * N) + 1;
        y = row_from_top_2 * N;
    elseif (x > 5 * N^2 && floor(x / N) == x / N)  % x in upper arm, rightmost column 
        row_from_bottom_2 = (x - 5 * N^2) / N;
        y = 5 * N^2 - N + row_from_bottom_2;
    else  
        y = x + 1;
    end
end

function y = left(x, N)
    if (x <= N^2 && floor((x - 1) / N) == (x - 1) / N) % x in lower half of lower arm, leftmost column
        row_from_bottom = (x - 1) / N + 1;
        y = 5 * N^2 - (row_from_bottom) * 3 * N + 1;
    elseif (x > N^2 && x <= 2 * N^2 && floor((x - 1) / N) == (x - 1) / N) % x in upper half of lower arm, rightmost column
        row_from_bottom_2 = ceil(x / N) - N;
        y = 2 * N^2 + row_from_bottom_2;
    elseif (x > 2 * N^2 && x <= 5 * N^2 && floor((x - 1 - 2 * N^2)/(3 * N)) == (x - 1 - 2 * N^2)/(3 * N)) % x in middle arm, leftmost column 
        row_from_top = ceil((5 * N^2 - x)/(3 * N));
        y = N * (row_from_top - 1) + 1;
    elseif (x > 5 * N^2 && floor((x - 1) / N) == (x - 1) / N)  % x in upper arm, rightmost column 
        row_from_top_2 = ceil((6 * N^2 - x) / N);
        y = 5 * N^2 - 3 * N + row_from_top_2;
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