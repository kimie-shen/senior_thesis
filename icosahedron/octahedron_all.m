% Calculate r values for various N
clear;
close all
tolerance = 1e-6;  % Tolerance for degenerate energy levels
corner_potential = 1;
%N_nums = [5, 10, 15];
N_nums = [5, 7, 11, 13 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173];
%N_nums = 101;

r_array = zeros(size(N_nums, 1), 6);
size_array = zeros(size(N_nums, 1), 6);
solvable_prop = zeros(size(N_nums, 1), 2);

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
    H = zeros((N + 1)^2, (N + 1)^2); % Hamiltonian matrix
    
    for i = 1:(N + 1)^2
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
            %fprintf([num2str(lower(i, N)) ' \n'])
            H(i, lower(i, N)) = -(1 + 2/sqrt(3))/L^2;
        end
    end

    %% Add corner potential
    faces = [1, 3, 4];
    for i = 1:3
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
    bottom_corner = index_from_coord(1, 1, 2, N);
    top_right_corner = index_from_coord(N, (N + 1) / 2, 2, N);
    top_left_corner = index_from_coord(1, (N + 1) / 2, 2, N);

    H(bottom_corner, bottom_corner) = corner_potential; 
    H(top_right_corner, top_right_corner) = corner_potential;
    H(top_left_corner, top_left_corner) = corner_potential;

    %fprintf([num2str(bottom_corner) ' ' num2str(top_right_corner) ' ' num2str(top_left_corner) '\n'])


    % Generate matrix for C3
    C3 = zeros((N + 1)^2, (N + 1)^2);
    
    for i = 1:((N + 1)^2)
        %fprintf([num2str(i) ' ' num2str(c3_index(i, N)) '\n'])
        C3(i, c3_index(i, N)) = 1;
    end
    
    % Generate matrix for C2
    C2 = zeros((N + 1)^2, (N + 1)^2);
    
    for i = 1:((N + 1)^2)
        %fprintf([num2str(i) ' ' num2str(c2_index(i, N)) '\n'])
        C2(i, c2_index(i, N)) = 1;
    end

    % Generate matrix for sigma_d symmetry 
    sigma_d = zeros((N + 1)^2, (N + 1)^2);
    for i = 1: (N + 1)^2
        %fprintf([num2str(i) ' ' num2str(sigma_d_index(i, N)) '\n'])
        sigma_d(i, sigma_d_index(i, N)) = 1;
    end
    
    % Calculate eigenvectors and eigenvalues
    [eigenvectors, eigenvalues] = eig(H + ((1 + sqrt(5))/2) * sigma_d);
    fprintf('evals and evecs done \n')
    toc
    
    
    %% Calculate characters of symmetries
    tic;
    % Preallocate data array
    sigma2_evals = zeros((N + 1)^2, 1);
    c3_evals = zeros((N + 1)^2, 1);
    c2_evals = zeros((N + 1)^2, 1);
    
    
    for i = 1:((N + 1)^2)
        sigma2_evals(i) = (eigenvectors(:, i).' * sigma_d * eigenvectors(:, i));
        eigenvalues(i, i) = eigenvalues(i, i) - ((1 + sqrt(5))/2) * sigma2_evals(i);
    end
    fprintf('sigma_d done \n')
    
    for i = 1:((N + 1)^2)
        c3_evals(i) = (eigenvectors(:, i).' * C3 * eigenvectors(:, i));
    end
    fprintf('c3 done \n')

    for i = 1:((N + 1)^2)
        c2_evals(i) = (eigenvectors(:, i).' * C2 * eigenvectors(:, i));
    end
    fprintf('c2 done \n')

    toc;
    
    %% Analyze Degeneracies
    % Reorder energies and evecs
    energy_levels = zeros((N + 1)^2, 6);
    [energies, id] = sort(diag(eigenvalues));
    eigenvectors = eigenvectors(:, id);
    sigma2_evals_sorted = sigma2_evals(id);
    c3_evals_sorted = c3_evals(id);
    c2_evals_sorted = c2_evals(id);
    
    % Analyze for degeneracies
    energy = energies(1);
    degeneracy = 1;
    index = 1;
    trace_sigma2 = sigma2_evals_sorted(1);
    trace_c3 = c3_evals_sorted(1);
    trace_c2 = c2_evals_sorted(1);
    
    for i = 2:((N + 1)^2)
        if (abs(energies(i) - energy) < tolerance) % Next energy is degenerate
            degeneracy = degeneracy + 1;
            trace_sigma2 = trace_sigma2 + sigma2_evals_sorted(i);
            trace_c3 = trace_c3 + c3_evals_sorted(i);
            trace_c2 = trace_c2 + c2_evals_sorted(i);
        else % Next energy is new
            % Record stats for previous energy level
            energy_levels(index, 1) = energy;
            energy_levels(index, 2) = degeneracy;
            energy_levels(index, 3) = trace_sigma2;
            energy_levels(index, 4) = trace_c3;
            energy_levels(index, 5) = trace_c2;
            energy_levels(index, 6) = i - 1;
    
            % Reset variables for new energy level
            energy = energies(i);
            degeneracy = 1;
            index = index + 1;
            trace_sigma2 = sigma2_evals_sorted(i);
            trace_c3 = c3_evals_sorted(i);
            trace_c2 = c2_evals_sorted(i);
        end
    
        % Record energy level if we reach the end
        if (i == (N + 1)^2)
            energy_levels(index, 1) = energy;
            energy_levels(index, 2) = degeneracy;
            energy_levels(index, 3) = trace_sigma2;
            energy_levels(index, 4) = trace_c3;
            energy_levels(index, 5) = trace_c2;
            energy_levels(index, 6) = i;
        end
    end
    
    energy_levels = energy_levels(1:index, :);
    
    %% Separate the representation classes
    % Allocate space for each irrep
    elevels_a1 = zeros((N + 1)^2, 6);
    elevels_a2 = zeros((N + 1)^2, 6);
    elevels_e = zeros((N + 1)^2, 6);
    elevels_t1 = zeros((N + 1)^2, 6);
    elevels_t2 = zeros((N + 1)^2, 6);
    
    index_a1 = 1;
    index_a2 = 1;
    index_e = 1;
    index_t1 = 1;
    index_t2 = 1;
    
    % Round the energy_levels characters
    energy_levels_rounded = zeros(size(energy_levels, 1), 6);
    energy_levels_rounded(:, 2:5) = round(energy_levels(:, 2:5));
    energy_levels_rounded(:, 1) = energy_levels(:, 1);
    energy_levels_rounded(:, 6) = energy_levels(:, 6);
    
    % Fill in elevel info in each irrep
    for i = 1:size(energy_levels_rounded, 1)
        trace_e = energy_levels_rounded(i, 2);
        trace_sigma_d = energy_levels_rounded(i, 3);
        trace_c3 = energy_levels_rounded(i, 4);
        trace_c2 = energy_levels_rounded(i, 5);
    
        traces = [trace_e, trace_sigma_d, trace_c3, trace_c2];
    
        if (isequal(traces, [1, 1, 1, 1])) % A1
            elevels_a1(index_a1, :) = energy_levels_rounded(i, :);
            index_a1 = index_a1 + 1;
        elseif (isequal(traces, [1, -1, 1, 1])) % A2
            elevels_a2(index_a2, :) = energy_levels_rounded(i, :);
            index_a2 = index_a2 + 1;
        elseif (isequal(traces, [2, 0, -1, 2])) % E
            elevels_e(index_e, :) = energy_levels_rounded(i, :);
            index_e = index_e + 1;
        elseif (isequal(traces, [3, -1, 0, -1])) % T1
            elevels_t1(index_t1, :) = energy_levels_rounded(i, :);
            index_t1 = index_t1 + 1;
        elseif (isequal(traces, [3, 1, 0, -1])) %T2
            elevels_t2(index_t2, :) = energy_levels_rounded(i, :);
            index_t2 = index_t2 + 1;
        %elseif (isequal(traces, [6, 0, 0, -2])) % Accidental degeneracy T1 + T2
            %elevels_t1(index_t1, :) = energy_levels_rounded(i, :);
            %index_t1 = index_t1 + 1;
            %elevels_t2(index_t2, :) = energy_levels_rounded(i, :); 
            %index_t2 = index_t2 + 1;
        %elseif (isequal(traces, [3, 1, 0, 3])) % Accidental degeneracy E + A1
            %elevels_e(index_e, :) = energy_levels_rounded(i, :);
            %index_e = index_e + 1;
            %elevels_a1(index_a1, :) = energy_levels_rounded(i, :); 
            %index_a1 = index_a1 + 1;
        elseif (isequal(traces, [3, -1, 0, 3])) % Accidental degeneracy E + A2
            elevels_e(index_e, :) = energy_levels_rounded(i, :);
            index_e = index_e + 1;
            elevels_a2(index_a2, :) = energy_levels_rounded(i, :); 
            index_a2 = index_a2 + 1;
        %elseif (isequal(traces, [6, 0, 0, 6])) % Accidental degeneracy 2E + A1 + A2
            %elevels_e(index_e, :) = energy_levels_rounded(i, :);
            %elevels_e(index_e + 1, :) = energy_levels_rounded(i, :);
            %index_e = index_e + 2;
            %elevels_a1(index_a1, :) = energy_levels_rounded(i, :); 
            %index_a1 = index_a1 + 1;
            %elevels_a2(index_a2, :) = energy_levels_rounded(i, :); 
            %index_a2 = index_a2 + 1;
        else 
            fprintf([num2str(energy_levels_rounded(i, 1)) ' ' num2str(energy_levels_rounded(i, 2)) ' ' ...
                num2str(energy_levels_rounded(i, 3)) ' ' num2str(energy_levels_rounded(i, 4)) ' ' ...
                num2str(energy_levels_rounded(i, 5)) '\n'])
        end
    end
    
    % Remove extra rows of zeros
    if (index_a1 > 1)
        elevels_a1 = elevels_a1(1:(index_a1 - 1), :);
    end 
    
    if (index_a2 > 1)
        elevels_a2 = elevels_a2(1:(index_a2 - 1), :);
    end
    
    if (index_e > 1)
        elevels_e = elevels_e(1:(index_e - 1), :);
    end
    
    if (index_t1 > 1)
        elevels_t1 = elevels_t1(1:(index_t1 - 1), :);
    end
    
    if (index_t2 > 1)
        elevels_t2 = elevels_t2(1:(index_t2 - 1), :);
    end
    
  
    
    %% Find energy level spacings
    % Make vector of energy level spacings
    spacings_a1 = zeros(size(elevels_a1, 1) - 1, 1);
    spacings_a2 = zeros(size(elevels_a2, 1) - 1, 1);
    spacings_e = zeros(size(elevels_e, 1) - 1, 1);
    spacings_t1 = zeros(size(elevels_t1, 1) - 1, 1);
    spacings_t2 = zeros(size(elevels_t2, 1) - 1, 1);
    
    for i = 1:(size(elevels_a1, 1) - 1)
        spacings_a1(i) = abs(elevels_a1(i, 1) - elevels_a1(i+1, 1));
    end
    
    for i = 1:(size(elevels_a2, 1) - 1)
        spacings_a2(i) = abs(elevels_a2(i, 1) - elevels_a2(i+1, 1));
    end
    
    for i = 1:(size(elevels_e, 1) - 1)
        spacings_e(i) = abs(elevels_e(i, 1) - elevels_e(i+1, 1));
    end
    
    for i = 1:(size(elevels_t1, 1) - 1)
        spacings_t1(i) = abs(elevels_t1(i, 1) - elevels_t1(i+1, 1));
    end
    
    for i = 1:(size(elevels_t2, 1) - 1)
        spacings_t2(i) = abs(elevels_t2(i, 1) - elevels_t2(i+1, 1));
    end
    
    % Compute level spacing ratios
    r_a1 = LSR(spacings_a1);
    r_a2 = LSR(spacings_a2);
    r_e = LSR(spacings_e);
    r_t1 = LSR(spacings_t1);
    r_t2 = LSR(spacings_t2);

    % Fill in r array values
    r_array(index_r, :) = [N, r_a1, r_a2, r_e, r_t1, r_t2];
   
    % Fill in size of irreps
    size_array(index_size, :) = [N, size(elevels_a1, 1), size(elevels_a2, 1), size(elevels_e, 1), ...
        size(elevels_t1, 1), size(elevels_t2, 1)];

    % Find proportion of solvable states
    solvable_prop(index_size, :) = [size(elevels_a1, 1), size(elevels_a2, 1)] / ((N + 1)^2);

    %% Plot the level spacings histogram and show r values
    % Plot histograms
    figure(index_n)
    tiledlayout(1, 5,'TileSpacing', 'tight','Padding','Tight')
    bin_factor = 5;
    
    nexttile
    histogram(spacings_a1, ceil(size(elevels_a1, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('A1')
    subtitle(['r = ' num2str(r_a1) '; total = ' num2str(size_array(index_size, 2))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_a2, ceil(size(elevels_a2, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('A2')
    subtitle(['r = ' num2str(r_a2) '; total = ' num2str(size_array(index_size, 3))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_e, ceil(size(elevels_e, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('E')
    subtitle(['r = ' num2str(r_e) '; total = ' num2str(size_array(index_size, 4))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_t1, ceil(size(elevels_t1, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('T1')
    subtitle(['r = ' num2str(r_t1) '; total = ' num2str(size_array(index_size, 5))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_t2, ceil(size(elevels_t2, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('T2')
    subtitle(['r = ' num2str(r_t2) '; total = ' num2str(size_array(index_size, 6))])
    xlabel('Energy spacing')
    ylabel('Counts')

    set(figure(index_n),'position',[0,100,1500,200])
    
    % Increment indices
    index_n = index_n + 1;
    index_size = index_size + 1;
    index_r = index_r + 1;
end

%% Plot r as a function of size in each irrep
figure(index_n)
tiledlayout(1, 5, 'TileSpacing', 'tight', 'Padding', 'Tight')
ylow = 0.3;
yhigh = 0.7;

nexttile
plot(N_nums.', r_array(:, 2).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('A1')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_array(:, 3).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('A2')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_array(:, 4).', 'Linestyle', '-','Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('E')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_array(:, 5).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('T1')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_array(:, 6).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('T2')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow, yhigh])

set(figure(index_n),'position',[0,100,1500,200])

%% Plot proportions of solvable states
figure(index_n + 1)
plot(N_nums, sum(solvable_prop, 2), 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Proportion of states which are solvable')
xlabel('N')
yline(1/12, '--', '1/12')
ylim([0.0, 0.1])
ylabel('Proportion')


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

        y = mod(index - 1, row_length) + 1;
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
    row_w = row_width(index, N);

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