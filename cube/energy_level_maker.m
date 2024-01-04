% Take energy eigenvalues and remove degeneracies 
tolerance = 1e-8;   % Max difference in energies considered degenerate
N = 75;

tic; 
%9.5 hours for N = 75
% Preallocate data array
sigma2_evals = zeros(6 * N^2, 1);
c3_evals = zeros(6 * N^2, 1);
inv_evals = zeros(6 * N^2, 1);
c2_evals = zeros(6 * N^2, 1);

% Subtract off sigma2 eigenvalues from eigenvalues 
% (START WITH UNSORTED EVALS AND EVECS)
for i = 1:(6 * N^2)
    inv_evals(i) = (eigenvectors(:, i).' * inversion * eigenvectors(:, i));
end
fprintf('inv done')

for i = 1:(6 * N^2)
    sigma2_evals(i) = (eigenvectors(:, i).' * sigma2 * eigenvectors(:, i));
    eigenvalues(i, i) = eigenvalues(i, i) - ((1 + sqrt(5))/2) * sigma2_evals(i);
end
fprintf('sigma_d done')

for i = 1:(6 * N^2)
    c3_evals(i) = (eigenvectors(:, i).' * C3 * eigenvectors(:, i));
end
fprintf('c3 done')

for i = 1:(6 * N^2)
    c2_evals(i) = (eigenvectors(:, i).' * C2 * eigenvectors(:, i));
end
fprintf('c2 done')

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

toc
