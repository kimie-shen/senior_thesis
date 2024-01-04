N = 75;

tic;
% Preallocate data array
c3_exp = zeros(6 * N^2, 1);

% Subtract off sigma2 eigenvalues from eigenvalues
for i = 1:(6 * N^2)
    c3_exp(i) = (eigenvectors(:, i).' * C3 * eigenvectors(:, i));
end

toc