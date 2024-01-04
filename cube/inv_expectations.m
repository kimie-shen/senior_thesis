N = 75;

tic;
% Preallocate data array
inv_exp = zeros(6 * N^2, 1);

% Subtract off sigma2 eigenvalues from eigenvalues
for i = 1:(6 * N^2)
    inv_exp(i) = (eigenvectors(:, i).' * C3 * eigenvectors(:, i));
    fprintf([num2str(i) '/' num2str(6 * N^2) '\n'])
end

toc