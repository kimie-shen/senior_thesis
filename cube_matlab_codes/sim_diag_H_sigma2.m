% Obtain eigenvectors and eigenvalues of H and sigma2
tic;
[eigenvectors, eigenvalues] = eig(H + ((1 + sqrt(5))/2)*sigma2);
toc

%  Takes 45 min for N=75


