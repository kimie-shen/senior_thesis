% Check if discretized wavefunction is eigenvector of tight-binding matrix

num_sites = 100;
n1 = 1;
n2 = 2;

pred_energy = 2 * pi^2 * (n1^2 + n2^2);

E_wvfn = zeros(1, 6 * num_sites^2);

for i = 1:(6 * num_sites^2)
    E_wvfn(1, i) = pred_energy * wavefn_discrete(i);
end

tic; % Takes about 20 sec for num_sites=100
test = H * wavefn_discrete - E_wvfn;
toc


