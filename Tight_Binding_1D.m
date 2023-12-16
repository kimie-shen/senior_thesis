clear
N = 100; % Total number of sites
L = 1; % Spacing
H = zeros(N, N); % Hamiltonian matrix

% Set matrix values
for i = 1:N
    H(i, i) = -2/L^2;
    H(i, mod(i-2, N) + 1) = 1/L^2;
    H(i, mod(i, N) + 1) = 1/L^2;
end

% Find eigenvalues
e = eig(H);

%spacings = zeros(N - 1, 1);
%for i=1:(N - 1)
    %spacings(i) = (abs(e(i)) - abs(e(i+1)));
%end

histogram(spacings)

% Plot eigenvalues
for i = 1:N
    yline(-e(i), '-r')
    hold on
end

hold on
for i = 0:(N/4)
    yline((2*pi*i)^2/(L*N)^2, '-b');
end

hold off
