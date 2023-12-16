% Compute spectral form factor functions
logt = linspace(-12, 2, 10000);

tic;
% Allocate variables
logkt_a1g = zeros(1, size(logt, 2));
logkt_a1u = zeros(1, size(logt, 2));
logkt_a2g = zeros(1, size(logt, 2));
logkt_a2u = zeros(1, size(logt, 2));
logkt_eg = zeros(1, size(logt, 2));
logkt_eu = zeros(1, size(logt, 2));
logkt_t1g = zeros(1, size(logt, 2));
logkt_t1u = zeros(1, size(logt, 2));
logkt_t2g = zeros(1, size(logt, 2));
logkt_t2u = zeros(1, size(logt, 2));

% Fill in log k(t) values
for i = 1:size(logt, 2)
    logkt_a1g(1, i) = SFF(elevels_a1g,  exp(logt(1, i)));
    logkt_a1u(1, i) = SFF(elevels_a1u,  exp(logt(1, i)));
    logkt_a2g(1, i) = SFF(elevels_a2g,  exp(logt(1, i)));
    logkt_a2u(1, i) = SFF(elevels_a2u,  exp(logt(1, i)));
    logkt_eg(1, i) = SFF(elevels_eg,  exp(logt(1, i)));
    logkt_eu(1, i) = SFF(elevels_eu,  exp(logt(1, i)));
    logkt_t1g(1, i) = SFF(elevels_t1g,  exp(logt(1, i)));
    logkt_t1u(1, i) = SFF(elevels_t1u,  exp(logt(1, i)));
    logkt_t2g(1, i) = SFF(elevels_t2g,  exp(logt(1, i)));
    logkt_t2u(1, i) = SFF(elevels_t2u,  exp(logt(1, i)));
end

% Plot log k(t) vs log t
figure(1)
tiledlayout(2,5)
nexttile
scatter(logt, logkt_a1g(1, :), '.')
xlabel('log t')
ylabel('log K(t)')
title('A1g')

nexttile
scatter(logt, logkt_a2g(1, :), '.')
xlabel('log t')
ylabel('log K(t)')
title('A2g')

nexttile
scatter(logt, logkt_eg(1, :), '.')
xlabel('log t')
ylabel('log K(t)')
title('Eg')

nexttile
scatter(logt, logkt_t1g(1, :), '.')
xlabel('log t')
ylabel('log K(t)')
title('T1g')

nexttile
scatter(logt, logkt_t2g(1, :), '.')
xlabel('log t')
ylabel('log K(t)')
title('T2g')

nexttile
scatter(logt, logkt_a1u(1, :), '.')
xlabel('log t')
ylabel('log K(t)')
title('A1u')

nexttile
scatter(logt, logkt_a2u(1, :), '.')
xlabel('log t')
ylabel('log K(t)')
title('A2u')

nexttile
scatter(logt, logkt_eu(1, :), '.')
xlabel('log t')
ylabel('log K(t)')
title('Eu')

nexttile
scatter(logt, logkt_t1u(1, :), '.')
xlabel('log t')
ylabel('log K(t)')
title('T1u')

nexttile
scatter(logt, logkt_t2u(1, :), '.')
xlabel('log t')
ylabel('log K(t)')
title('T2u')

set(figure(1),'position',[0,100,1500,400])

toc

function logkt = SFF(elevels, t)
    num_elevels = size(elevels, 1);
    sum = 0;

    for i = 1:num_elevels
        sum = sum + exp(1i * elevels(i, 1) * t);
    end

    logkt = log(abs(sum)^2);
end
