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

function r = LSR(spacings)
    r_n = zeros(size(spacings, 1) - 1, 1);
    
    for i = 1:(size(spacings, 1) - 1)
        r_n(i) = min(spacings(i), spacings(i + 1)) ...
            / max(spacings(i), spacings(i + 1));
    end

    r = mean(r_n);
end
