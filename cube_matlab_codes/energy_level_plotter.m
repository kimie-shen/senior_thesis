
% Plot energy levels
figure(1)
for i = 1:length(energy_levels(:,1))
    x1 = [0.5, 1.5];
    x2 = [1.5, 2.5];
    x3 = [2.5, 3.5];
    x4 = [3.5, 4.5];
    y = [-energy_levels(i, 1), -energy_levels(i, 1)];

    if (energy_levels(i, 2) == 1)
        plot(x1, y,'-r')
        hold on
    end

    if (energy_levels(i, 2) == 2)
        plot(x2, y, 'g')
        hold on
    end

    if (energy_levels(i, 2) == 3)
        plot(x3, y, 'b')
        hold on
    end

    if (energy_levels(i, 2) > 3)
        plot (x4, y, 'm')
    end

end

%for i = 0:30
    %for j = i:30
        %E = 2 * pi^2 * (i^2 + j^2);
        %x = [0, 1];
        %y = [E, E];
        %plot(x, y, '-k');
        %hold on
    %end
%end

hold off