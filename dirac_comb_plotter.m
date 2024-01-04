x1=0.6;
x2=0.3;

for n1 = -3:1:3
    for n2 = -3:1:3
        for j = 0:1:1
            plot(-x1-(2*n1+j),-x2-(2*n2+j),'o','MarkerFaceColor','r', 'MarkerEdgeColor', 'r')
            plot(x1+(2*n1+j),x2 + (2*n2+j),'o','MarkerFaceColor','m', 'MarkerEdgeColor', 'm')
            plot(-x2-(2*n1+j),x1+(2*n2+j),'o','MarkerFaceColor','g', 'MarkerEdgeColor', 'g')
            plot(x2+(2*n1+j), -x1-(2*n2+j),'o','MarkerFaceColor','b', 'MarkerEdgeColor', 'b')
            hold on
        end
    end
end

% Plot outline of unfolded cube
yline([-1,0,1,2])
xline([-1,0,1,2])

xlim([-2,2])
ylim([-2,2])
hold off