n1 = 2;
n2 = 2;

f1 = @(x,y) sin(pi*(x*n1+y*n2))*sin(pi*(y*n1-x*n2));
fsurf(f1,[-2 2 -1 0],'edgeColor','none')
hold on
f2 = @(x,y) sin(pi*(x*n1+y*n2))*sin(pi*(y*n1-x*n2));
fsurf(f2,[-1 0 0 1],'edgeColor','none')
f3 = @(x,y) sin(pi*(x*n1+y*n2))*sin(pi*(y*n1-x*n2));
fsurf(f3,[-1 0 -2 -1],'edgeColor','none')
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
hold off