n1 = 3;
n2 = 0;

f1 = @(x,y) cos(x*2*pi)+cos(y*2*pi);
fsurf(f1,[-2 2 -1 0],'edgeColor','none')
hold on
f2 = @(x,y) cos(x*2*pi)+cos(y*2*pi);
fsurf(f2,[-1 0 0 1],'edgeColor','none')
f3 = @(x,y) cos(x*2*pi)+cos(y*2*pi);
fsurf(f3,[-1 0 -2 -1],'edgeColor','none')
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
hold off
