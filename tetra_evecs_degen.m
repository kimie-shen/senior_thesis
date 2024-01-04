n = 2;
m = 3;

k1=2*pi*n;
k2=2*pi*(2*m-n)/sqrt(3);

P1 = [0,0];
P2 = [2,0];
P3 = [1,sqrt(3)];


% Create triangular mesh of triangla P1,P2,P3
n=200;
[I,J]=meshgrid(linspace(0,1,n));
keep=I+J <= 1+eps;
IJ = [I(keep),J(keep)];
K = 1-sum(IJ,2);
F=delaunay(IJ);
IJK = [IJ,K];
XY=IJK*[P1;P2;P3];
x=XY(:,1);
y=XY(:,2);

tiledlayout(1,4, 'TileSpacing','Tight')
nexttile
zfun = @(x1,x2) cos(pi*(m+n)*x1+pi*(m-n)*x2/sqrt(3)) ...
    + cos(pi*((-m-n)+m)*x1+pi*((-m-n)-m)*x2/sqrt(3)) ...
    + cos(pi*(n+(-m-n))*x1+pi*(n-(-m-n))*x2/sqrt(3));
trisurf(F,x,y,zfun(x,y),'Linestyle','none')
xlabel('$x_1$','Interpreter','Latex');
ylabel('$x_2$','Interpreter','Latex');
title('$\psi_{n,m} + \psi_{-m-n,m}+\psi_{n,-m-n}$','Interpreter','Latex');
view(0,90);
daspect([1 1 0.2])
ylim([0,sqrt(3)])

nexttile
zfun = @(x1,x2) cos(pi*(n+m)*x1+pi*(n-m)*x2/sqrt(3)) ...
    + cos(pi*((-n-m)+n)*x1+pi*((-n-m)-n)*x2/sqrt(3)) ...
    + cos(pi*(m+(-n-m))*x1+pi*(m-(-n-m))*x2/sqrt(3));
trisurf(F,x,y,zfun(x,y),'Linestyle','none')
xlabel('$x_1$','Interpreter','Latex');
ylabel('$x_2$','Interpreter','Latex');
title('$(n \leftrightarrow m)$','Interpreter','Latex');
view(0,90);
daspect([1 1 0.2])
ylim([0,sqrt(3)])

nexttile
zfun = @(x1,x2) cos(pi*(m+n)*x1+pi*(m-n)*x2/sqrt(3)) ...
    + cos(pi*((-m-n)+m)*x1+pi*((-m-n)-m)*x2/sqrt(3)) ...
    + cos(pi*(n+(-m-n))*x1+pi*(n-(-m-n))*x2/sqrt(3)) ...
    + cos(pi*(n+m)*x1+pi*(n-m)*x2/sqrt(3)) ...
    + cos(pi*((-n-m)+n)*x1+pi*((-n-m)-n)*x2/sqrt(3)) ...
    + cos(pi*(m+(-n-m))*x1+pi*(m-(-n-m))*x2/sqrt(3));
trisurf(F,x,y,zfun(x,y),'Linestyle','none')
xlabel('$x_1$','Interpreter','Latex');
ylabel('$x_2$','Interpreter','Latex');
title('$\psi_{n,m} + \psi_{-m-n,m}+\psi_{n,-m-n} + (n \leftrightarrow m)$','Interpreter','Latex');
view(0,90);
daspect([1 1 0.2])
ylim([0,sqrt(3)])

nexttile
zfun = @(x1,x2) cos(pi*(m+n)*x1+pi*(m-n)*x2/sqrt(3)) ...
    + cos(pi*((-m-n)+m)*x1+pi*((-m-n)-m)*x2/sqrt(3)) ...
    + cos(pi*(n+(-m-n))*x1+pi*(n-(-m-n))*x2/sqrt(3)) ...
    - cos(pi*(n+m)*x1+pi*(n-m)*x2/sqrt(3)) ...
    - cos(pi*((-n-m)+n)*x1+pi*((-n-m)-n)*x2/sqrt(3)) ...
    - cos(pi*(m+(-n-m))*x1+pi*(m-(-n-m))*x2/sqrt(3));
trisurf(F,x,y,zfun(x,y),'Linestyle','none')
xlabel('$x_1$','Interpreter','Latex');
ylabel('$x_2$','Interpreter','Latex');
title('$\psi_{n,m} + \psi_{-m-n,m}+\psi_{n,-m-n} - (n \leftrightarrow m)$','Interpreter','Latex');
view(0,90);
daspect([1 1 0.2])
ylim([0,sqrt(3)])



