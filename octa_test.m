n = 4;
m = 2;

k1 = pi*2*m/3;
k2 = 2*pi*n/(sqrt(3));

P1 = [-3/2,-sqrt(3)/2];
P2 = [1/2,-sqrt(3)/2];
P3 = [-1/2,sqrt(3)/2];
zfun = @(x1, x2) imag((exp(1i * (k1 * x1 + k2 * x2)) + exp(1i * ((-(k1 / 2) - sqrt(3) * k2 / 2) * x1 + ((sqrt(3) * k1 / 2) - k2 / 2) * x2))) ...
    + exp(1i * ((-(k1 / 2) + sqrt(3) * k2 / 2) * x1 - ((sqrt(3) * k1 / 2) + k2 / 2) * x2)));


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
trisurf(F,x,y,zfun(x,y),'Linestyle','none')

hold on
P1 = [1,0];
P2 = [0,0];
P3 = [1/2,-sqrt(3)/2];
% Create triangular mesh of triangla P1,P2,P3
[I,J]=meshgrid(linspace(0,1,n));
keep=I+J <= 1+eps;
IJ = [I(keep),J(keep)];
K = 1-sum(IJ,2);
F=delaunay(IJ);
IJK = [IJ,K];
XY=IJK*[P1;P2;P3];
x=XY(:,1);
y=XY(:,2);
trisurf(F,x,y,zfun(x,y),'Linestyle','none')

hold on
P1 = [1,0];
P2 = [3/2,-sqrt(3)/2];
P3 = [1/2,-sqrt(3)/2];
% Create triangular mesh of triangla P1,P2,P3
[I,J]=meshgrid(linspace(0,1,n));
keep=I+J <= 1+eps;
IJ = [I(keep),J(keep)];
K = 1-sum(IJ,2);
F=delaunay(IJ);
IJK = [IJ,K];
XY=IJK*[P1;P2;P3];
x=XY(:,1);
y=XY(:,2);
trisurf(F,x,y,zfun(x,y),'Linestyle','none')

hold on
P1 = [0,-sqrt(3)];
P2 = [-1/2,-sqrt(3)/2];
P3 = [1/2,-sqrt(3)/2];
% Create triangular mesh of triangla P1,P2,P3
[I,J]=meshgrid(linspace(0,1,n));
keep=I+J <= 1+eps;
IJ = [I(keep),J(keep)];
K = 1-sum(IJ,2);
F=delaunay(IJ);
IJK = [IJ,K];
XY=IJK*[P1;P2;P3];
x=XY(:,1);
y=XY(:,2);
trisurf(F,x,y,zfun(x,y),'Linestyle','none')

hold on
P1 = [-1,0];
P2 = [-2,0];
P3 = [-3/2,-sqrt(3)/2];
% Create triangular mesh of triangla P1,P2,P3
[I,J]=meshgrid(linspace(0,1,n));
keep=I+J <= 1+eps;
IJ = [I(keep),J(keep)];
K = 1-sum(IJ,2);
F=delaunay(IJ);
IJK = [IJ,K];
XY=IJK*[P1;P2;P3];
x=XY(:,1);
y=XY(:,2);
trisurf(F,x,y,zfun(x,y),'Linestyle','none')

xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');


