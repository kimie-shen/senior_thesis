n1=1;
n2=-2;

P1 = [-1/2,-sqrt(3)/2];
P2 = [3/2,-sqrt(3)/2];
P3 = [1/2,sqrt(3)/2];

zfun = @(x,y) cos(2*pi*(n1+n2)*x+2*sqrt(3)*pi*(-n1+n2/3)*y)+cos(2*pi*(n1-n2)*x+2*sqrt(3)*pi*(n1+n2/3)*y)+cos(4*pi*n1*x+4*pi*n2*y/sqrt(3));

% Create triangular mesh of triangles P1,P2,P3
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
P1 = [-5/2,-sqrt(3)/2];
P2 = [-1/2,-sqrt(3)/2];
P3 = [-3/2,sqrt(3)/2];
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
P1 = [-5/2,-sqrt(3)/2];
P2 = [-3/2,-sqrt(3)/2];
P3 = [-2,-sqrt(3)];
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
P1 = [-3/2,-sqrt(3)/2];
P2 = [-1/2,-sqrt(3)/2];
P3 = [-1,-sqrt(3)];
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
P1 = [-1/2,-sqrt(3)/2];
P2 = [1/2,-sqrt(3)/2];
P3 = [0,-sqrt(3)];
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
P1 = [1/2,-sqrt(3)/2];
P2 = [3/2,-sqrt(3)/2];
P3 = [1,-sqrt(3)];
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
P1 = [3/2,-sqrt(3)/2];
P2 = [5/2,-sqrt(3)/2];
P3 = [2,-sqrt(3)];
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
P2 = [0,0];
P3 = [-1/2,-sqrt(3)/2];
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
P2 = [2,0];
P3 = [3/2,-sqrt(3)/2];
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
P1 = [2,0];
P2 = [3,0];
P3 = [5/2,-sqrt(3)/2];
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
P1 = [2,0];
P2 = [1/2,-sqrt(3)/2];
P3 = [3/2,-sqrt(3)/2];
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
P2 = [0,0];
P3 = [-1/2,sqrt(3)/2];
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
P1 = [2,0];
P2 = [1,0];
P3 = [3/2,sqrt(3)/2];
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
P1 = [3,0];
P2 = [2,0];
P3 = [5/2,sqrt(3)/2];
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
P1 = [2,0];
P2 = [3/2,-sqrt(3)/2];
P3 = [5/2,-sqrt(3)/2];
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


