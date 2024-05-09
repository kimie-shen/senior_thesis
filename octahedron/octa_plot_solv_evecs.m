n1 = 4/3;
n2 = -2;


P1 = [1,0];
P2 = [0,sqrt(3)];
P3 = [-1,0];

%zfun = @(x,y) cos(pi*(n1+n2)*x+sqrt(3)*pi*(-n1+n2/3)*y);
%zfun = @(x,y) cos(-2*pi*n1*x-2*pi*n2*y/sqrt(3));
%zfun = @(x,y) cos(pi*(n1-n2)*x+sqrt(3)*pi*(n1+n2/3)*y);
%zfun = @(x,y) cos(pi*(n1+n2)*x+sqrt(3)*pi*(-n1+n2/3)*y)+cos(pi*(n1-n2)*x+sqrt(3)*pi*(n1+n2/3)*y)+cos(-2*pi*n1*x-2*pi*n2*y/sqrt(3));
%zfun = @(x,y) cos(pi*(n1+n2)*x+sqrt(3)*pi*(-n1+n2/3)*y)+cos(pi*(n1-n2)*x+sqrt(3)*pi*(n1+n2/3)*y)+cos(2*pi*n1*x+2*pi*n2*y/sqrt(3)) +...
   % cos(pi*(n1-n2)*x+sqrt(3)*pi*(-n1-n2/3)*y)+cos(pi*(n1+n2)*x+sqrt(3)*pi*(n1-n2/3)*y)+cos(2*pi*n1*x-2*pi*n2*y/sqrt(3));

zfun = @(x,y) sin(pi*(n1+n2)*x+sqrt(3)*pi*(-n1+n2/3)*y)+sin(pi*(n1-n2)*x+sqrt(3)*pi*(n1+n2/3)*y)+sin(-2*pi*n1*x-2*pi*n2*y/sqrt(3));
%zfun = @(x,y) sin(2*pi*n1*x+2*pi*n2*y/sqrt(3));

% Create triangular mesh of triangles P1,P2,P3

figure(3)

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
P1 = [-1,0];
P2 = [-3/2, sqrt(3)/2];
P3 = [-1/2, sqrt(3)/2];
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
P1 = [0, 0];
P2 = [1, 0];
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
P2 = [3/2, sqrt(3)/2];
P3 = [1/2,sqrt(3)/2];
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

ylim([-sqrt(3)/2, sqrt(3)])
view(0, 90);
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');


