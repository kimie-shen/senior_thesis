m = 4;
n = 0;

%m and n must have same parity

k1=2*pi*n;
k2=2*pi*m/sqrt(3);

P1 = [-3/2,-sqrt(3)/2];
P2 = [1/2,-sqrt(3)/2];
P3 = [-1/2,sqrt(3)/2];

% Sin wavefunction
%zfun = @(x1,x2) sin(k1*x1+k2*x2) + sin((-k1-sqrt(3)*k2)*x1/2 + (k1*sqrt(3)-k2)*x2/2)+ sin((-k1+sqrt(3)*k2)*x1/2 + (-k1*sqrt(3)-k2)*x2/2);

% Cos wavefunction
zfun = @(x1,x2) cos(k1*x1+k2*x2) + cos((-k1-sqrt(3)*k2)*x1/2 + (k1*sqrt(3)-k2)*x2/2)  + cos((-k1+sqrt(3)*k2)*x1/2 + (-k1*sqrt(3)-k2)*x2/2);

% Cosine degeneracy
%zfun = @(x1,x2) cos(k1*x1+k2*x2) + cos((-k1-sqrt(3)*k2)*x1/2 + (k1*sqrt(3)-k2)*x2/2) + cos((-k1+sqrt(3)*k2)*x1/2 + (-k1*sqrt(3)-k2)*x2/2) ...
%- (cos(k1*x1-k2*x2) + cos((-k1+sqrt(3)*k2)*x1/2 + (k1*sqrt(3)+k2)*x2/2) + cos((-k1-sqrt(3)*k2)*x1/2 + (-k1*sqrt(3)+k2)*x2/2));

% Sine degneracy
%zfun = @(x1,x2) sin(k1*x1+k2*x2) + sin((-k1-sqrt(3)*k2)*x1/2 + (k1*sqrt(3)-k2)*x2/2) + sin((-k1+sqrt(3)*k2)*x1/2 + (-k1*sqrt(3)-k2)*x2/2) ...
%- (sin(k1*x1-k2*x2) + sin((-k1+sqrt(3)*k2)*x1/2 + (k1*sqrt(3)+k2)*x2/2) + sin((-k1-sqrt(3)*k2)*x1/2 + (-k1*sqrt(3)+k2)*x2/2));

% Full Degeneracy
%zfun = @(x1,x2) cos(k1*x1+k2*x2) + cos((-k1-sqrt(3)*k2)*x1/2 + (k1*sqrt(3)-k2)*x2/2) + cos((-k1+sqrt(3)*k2)*x1/2 + (-k1*sqrt(3)-k2)*x2/2) ...
%- (cos(k1*x1-k2*x2) + cos((-k1+sqrt(3)*k2)*x1/2 + (k1*sqrt(3)+k2)*x2/2) + cos((-k1-sqrt(3)*k2)*x1/2 + (-k1*sqrt(3)+k2)*x2/2)) + ...
%+ (sin(k1*x1+k2*x2) + sin((-k1-sqrt(3)*k2)*x1/2 + (k1*sqrt(3)-k2)*x2/2) + sin((-k1+sqrt(3)*k2)*x1/2 + (-k1*sqrt(3)-k2)*x2/2))  ...
%- (sin(k1*x1-k2*x2) + sin((-k1+sqrt(3)*k2)*x1/2 + (k1*sqrt(3)+k2)*x2/2) + sin((-k1-sqrt(3)*k2)*x1/2 + (-k1*sqrt(3)+k2)*x2/2));

% Two degeneracy
%zfun = @(x1,x2) sin(k1*x1+k2*x2) + sin((-k1-sqrt(3)*k2)*x1/2 + (k1*sqrt(3)-k2)*x2/2)+ sin((-k1+sqrt(3)*k2)*x1/2 + (-k1*sqrt(3)-k2)*x2/2)...
%    - (cos(k1*x1+k2*x2) + cos((-k1-sqrt(3)*k2)*x1/2 + (k1*sqrt(3)-k2)*x2/2)  + cos((-k1+sqrt(3)*k2)*x1/2 + (-k1*sqrt(3)-k2)*x2/2));



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
view(0, 90);
ylim([-2*sqrt(3)/2, sqrt(3)/2])
daspect([1 1 0.2])

