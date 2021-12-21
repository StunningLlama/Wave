function main()

%schrodinger1D();
schrodinger2D();


%preparesimulationparameters2D([0 1], [0 1], 0.2, 1);
end

function saveimg(name)
%saveas(gcf, "C:\Users\Brandon\Desktop\2218dia\" + name);
end

function [X A L p intnodes] = preparesimulationparameters1D(xint, ds)
X = [xint(1): ds: xint(2)];
n = size(X, 2);
grid = ones(1,n);
grid(1) = 0;
grid(n) = 0;
p = find(grid);
intnodes = length(p);
L = zeros(1,n);
L(2:n-1) = 1:n-2;
e = ones(n-2,1);
A = spdiags([e -2*e e],[-1 0 1],n-2,n-2);
end

function [X Y A L Ldisp p intnodes] = preparesimulationparameters2D(xrange, yrange, ds, domain)
xvec = [xrange(1): ds: xrange(2)];
yvec = [xrange(1): ds: xrange(2)];
[X Y] = meshgrid(xvec, yvec);

n = size(xvec, 2);
m = size(yvec, 2);
grid = zeros(m,n);
grid(2:m-1, 2:n-1) = ones(m-2,n-2);

grid = grid.*(~((X < 2.1 & X > 1.9) & (Y < 1.7 | (Y > 1.8 & Y < 2.2) | Y > 2.)));

p = find(grid);
intnodes = length(p);
L = zeros(n);
L(p) = 1:intnodes;
Ldisp = L;
A = -delsq(L);
%L
if (domain == 1)
A = joinboundary(A, L(2,2:n-1),L(m-1,2:n-1));
A = joinboundary(A, L(2:m-1,2),L(2:m-1,n-1));
Ldisp(1,2:n-1) = Ldisp(m-1,2:n-1);
Ldisp(m,2:n-1) = Ldisp(2,2:n-1);
Ldisp(2:m-1,1) = Ldisp(2:m-1,n-1);
Ldisp(2:m-1,n) = Ldisp(2:m-1,2);
Ldisp(1,1)=Ldisp(2,2);
Ldisp(m,1)=Ldisp(m-1,2);
Ldisp(1,n)=Ldisp(2,n-1);
Ldisp(m,n)=Ldisp(m-1,n-1);
end
%full(A)
end

function A = joinboundary(A, boundA, boundB)
boundA = boundA(boundA~=0);
boundB = boundB(boundB~=0);
A(sub2ind(size(A), boundA, boundB)) = 1;
A(sub2ind(size(A), boundB, boundA)) = 1;
end

function plot1Dsol(X, L, u, titlestr, xlabelstr, ylabelstr)
utmp = [0; u];
U = utmp(L+1);
plot(X.',U);
title(titlestr);
xlabel(xlabelstr);
ylabel(ylabelstr);
end

function plot1Dsolwithboundary(X, L, u, bv, titlestr, xlabelstr, ylabelstr)
utmp = [0; u];
U = utmp(L+1) + bv;
plot(X,U);
title(titlestr);
xlabel(xlabelstr);
ylabel(ylabelstr);
end

function plot2Dsol(X, Y, L, u, titlestr)
utmp = [0; u];
U = utmp(L+1);
[c, h] = contourf(X,Y,U); clabel(c,h);  colorbar;%caxis([0 1]);
C=caxis;
caxis([C(1),C(2)])
hold on;
contour(X, Y, L==0);
hold off;
%surf(X,Y,U);
title(titlestr);
end

function plot2Dsolwithboundary(X, Y, L, u, bv, titlestr)
utmp = [0; u];
U = utmp(L+1) + bv;
[c, h] = contourf(X,Y,U); clabel(c,h), colorbar;
title(titlestr);
end

function schrodinger1D()
ds = 0.01;
dt = ds*ds*0.2;

hbar = 1;
m = 1;

Tfinal = 2.0;
frames = 500;

sigma =(dt/ds^2);
[X A L p intnodes] = preparesimulationparameters1D([0 2], ds);
xvals = X(p);
%b = 1000.0*(xvals > 1.0).';
b = (1000.0*(xvals-1.0).^2).';
u = exp(-100.0*(xvals-0.5).^2 + 0*i*(xvals-0.5)).';
%u = u/(sum(conj(u).*u)*ds);
A(1, length(p)) = 1;
A(length(p),1) = 1;
Aprime = eye(length(p)) + i*hbar*dt*diag(b) - (i*hbar/(2*m))*sigma*A;
Aprime = sparse(Aprime);
set(gcf,'CurrentCharacter', ' ');
for ii = [1:1:frames]
    for j = [0:dt:Tfinal/(frames)]
        u = Aprime\u;
        u = u/sqrt(sum(conj(u).*u)*ds);
    end
    t = ii/frames*Tfinal;
    plot1Dsol(X, L, abs(u), "1(b) Solution, u(x) for t = " + num2str(t), "x", "u");
    axis([0 2 0 10]);
    if (get(gcf,'CurrentCharacter') == 's')
        break;
    end
    pause(0.1);
end
end


function schrodinger2D()
ds = 0.05;
dt = ds*ds*0.25;

hbar = 1;
m = 1;

Tfinal = 2.0;
frames = 500;

sigma =(dt/ds^2);
[X Y A L Ldisp p intnodes] = preparesimulationparameters2D([0, 4], [0, 4], ds, 0);
xvals = X(p);
yvals = Y(p);
u = exp(-10.0*((xvals-0.5).^2 + (yvals-2.0).^2) + 50*i*(xvals-0.5));
b = 000.0*(xvals > 1);
Aprime = eye(length(p)) + i*hbar*dt*diag(b) - (i*hbar/(2*m))*sigma*A;
Aprime = sparse(Aprime);
spy(Aprime)
pause();
set(gcf,'CurrentCharacter', ' ');
for ii = [1:1:frames]
    for j = [0:dt:Tfinal/(frames)]
        u = Aprime\u;
        u = u/sqrt(sum(conj(u).*u)*ds);
    end
    t = ii/frames*Tfinal;
    plot2Dsol(X, Y, Ldisp, abs(u), "2(b) Solution, u(x,y) for t = " + num2str(t));
    
    if (get(gcf,'CurrentCharacter') == 's')
        break;
    end
    pause(0.1);
end
end
