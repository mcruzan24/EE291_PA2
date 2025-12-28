%section 2: boundary value problems

clear; clc; close all;

%inital conditions
w = 4;
T = 10000;
y0 = 25;
yL = 25;

%define ode
function dydx = odefun(x, y)
    w = 4; T = 10000;
    dydx = [y(2); w/T * sqrt(1 + y(2)^2)];
end

%define boundary conditions
function res = bcfun(ya, yb)
    y0 = 25; yL = 25;
    res = [ya(1)-y0; yb(1)-yL];
end

%q1: max wire drop when posts are 100 m apart
L1 = 100;
xmesh = linspace(0, L1, 10);
solinit = bvpinit(xmesh, [y0 0]);

%solve bvp
sol1 = bvp4c(@odefun, @bcfun, solinit);

%get solution
x1 = linspace(0, L1, 100);
y1 = deval(sol1, x1); %deval evaluates differential equations 

%calculate max drop
drop1 = max(y1(1,:)) - min(y1(1,:));
disp("Wire drop distance for 100 m: ");
disp(drop1);

%q2: max distance for 2 posts
L = 100;
step = 1;
while true
    xmesh_1 = linspace(0, L, 10);
    solinit_1 = bvpinit(xmesh_1, [y0 0]);
    sol = bvp4c(@odefun, @bcfun, solinit_1);
    ytemp = deval(sol, linspace(0,L,100));
    max_drop = max(ytemp(1,:)) - min(ytemp(1,:));
    if max_drop > 3 %if max drop is greater than 3, 
        Lmax = L - step; %its the previous value, break
        break
    end
    L = L + step; %if not, keep going
end
disp("Maximum post distance for drop <= 3 m: ");
disp(Lmax);

x2 = linspace(0, Lmax, 100);
xmesh_2 = linspace(0, Lmax, 10);
solinit_2 = bvpinit(xmesh_2, [y0 0]);
sol2 = bvp4c(@odefun, @bcfun, solinit_2);
y2 = deval(sol2, x2);

%q3: plot
hold on;
plot(x1, y1(1,:)); 
plot(x2, y2(1,:));
hold off;

xlabel("Horizontal Span (m)");
ylabel("Vertical Span (m)");

%q4: compare to PA2 solution
disp("The value I got for max drop and max distance in PA2 and PA3 were virtually the same.");
disp("The value I got from PA2 using a numeric ODE solver as well as the shooting method");
disp("was 0.5009 m and 244.8547 m, while the value I got from PA3 using the built in");
disp("MATLAB functions was 0.4999 m and 244 m. This shows how effective the numerical methods");
disp("implemented in PA2 are.")
