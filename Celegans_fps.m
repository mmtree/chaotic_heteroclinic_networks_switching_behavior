%% create a grid of fps for C. elegans

clear all; close all; clc;

% linearized dynamics around fixed point 1
A1 = [-1 0;
    0 1/2];

% linearized dynamics around fixed point 2
A2 = [-1 0;
    0 1/2];

% location of fp 1
fp1_loc = [-1; -1];

% location of fp 2
fp2_loc = [1; 1];

b1 = -A1*fp1_loc;
b2 = -A2*fp2_loc;

%%

% heteroclinic network dynamics
x0_vec = 4*(rand(2,100)-[0.5; 0.5])+[0;0];
dt = 0.05;
tspan = 0:dt:10;

figure(); hold on;
for j = 1:length(x0_vec)
    x0 = x0_vec(:,j);
    [t,x] = ode45(@(t,x) F(t,x,A1,A2, b1,b2), tspan, x0);
    plot(x(1,1),x(1,2),'color','c');
    plot(x(:,1),x(:,2),'color','c');
end

scatter(fp1_loc(1), fp1_loc(2), 100,'k','*','LineWidth',1.5)
scatter(fp2_loc(1), fp2_loc(2), 100,'k','*','LineWidth',1.5)
xlabel('x')
ylabel('y')
title('heteroclinic network connecting 2 fixed points')

xlim([-3 3]);
ylim([-3 3]);


%%

% weighting function
figure();
x = -5:0.1:5;
y = -5:0.1:5;
[X,Y] = meshgrid(x,y);
Z = W11(X,Y);
surf(X,Y,Z);
xlabel('x');
ylabel('y');
title('weighting function')



%% functions


function w1 = W1(x,y)  % weighting fp1
    s = 10;
    w1 = -(1/4)*(tanh(s*y)-tanh(s*(y+2))).*(-tanh(s*x)-tanh(-s*(x+2)));
end

function w2 = W2(x,y) % weighting fp2
    s = 10;
    w2 = (1/4)*(tanh(s*y)-tanh(s*(y-2))).*(tanh(s*x)-tanh(s*(x-2)));
end

function w12 = W12(x,y) % weighting for fp1 -> fp2
    s = 10;
    w12 = -(1/4)*(tanh(s*y)-tanh(s*(y-2))).*(tanh(s*x)-tanh(s*(x+2)));
end

function w21 = W21(x,y) % weighting for fp2 -> fp1
    s = 10;
    w21 = -(1/4)*(tanh(s*y)-tanh(s*(y+2))).*(tanh(s*x)-tanh(s*(x-2)));
end

function w22 = W22(x,y) % weighting for fp2 -> fp1
    s = 10;
    w22 = (1/4)*(tanh(s*(y-2))-tanh(s*(y-4))).*(tanh(s*(x-2))-tanh(s*(x-4))) +....
        (1/4)*(tanh(s*(y-2))-tanh(s*(y-4))).*(tanh(s*(x))-tanh(s*(x-2))) +....
        (1/4)*(tanh(s*(y))-tanh(s*(y-2))).*(tanh(s*(x-2))-tanh(s*(x-4)));
end


function w11 = W11(x,y) % weighting for fp2 -> fp1
    s = 10;
    w11 = (1/4)*(tanh(s*(y+2))-tanh(s*(y+4))).*(tanh(s*(x+2))-tanh(s*(x+4))) +....
        (1/4)*(tanh(s*(y+2))-tanh(s*(y+4))).*(tanh(s*(x))-tanh(s*(x+2))) +....
        (1/4)*(tanh(s*(y))-tanh(s*(y+2))).*(tanh(s*(x+2))-tanh(s*(x+4)));
end



function dxdt = F(t,x,A1,A2,b1,b2)
    alpha1 = W1(x(1),x(2));
    alpha2 = W2(x(1),x(2));
    
    alpha12 = W12(x(1),x(2));
    alpha21 = W21(x(1),x(2));
    alpha22 = W22(x(1),x(2));
    alpha11 = W11(x(1),x(2));
    
    dxdt = alpha1*(A1*x+b1) + alpha2*(A2*x+b2) + alpha12*curve12(x) + alpha21*curve12(x) + ....
        alpha22*curve22(x) + alpha11*curve11(x);
end

function dXdt = curve12(x)
    r = sqrt(x(1)^2 + x(2)^2);
    theta = atan2(x(2),x(1));
    drdt = 3*r*(1-r);
    dthetadt= -1;
    dxdt = -r*sin(theta)*dthetadt+cos(theta)*drdt;
    dydt = r*cos(theta)*dthetadt + sin(theta)*drdt;
    dXdt = [dxdt; dydt];  
end

function dXdt = curve22(x)
    cloc = [2; 2];
    r = sqrt((x(1)-cloc(1))^2 + (x(2)-cloc(2))^2);
    theta = atan2(x(2)-cloc(2),x(1)-cloc(1));
    drdt = 3*r*(1-r);
    dthetadt= -1;
    dxdt = -r*sin(theta)*dthetadt+cos(theta)*drdt;
    dydt = r*cos(theta)*dthetadt + sin(theta)*drdt;
    dXdt = [dxdt; dydt];  
end

function dXdt = curve11(x)
    cloc = [-2; -2];
    r = sqrt((x(1)-cloc(1))^2 + (x(2)-cloc(2))^2);
    theta = atan2(x(2)-cloc(2),x(1)-cloc(1));
    drdt = 3*r*(1-r);
    dthetadt= -1;
    dxdt = -r*sin(theta)*dthetadt+cos(theta)*drdt;
    dydt = r*cos(theta)*dthetadt + sin(theta)*drdt;
    dXdt = [dxdt; dydt];  
end
