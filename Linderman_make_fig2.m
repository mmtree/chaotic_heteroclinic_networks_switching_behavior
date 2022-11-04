%% create a grid of fps for C. elegans, use a little bit of noise (SDE)

clear all; close all; clc;

% determine leave time
slow = 1/5;
med = 1/4;
fast = 1/2;

A1 = [-1 0; % blue: slow
    0 slow];

A2 = [-1 0; %purple: fast
    0 fast];

A2p2 = [slow 0; %purple2, crimson: slow
    0 -1];

A3 = [-1 0; %green: fast
    0 fast];

A5 = [fast 0; %red: use to be slow, now fast
    0 -1];

A7 = [med 0; %orange: use to be fast, now med
    0 -1];

A8 = [-1 0; %brown: was fast, now slow
    0 slow];

fp1_loc = [2; -1];    % blue
fp2_loc = [4; 1];     % purple
fp2p2_loc = [-11; 0];  % purple2
fp3_loc = [-5; -7];   % green
fp5_loc = [-9; -3];   % red
fp7_loc = [-5; -5];   % orange
fp8_loc = [-5; -3];   % brown

b1 = -A1*fp1_loc;
b2 = -A2*fp2_loc;
b2p2 = -A2p2*fp2p2_loc;
b3 = -A3*fp3_loc;
b5 = -A5*fp5_loc;
b7 = -A7*fp7_loc;
b8 = -A8*fp8_loc;

%%

x0_vec = 22*(rand(2,800)-[0.5; 0.5])+[-3;-3];
dt = 0.05;
tspan = 0:dt:10;

figure(); hold on;

for i = 1:length(x0_vec)
    x0 = x0_vec(:,i);
    x = NaN(length(tspan),2);
    x(1,:) = x0;
    for j = 2:length(tspan)
        dxdt = F_perturb(dt,x(j-1,:)',A1,A2,A2p2,A3,A5,A7,A8,b1,b2,b2p2,b3,b5,b7,b8);
        x(j,:) = x(j-1,:) + dt*dxdt';
    end
    %plot(x(1,1),x(1,2),'color','k');
    %plot(x(:,1),x(:,2),'.','color','c');
    plot(x(:,1),x(:,2),'color',[0.5 0.5 0.5]);
end

xlim([-14 6]);
ylim([-10 4]);

% plot single trajectory
tspan = 0:dt:1000;  %10000  %5000 (most recent)
x0 = [2; 0];
x = NaN(length(tspan),2);
x(1,:) = x0;
for j = 2:length(tspan)
    dxdt = F_perturb(dt,x(j-1,:)',A1,A2,A2p2,A3,A5,A7,A8,b1,b2,b2p2,b3,b5,b7,b8);
    x(j,:) = x(j-1,:) + dt*dxdt';
end

%%%%%%%%%%
plot(x(:,1),x(:,2),'.','color','k');

% plot fp locations with color
scatter(fp1_loc(1),fp1_loc(2),50,'b','linewidth',4); %fp1: blue
scatter(fp2_loc(1),fp2_loc(2),50,[0.4941 0.1843 0.5569],'linewidth',4); %fp2: purple
scatter(fp2p2_loc(1),fp2p2_loc(2),50,[0.4941 0.1843 0.5569],'linewidth',4); %fp2: purple
scatter(fp3_loc(1),fp3_loc(2),50,'g','linewidth',4); %fp3: green
scatter(fp5_loc(1),fp5_loc(2),50,'r','linewidth',4); %fp5: red
scatter(fp7_loc(1),fp7_loc(2),50,[ 0.9100 0.4100 0.1700],'linewidth',4); %fp7: orange
scatter(fp8_loc(1),fp8_loc(2),50,[0.2 0 0],'linewidth',4); %fp8: brown
scatter(-11,-1.5,50,[0.6350, 0.0780, 0.1840],'linewidth',4); %fp6: dark red
axis off

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
%print(fig,'figures/Linderman_flowlines','-dpdf','-r800')

%{
figure();
plot(x(:,1),'color','k');
%}

%% compute movement between fixed points

fp_vec = [];
for j = 1:length(x)
    d = Dist(x(j,:)');
    
    if d(1)< 0.6
        fp_vec = [fp_vec 1];
    elseif d(2) < 0.6
        fp_vec = [fp_vec 2];
    elseif d(3) < 0.6
        fp_vec = [fp_vec 3];
    elseif d(4) < 0.2
        fp_vec = [fp_vec 4];
    elseif d(5) < 0.6
        fp_vec = [fp_vec 5];
    elseif d(6) < 0.6  % 0.6 crimson
        fp_vec = [fp_vec 6];
    elseif d(7) < 0.6
        fp_vec = [fp_vec 7];
    elseif d(8) < 0.6
        fp_vec = [fp_vec 8];
    elseif d(9) < 0.6   % 0.2 % other purple
        fp_vec = [fp_vec 9];
    else
        fp_vec = [fp_vec 0];
    end
    if mod(j,1000)==0
        disp(j)
    end
end

fp_vec(fp_vec==0) = [];


idx9 = find(fp_vec==9);
fp_vec(idx9) = 2*ones(1,length(idx9));

fp_ts = fp_vec;

%{
%% compute leave times

fp_ts = fp_vec;

fp_vec = fp_vec(logical([1 diff(fp_vec,1,2)~=0]));

test= find(fp_ts==1);

diff(find((diff(test)~=1)==1))

for j = [1 2 3 4 5 6 7 8]  
    
end
%}

%%% compute Markov matrix

fp_vec = fp_vec(logical([1 diff(fp_vec,1,2)~=0]));
fp_vec = [fp_vec 0];
%{
figure(); hold on;
plot(fp_vec);
plot(fp_vec,'*');
%}

idx1 = find(fp_vec==1);
fp1_after = fp_vec(idx1+1);
fp1_after(fp1_after==0) = [];

idx2 = find(fp_vec==2);
fp2_after = fp_vec(idx2+1);
fp2_after(fp2_after==0) = [];

idx3 = find(fp_vec==3);
fp3_after = fp_vec(idx3+1);
fp3_after(fp3_after==0) = [];

idx4 = find(fp_vec==4);
fp4_after = fp_vec(idx4+1);
fp4_after(fp4_after==0) = [];

idx5 = find(fp_vec==5);
fp5_after = fp_vec(idx5+1);
fp5_after(fp5_after==0) = [];

idx6 = find(fp_vec==6);
fp6_after = fp_vec(idx6+1);
fp6_after(fp6_after==0) = [];

idx7 = find(fp_vec==7);
fp7_after = fp_vec(idx7+1);
fp7_after(fp7_after==0) = [];

idx8 = find(fp_vec==8);
fp8_after = fp_vec(idx8+1);
fp8_after(fp8_after==0) = [];

% idx22 = find(fp_vec==9);
% fp22_after = fp_vec(idx22+1);
% fp22_after(fp22_after==0) = [];

%%
Markov_mat = NaN(8,8);

ct1 = length(fp1_after);
ct2 = length(fp2_after);
ct3 = length(fp3_after);
ct4 = length(fp4_after);
ct5 = length(fp5_after);
ct6 = length(fp6_after);
ct7 = length(fp7_after);
ct8 = length(fp8_after);

%ct22 = length(fp22_after);

Markov_mat(1,:) = [sum(fp1_after==1)/ct1 sum(fp1_after==2)/ct1 sum(fp1_after==3)/ct1 sum(fp1_after==4)/ct1 sum(fp1_after==5)/ct1 sum(fp1_after==6)/ct1 sum(fp1_after==7)/ct1 sum(fp1_after==8)/ct1];
Markov_mat(2,:) = [sum(fp2_after==1)/ct2 sum(fp2_after==2)/ct2 sum(fp2_after==3)/ct2 sum(fp2_after==4)/ct2 sum(fp2_after==5)/ct2 sum(fp2_after==6)/ct2 sum(fp2_after==7)/ct2 sum(fp2_after==8)/ct2];
%Markov_mat(2,:) = [sum(fp2_after==1)/ct2 (sum(fp2_after==2)+sum(fp2_after==9))/ct2 sum(fp2_after==3)/ct2 sum(fp2_after==4)/ct2 sum(fp2_after==5)/ct2 sum(fp2_after==6)/ct2 sum(fp2_after==7)/ct2 sum(fp2_after==8)/ct2];

Markov_mat(3,:) = [sum(fp3_after==1)/ct3 sum(fp3_after==2)/ct3 sum(fp3_after==3)/ct3 sum(fp3_after==4)/ct3 sum(fp3_after==5)/ct3 sum(fp3_after==6)/ct3 sum(fp3_after==7)/ct3 sum(fp3_after==8)/ct3];
Markov_mat(4,:) = [sum(fp4_after==1)/ct4 sum(fp4_after==2)/ct4 sum(fp4_after==3)/ct4 sum(fp4_after==4)/ct4 sum(fp4_after==5)/ct4 sum(fp4_after==6)/ct4 sum(fp4_after==7)/ct4 sum(fp4_after==8)/ct4];
Markov_mat(5,:) = [sum(fp5_after==1)/ct5 sum(fp5_after==2)/ct5 sum(fp5_after==3)/ct5 sum(fp5_after==4)/ct5 sum(fp5_after==5)/ct5 sum(fp5_after==6)/ct5 sum(fp5_after==7)/ct5 sum(fp5_after==8)/ct5];
Markov_mat(6,:) = [sum(fp6_after==1)/ct6 sum(fp6_after==2)/ct6 sum(fp6_after==3)/ct6 sum(fp6_after==4)/ct6 sum(fp6_after==5)/ct6 sum(fp6_after==6)/ct6 sum(fp6_after==7)/ct6 sum(fp6_after==8)/ct6];
Markov_mat(7,:) = [sum(fp7_after==1)/ct7 sum(fp7_after==2)/ct7 sum(fp7_after==3)/ct7 sum(fp7_after==4)/ct7 sum(fp7_after==5)/ct7 sum(fp7_after==6)/ct7 sum(fp7_after==7)/ct7 sum(fp7_after==8)/ct7];
Markov_mat(8,:) = [sum(fp8_after==1)/ct8 sum(fp8_after==2)/ct8 sum(fp8_after==3)/ct8 sum(fp8_after==4)/ct8 sum(fp8_after==5)/ct8 sum(fp8_after==6)/ct8 sum(fp8_after==7)/ct8 sum(fp8_after==8)/ct8];

%% plot Markov matrix
%fig = figure('position', [0, 0, 600, 400]); hold on;
fig = figure();

xlab = [1 2 3 4 5 6 7 8];
ylab = [1 2 3 4 5 6 7 8];

colormap gray

cmap = colormap;
cmap = flip(cmap);
colormap(cmap)

imagesc(xlab,ylab,Markov_mat);
title('Markov transition matrix')

pbaspect([1 1 1]);
set(gca,'xaxisLocation','top')
xticks([1 2 3 4 5 6 7 8]);
yticks([1 2 3 4 5 6 7 8]);
colorbar
caxis([0 1]);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
%print(fig,'figures/Linderman_markov','-dpdf','-r800')

%% compute state timeseries and plot


hex_color_array=["3878bf","816286","7db171","feb30c","e40001","8b0011","f57506","b46952"];

hex_color_array = reshape(hex_color_array,[],1);  % reshape to Nx1 color array
hex_color_array = convertStringsToChars(hex_color_array);  % convert string array to character array
rgb_color_array = hex2rgb(hex_color_array);

mymap = rgb_color_array;

%{
mymap = [0    0.4471    0.7412;
    0.4314    0.3529    0.6000;
    0.3608    0.7098    0.3569;
    0.9608    0.8275    0.0863;
    0.8510    0.0353    0.0353;
    0.4902    0.0392    0.1216;
    0.9098    0.5294         0;
    0.6000    0.3725    0.2275;
    0.4314    0.3529    0.6000];
%}
    

    
%span = 20000:30000;
%span = 8000:15000;
%span = 50000:57000;

span = 1:10000;
%span = 40000:50000;
%span = 50000:60000;
%span = 80000:90000;
xvec = dt*span;
yvec = 1;
fig = figure('position', [0, 0, 500, 80]); hold on;
image(xvec,yvec,fp_ts(span))
title('behavioral states timeseries')

%imagesc(fp_ts(span));
colormap(gca,mymap);
xlim([xvec(1) xvec(end)]);
axis off

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
%print(fig,'figures/Linderman_example_fp_ts1','-dpdf','-r800')
%print(fig,'figures/Linderman_example_fp_ts2','-dpdf','-r800')
%print(fig,'figures/Linderman_example_fp_ts3','-dpdf','-r800')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot wave function  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xtest = 0:0.001:1;
ytest=[];
bias = 0.0; %0.4 % -1 <-> +1
s1 = 10;
s2 = 10;
endlen = 0.2;
for j = xtest
    ytest = [ytest wave3(j,bias,s1,s2,endlen)];
end


fig = figure('position', [0, 0, 90, 120]); hold on;
plot(xtest,0*ytest,'linewidth',2,'color',[0.75 0.75 0.75]);
plot(xtest,ytest,'linewidth',2,'color',[0.5 0.5 0.5]);
scatter(xtest(1:15:end),ytest(1:15:end),60,'.','k');
set(gca,'xticklabel',{[]},'yticklabel',{[]})
axis off
title('g(x)')

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
%print(fig,'figures/Linderman_example_wave','-dpdf','-r800')



%% plot weights for functions

figure();
x = -14:0.1:5;
y = -10:0.1:5;
[X,Y] = meshgrid(x,y);
%Z = W712(X,Y);

Z = W1(X,Y) + W2(X,Y) + W22(X,Y) + W225(X,Y) + W2252(X,Y) + W12(X,Y) + W21(X,Y) + W24(X,Y) + W4(X,Y) + W45(X,Y) +....
    W452(X,Y) + W5(X,Y) + W56(X,Y) + W6(X,Y) + W8(X,Y) +  W81(X,Y) + W13(X,Y) + W7(X,Y) +....
    W3(X,Y) + W35(X,Y) + W352(X,Y) + W353(X,Y) + W71(X,Y) + W712(X,Y) + W78(X,Y) +....
    + W58(X,Y) + W812(X,Y);

%Z =  W1(X,Y) + W13(X,Y);


surf(X,Y,Z);
xlabel('x');
ylabel('y');
title('weighting functions')



%% functions

function dxdt = F(x,A1,A2,A2p2,A3,A5,A7,A8,b1,b2,b2p2,b3,b5,b7,b8)
    alpha1 = W1(x(1),x(2));
    alpha2 = W2(x(1),x(2));
    alpha22 = W22(x(1),x(2));
    alpha12 = W12(x(1),x(2));
    alpha21 = W21(x(1),x(2));
    alpha24 = W24(x(1),x(2));
    alpha225 = W225(x(1),x(2));
    alpha2252 = W2252(x(1),x(2));  
    alpha4 = W4(x(1),x(2));
    alpha45 = W45(x(1),x(2));
    alpha452 = W452(x(1),x(2));
    alpha5 = W5(x(1),x(2));
    alpha56 = W56(x(1),x(2));
    alpha6 = W6(x(1),x(2));
    alpha8 = W8(x(1),x(2));
    alpha81 = W81(x(1),x(2));
    alpha13 = W13(x(1),x(2));
    alpha7 = W7(x(1),x(2));
    alpha3 = W3(x(1),x(2));
    alpha35 = W35(x(1),x(2));
    alpha352 = W352(x(1),x(2));
    alpha353 = W353(x(1),x(2));
    alpha71 = W71(x(1),x(2));
    alpha712 = W712(x(1),x(2));
    alpha78 = W78(x(1),x(2));
    alpha58 = W58(x(1),x(2));
    alpha812 = W812(x(1),x(2));
    
    fp6_dyn = [-0.1*(x(1)+11);
                (x(2)+1.5)^2+0.1];
    
    dxdt = alpha1*(A1*x+b1) + alpha2*(A2*x+b2) + alpha22*(A2p2*x+b2p2) + alpha12*curve12(x) + alpha21*curve12(x) + ....
           alpha24*curve24(x) + alpha225*curve25(x) + alpha2252*trans2252(x) + alpha4*trans4(x) + alpha45*curve45(x) + ....
         alpha452*trans452(x) + alpha5*(A5*x+b5) + alpha56*curve56(x) + alpha6*fp6_dyn + alpha8*(A8*x+b8) +....
         alpha81*curve81(x) + alpha13*curve13(x) + alpha7*(A7*x+b7) + alpha3*(A3*x+b3) +....
         alpha35*curve35(x) + alpha352*trans352(x) + alpha353*curve353(x) + alpha71*curve71(x) +....
         alpha712*trans712(x) + alpha78*curve78(x) + alpha58*trans58(x) + alpha812*trans812(x);
 
end


function dxdt = F_perturb(dt,x,A1,A2,A2p2,A3,A5,A7,A8,b1,b2,b2p2,b3,b5,b7,b8)
    alpha1 = W1(x(1),x(2));
    alpha2 = W2(x(1),x(2));
    alpha22 = W22(x(1),x(2));
    alpha12 = W12(x(1),x(2));
    alpha21 = W21(x(1),x(2));
    alpha24 = W24(x(1),x(2));
    alpha225 = W225(x(1),x(2));
    alpha2252 = W2252(x(1),x(2));
    alpha4 = W4(x(1),x(2));
    alpha45 = W45(x(1),x(2));
    alpha452 = W452(x(1),x(2));
    alpha5 = W5(x(1),x(2));
    alpha56 = W56(x(1),x(2));
    alpha6 = W6(x(1),x(2));
    alpha8 = W8(x(1),x(2));
    alpha81 = W81(x(1),x(2));
    alpha13 = W13(x(1),x(2));
    alpha7 = W7(x(1),x(2));
    alpha3 = W3(x(1),x(2));
    alpha35 = W35(x(1),x(2));
    alpha352 = W352(x(1),x(2));
    alpha353 = W353(x(1),x(2));
    alpha71 = W71(x(1),x(2));
    alpha712 = W712(x(1),x(2));
    alpha78 = W78(x(1),x(2));
    alpha58 = W58(x(1),x(2));
    alpha812 = W812(x(1),x(2));
    
    fp6_dyn = [-0.1*(x(1)+11);
                (x(2)+1.5)^2+0.1];
    
    dxdt = alpha1*(A1*x+b1) + alpha2*(A2*x+b2) + alpha22*(A2p2*x+b2p2)  + alpha12*curve12(x) + alpha21*curve12(x) + ....
          + alpha24*curve24(x) + alpha225*curve25(x) + alpha2252*trans2252(x) + alpha4*trans4(x) + alpha45*curve45(x) + ....
         alpha452*trans452(x) + alpha5*(A5*x+b5) + alpha56*curve56(x) + alpha6*fp6_dyn + alpha8*(A8*x+b8) +....
         alpha81*curve81(x) + alpha13*curve13(x) + alpha7*(A7*x+b7) + alpha3*(A3*x+b3) +....
         alpha35*curve35(x) + alpha352*trans352(x) + alpha353*curve353(x) + alpha71*curve71(x) +....
         alpha712*trans712(x) + alpha78*curve78(x) + alpha58*trans58(x) + alpha812*trans812(x);
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%  Universal parameters  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    amp = 0.2; % 2
    bias = 0.0; %0.4 % -1 <-> +1
    s1 = 10; % use to be 20
    s2 = 10;
    endlen = 0.2; %use to be 0.1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% TO FP1 %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bias = 0.18; %0.4 % -1 <-> +1
    
    %%% sine kink fp2 to fp1 %%%
    xs1 = [3; -1]; % starting point for sine wave
    dxdt_tmp = F(xs1,A1,A2,A2p2,A3,A5,A7,A8,b1,b2,b2p2,b3,b5,b7,b8);
    xe1 = xs1 + dt*dxdt_tmp; % ending point for sine wave
    if (x(1)>min(xs1(1),xe1(1))) && (x(1)<max(xs1(1),xe1(1))) && (x(2)>min(xs1(2),xe1(2))-0.5) && (x(2)<max(xs1(2),xe1(2))+0.5)
        L = abs(xs1(1)-xe1(1));
        extra_x = 0;
        %extra_y = amp*(sin((1/L)*2*pi*(x(1)-xs1(1))+bias)-bias);
        extra_y = amp*(wave3((1/L)*(x(1)-xs1(1)),bias,s1,s2,endlen)); % use wave instead
        
        dxdt = dxdt+[extra_x; extra_y];
    end
    
    %%% sine kink fp8 to fp1 %%%
    xs1 = [1; -1]; % starting point for sine wave
    dxdt_tmp = F(xs1,A1,A2,A2p2,A3,A5,A7,A8,b1,b2,b2p2,b3,b5,b7,b8);
    xe1 = xs1 + dt*dxdt_tmp; % ending point for sine wave
    if (x(1)>min(xs1(1),xe1(1))) && (x(1)<max(xs1(1),xe1(1))) && (x(2)>min(xs1(2),xe1(2))-0.5) && (x(2)<max(xs1(2),xe1(2))+0.5)
        L = abs(xs1(1)-xe1(1));
        extra_x = 0;
        %extra_y = amp*(sin((1/L)*2*pi*(x(1)-xs1(1))+bias)-bias);
        extra_y = amp*(wave3((1/L)*(x(1)-xs1(1)),bias,s1,s2,endlen)); % use wave instead
        dxdt = dxdt+[extra_x; extra_y];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% TO FP2 %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bias = -0.27; % -1 <-> +1 %-0.18
    
    %%% sine kink fp1 to fp2 %%%
    xs1 = [3; 1]; % starting point for sine wave
    dxdt_tmp = F(xs1,A1,A2,A2p2,A3,A5,A7,A8,b1,b2,b2p2,b3,b5,b7,b8);
    xe1 = xs1 + dt*dxdt_tmp; % ending point for sine wave
    if (x(1)>min(xs1(1),xe1(1))) && (x(1)<max(xs1(1),xe1(1))) && (x(2)>min(xs1(2),xe1(2))-0.5) && (x(2)<max(xs1(2),xe1(2))+0.5)
        L = abs(xs1(1)-xe1(1));
        extra_x = 0;
        %extra_y = amp*(sin((1/L)*2*pi*(x(1)-xs1(1))+bias)-bias);
        extra_y = amp*(wave3((1/L)*(x(1)-xs1(1)),bias,s1,s2,endlen)); % use wave instead
        dxdt = dxdt+[extra_x; extra_y];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% TO FP2 point 2 %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    amp = 1; 
    bias = 0.08; % -1 <-> +1  % fix to make more crimson to red
    
    %%% sine kink fp6 to fp2p2 %%%
    xs1 = [-11; -1]; % starting point for sine wave
    dxdt_tmp = F(xs1,A1,A2,A2p2,A3,A5,A7,A8,b1,b2,b2p2,b3,b5,b7,b8);
    xe1 = xs1 + dt*dxdt_tmp; % ending point for sine wave
    if (x(1)>min(xs1(1),xe1(1))-0.5) && (x(1)<max(xs1(1),xe1(1))+0.5) && (x(2)>min(xs1(2),xe1(2))) && (x(2)<max(xs1(2),xe1(2)))
        L = abs(xs1(2)-xe1(2));
        %extra_x = amp*(sin((1/L)*2*pi*(x(2)-xs1(2))+bias)-bias);
        extra_x = amp*(wave3((1/L)*(x(2)-xs1(2)),bias,s1,s2,endlen)); % use wave instead
        extra_y = 0;
        dxdt = dxdt+[extra_x; extra_y];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% TO FP3 %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    amp = 0.5; % 2
    bias = 0.21; % -1 <-> +1
    
    %%% sine kink fp1 to fp3 %%%
    xs1 = [-3; -7]; % starting point for sine wave
    dxdt_tmp = F(xs1,A1,A2,A2p2,A3,A5,A7,A8,b1,b2,b2p2,b3,b5,b7,b8);
    xe1 = xs1 + dt*dxdt_tmp; % ending point for sine wave
    if (x(1)>min(xs1(1),xe1(1))) && (x(1)<max(xs1(1),xe1(1))) && (x(2)>min(xs1(2),xe1(2))-0.5) && (x(2)<max(xs1(2),xe1(2))+0.5)
        L = abs(xs1(1)-xe1(1));
        extra_x = 0;
        %extra_y = amp*(sin((1/L)*2*pi*(x(1)-xs1(1))+bias)-bias);
        extra_y = amp*(wave3((1/L)*(x(1)-xs1(1)),bias,s1,s2,endlen)); % use wave instead
        dxdt = dxdt+[extra_x; extra_y];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% TO FP5 %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    amp = 0.1;
    bias = 0.15; % -1 <-> +1
    
    %%% sine kink fp4 to fp5 %%% heteroclinic orbit from 4 to 5
    xs1 = [-9; -2]; % starting point for sine wave
    dxdt_tmp = F(xs1,A1,A2,A2p2,A3,A5,A7,A8,b1,b2,b2p2,b3,b5,b7,b8);
    xe1 = xs1 + dt*dxdt_tmp; % ending point for sine wave
    if (x(1)>min(xs1(1),xe1(1))-0.5) && (x(1)<max(xs1(1),xe1(1))+0.5) && (x(2)>min(xs1(2),xe1(2))) && (x(2)<max(xs1(2),xe1(2)))
        L = abs(xs1(2)-xe1(2));
        %extra_x = amp*(sin((1/L)*2*pi*(x(2)-xs1(2))+bias)-bias);
        extra_x = amp*(wave3((1/L)*(x(2)-xs1(2)),bias,s1,s2,endlen)); % use wave instead
        extra_y = 0;
        dxdt = dxdt+[extra_x; extra_y];
    end
    
    %%% sine kink fp3 to fp5 %%% heteroclinic orbit from 3 to 5
    xs1 = [-9; -8]; % starting point for sine wave
    dxdt_tmp = F(xs1,A1,A2,A2p2,A3,A5,A7,A8,b1,b2,b2p2,b3,b5,b7,b8);
    xe1 = xs1 + dt*dxdt_tmp; % ending point for sine wave
    if (x(1)>min(xs1(1),xe1(1))-0.5) && (x(1)<max(xs1(1),xe1(1))+0.5) && (x(2)>min(xs1(2),xe1(2))) && (x(2)<max(xs1(2),xe1(2)))
        L = abs(xs1(2)-xe1(2));
        %extra_x = amp*(sin((1/L)*2*pi*(x(2)-xs1(2))+bias)-bias);
        extra_x = amp*(wave3((1/L)*(x(2)-xs1(2)),bias,s1,s2,endlen)); % use wave instead
        extra_y = 0;
        dxdt = dxdt+[extra_x; extra_y];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% TO FP7 %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    amp = 1; 
    bias = -0.2; % -1 <-> +1
    
    %%% sine kink fp3 to fp7 %%% 
    xs1 = [-5; -6]; % starting point for sine wave
    dxdt_tmp = F(xs1,A1,A2,A2p2,A3,A5,A7,A8,b1,b2,b2p2,b3,b5,b7,b8);
    xe1 = xs1 + dt*dxdt_tmp; % ending point for sine wave
    if (x(1)>min(xs1(1),xe1(1))-0.5) && (x(1)<max(xs1(1),xe1(1))+0.5) && (x(2)>min(xs1(2),xe1(2))) && (x(2)<max(xs1(2),xe1(2)))
        L = abs(xs1(2)-xe1(2));
        %extra_x = amp*(sin((1/L)*2*pi*(x(2)-xs1(2))+bias)-bias);
        extra_x = amp*(wave3((1/L)*(x(2)-xs1(2)),bias,s1,s2,endlen)); % use wave instead
        extra_y = 0;
        dxdt = dxdt+[extra_x; extra_y];
    end
    
    %%% sine kink fp8 to fp7 %%% 
    xs1 = [-5; -4]; % starting point for sine wave
    dxdt_tmp = F(xs1,A1,A2,A2p2,A3,A5,A7,A8,b1,b2,b2p2,b3,b5,b7,b8);
    xe1 = xs1 + dt*dxdt_tmp; % ending point for sine wave
    if (x(1)>min(xs1(1),xe1(1))-0.5) && (x(1)<max(xs1(1),xe1(1))+0.5) && (x(2)>min(xs1(2),xe1(2))) && (x(2)<max(xs1(2),xe1(2)))
        L = abs(xs1(2)-xe1(2));
        %extra_x = amp*(sin((1/L)*2*pi*(x(2)-xs1(2))+bias)-bias);
        extra_x = amp*(wave3((1/L)*(x(2)-xs1(2)),bias,s1,s2,endlen)); % use wave instead
        extra_y = 0;
        dxdt = dxdt+[extra_x; extra_y];
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% TO FP8 %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    amp = 0.2; 
    bias = 0.2; % -1 <-> +1
    
    %%% sine kink fp5 to fp8 %%%
    xs1 = [-6; -3]; % starting point for sine wave
    dxdt_tmp = F(xs1,A1,A2,A2p2,A3,A5,A7,A8,b1,b2,b2p2,b3,b5,b7,b8);
    xe1 = xs1 + dt*dxdt_tmp; % ending point for sine wave
    if (x(1)>min(xs1(1),xe1(1))) && (x(1)<max(xs1(1),xe1(1))) && (x(2)>min(xs1(2),xe1(2))-0.5) && (x(2)<max(xs1(2),xe1(2))+0.5)
        L = abs(xs1(1)-xe1(1));
        extra_x = 0;
        %extra_y = amp*(sin((1/L)*2*pi*(x(1)-xs1(1))+bias)-bias);
        extra_y = amp*(wave3((1/L)*(x(1)-xs1(1)),bias,s1,s2,endlen)); % use wave instead
        dxdt = dxdt+[extra_x; extra_y];
    end
end


%% Distance
function d = Dist(x)
    fp1 = [2; -1]; % blue
    
    fp2 = [4; 1];   % purple
    
    fp3 = [-5; -7];   % green
    
    %fp4 = [-1; 3];    % yellow
  
    fp5 = [-9; -3];   % red
    
    fp6 = [-11; 0];  % dark red %[-11; -1.5];
    fp7 = [-5; -5];   % orange
    fp8 = [-5; -3];   % brown
    
    fp22 = [-13; 1];  % purple 2  %[-12; 0];

    dfp1 = norm(x-fp1);
    dfp2 = norm(x-fp2);
    
    %dfp3 = norm(x-fp3);
    
    % green region, L-shape
    if ((x(1)>2-0.5 && x(1)<2+0.5) && (x(2)>-7-0.5 && x(2) < -3)) | ((x(1)>-5 && x(1)<3) && (x(2)>-7-0.5 && x(2) < -7+0.5))
        dfp3=0;
    else
        dfp3=2;
    end
   
    
    %dfp4 = norm(x-fp4);
    dfp4 = abs(x(2)-3);  % yellow region
    if dfp4<0.6 % yellow
    elseif dfp4>0.6 && ((x(2)>0.2 && x(2)<3.5) && (x(1)>-9-0.5 && x(1) < -9+0.5))  % add a curve to yellow region
        dfp4 = 0;
    else
        dfp4 = 2;
    end

    
    dfp5 = norm(x-fp5);
    dfp6 = norm(x-fp6);
    
    
    dfp7 = norm(x-fp7); % orange
    if dfp7<0.6 % orange
    elseif dfp7>0.6 && (((x(2)>-5-0.5 && x(2)<-5+0.5) && (x(1)>-5 && x(1) < 0)) | ((x(1)>-1-0.5 && x(1)<-1+0.5) && (x(2)>-5-0.5 && x(2) < -2)))
        dfp7=0;
    else
        dfp7=2;
    end
    
    
    dfp8 = norm(x-fp8); % tan
    
    if dfp8<0.6 % tan
        
    elseif dfp8>0.6 && (((x(1)>-5-0.5 && x(1)<-5+0.5) && (x(2)>-5 && x(2) < 0)) | ((x(2)>-1-0.5 && x(2)<-1+0.5) && (x(1)>-5-0.5 && x(1) < -2)))
        dfp8=0;
    else
        dfp8=2;
    end
    
    dfp22 = norm(x-fp22);
    
    d = [dfp1; dfp2; dfp3; dfp4; dfp5; dfp6; dfp7; dfp8; dfp22];
end

%% Weighting functions

% Weighting functions
function w1 = W1(x,y)  % weighting fp1
    s = 10;
    w1 = -(1/4)*(tanh(s*y)-tanh(s*(y+6))).*(-tanh(s*(x-3))-tanh(-s*(x-1)));
end

function w2 = W2(x,y) % weighting fp2
    s = 10;
    w2 = (1/4)*(tanh(s*y)-tanh(s*(y-2))).*(tanh(s*(x-3))-tanh(s*(x-5)));
end

function w22 = W22(x,y) % weighting fp2
    s = 10;
    w22 = (1/4)*(tanh(s*(y+1))-tanh(s*(y-1))).*(tanh(s*(x+12))-tanh(s*(x+10)));
end

function w225 = W225(x,y) % weighting fp2
    s = 10;
    w225 = (1/4)*(tanh(s*(y-3))-tanh(s*(y+1))).*(tanh(s*(x+12))-tanh(s*(x+14)));
end

function w2252 = W2252(x,y) % weighting fp2
    s = 10;
    w2252 = (1/4)*(tanh(s*(y-3))-tanh(s*(y-1))).*(tanh(s*(x+10))-tanh(s*(x+12)));
end

function w12 = W12(x,y) % weighting for fp1 -> fp2
    s = 10;
    w12 = -(1/4)*(tanh(s*y)-tanh(s*(y-2))).*(tanh(s*(x-3))-tanh(s*(x-1)));
end

function w21 = W21(x,y) % weighting for fp2 -> fp1
    s = 10;
    w21 = -(1/4)*(tanh(s*y)-tanh(s*(y+2))).*(tanh(s*(x-3))-tanh(s*(x-5)));
end

function w24 = W24(x,y) % weighting for fp2 -> fp4 (not really fp, yellow)
    s = 10;
    w24 = (1/4)*(tanh(s*(y-2))-tanh(s*(y-4))).*(tanh(s*(x-3))-tanh(s*(x-5)));
end

function w4 = W4(x,y) % weighting for fp2 -> fp4 (not really fp, yellow)
    s = 10;
    w4 = -(1/4)*(tanh(s*(y-2))-tanh(s*(y-4))).*(tanh(s*(x-3))-tanh(s*(x+8)));
end

function w45 = W45(x,y) % weighting for fp2 -> fp4 (not really fp, yellow)
    s = 10;
    w45 = -(1/4)*(tanh(s*(y-2))-tanh(s*(y-4))).*(tanh(s*(x+8))-tanh(s*(x+10)));
end

function w452 = W452(x,y) % weighting for fp2 -> fp4 (not really fp, yellow)
    s = 10;
    w452 = (1/4)*(tanh(s*(y+2))-tanh(s*(y-2))).*(tanh(s*(x+10))-tanh(s*(x+8)));
end

function w5 = W5(x,y)  % weighting fp1
    s = 10;
    w5 = -(1/4)*(tanh(s*(y+8))-tanh(s*(y+2))).*(-tanh(s*(x+10))-tanh(-s*(x+8)));
end

function w56 = W56(x,y)  % weighting fp1
    s = 10;
    w56 = (1/4)*(tanh(s*(y+4))-tanh(s*(y+2))).*(-tanh(s*(x+10))-tanh(-s*(x+12)));
end

function w6 = W6(x,y)  % weighting fp1
    s = 10;
    w6 = (1/4)*(tanh(s*(y+2))-tanh(s*(y+1))).*(-tanh(s*(x+10))-tanh(-s*(x+12)));
end

function w8 = W8(x,y)  % weighting fp1
    s = 10;
    w8 = -(1/4)*(tanh(s*(y+4))-tanh(s*(y+2))).*(-tanh(s*(x+6))-tanh(-s*(x+4)));
end

function w81 = W81(x,y)  % weighting fp1
    s = 10;
    w81 = -(1/4)*(tanh(s*(y+2))-tanh(s*(y))).*(-tanh(s*(x+6))-tanh(-s*(x+4)));
end

function w13 = W13(x,y)  % weighting fp1
    s = 10;
    w13 = -(1/4)*(tanh(s*(y+8))-tanh(s*(y+6))).*(-tanh(s*(x-1))-tanh(-s*(x-3)));
end

function w7 = W7(x,y)  % weighting fp1
    s = 10;
    w7 = -(1/4)*(tanh(s*(y+6))-tanh(s*(y+4))).*(-tanh(s*(x+6))-tanh(-s*(x+2)));
end

function w3 = W3(x,y)  % weighting fp1
    s = 10;
    w3 = -(1/4)*(tanh(s*(y+8))-tanh(s*(y+6))).*(-tanh(s*(x+6))-tanh(-s*(x-1)));
end

function w35 = W35(x,y)  % weighting fp1
    s = 10;
    w35 = -(1/4)*(tanh(s*(y+10))-tanh(s*(y+8))).*(-tanh(s*(x+6))-tanh(-s*(x+4)));
end

function w352 = W352(x,y)  % weighting fp1
    s = 10;
    w352 = -(1/4)*(tanh(s*(y+10))-tanh(s*(y+8))).*(-tanh(s*(x+8))-tanh(-s*(x+6)));
end

function w353 = W353(x,y)  % weighting fp1
    s = 10;
    w353 = -(1/4)*(tanh(s*(y+10))-tanh(s*(y+8))).*(-tanh(s*(x+10))-tanh(-s*(x+8)));
end

function w71 = W71(x,y)  % weighting fp1
    s = 10;
    w71 = -(1/4)*(tanh(s*(y+6))-tanh(s*(y+4))).*(-tanh(s*(x+2))-tanh(-s*(x-0)));
end

function w712 = W712(x,y)  % weighting fp1
    s = 10;
    w712 = -(1/4)*(tanh(s*(y+4))-tanh(s*(y+2))).*(-tanh(s*(x+4))-tanh(-s*(x-1)));
end

function w78 = W78(x,y)  % weighting fp1
    s = 10;
    w78 = -(1/4)*(tanh(s*(y+6))-tanh(s*(y+4))).*(-tanh(s*(x+8))-tanh(-s*(x+6)));
end

function w58 = W58(x,y)  % weighting fp1
    s = 10;
    w58 = -(1/4)*(tanh(s*(y+4))-tanh(s*(y+2))).*(-tanh(s*(x+8))-tanh(-s*(x+6)));
end

function w812 = W812(x,y)  % weighting fp1
    s = 10;
    w812 = -(1/4)*(tanh(s*(y+2))-tanh(s*(y+0))).*(-tanh(s*(x+4))-tanh(-s*(x-1)));
end

%% Connecting functions (heteroclinic orbits)

function dXdt = curve12(x)
    cloc = [3; 0];
    dXdt = polar_cw(x,cloc);
end

function dXdt = curve24(x)
    cloc = [3; 2];
    dXdt = polar_ccw(x,cloc);
end

function dXdt = curve25(x) %%%% slow down purple turn???
    cloc = [-12; 1];
    dXdt = 0.1*polar_cw(x,cloc);
end

function dXdt = trans2252(x)
    dxdt = 1;
    dydt = -3*(x(2)-2);
    dXdt = [dxdt; dydt];  
end

function dXdt = trans4(x)  %yellow
    dxdt = -0.25;  % use to be -0.8, then -0.5
    dydt = -3*(x(2)-3);
    dXdt = [dxdt; dydt];  
end

function dXdt = curve45(x)
    cloc = [-8; 2];
    dXdt = polar_ccw(x,cloc);
end

function dXdt = trans452(x)
    dxdt = -6*(x(1)+9);
    dydt = -0.2;
    dXdt = [dxdt; dydt];  
end

function dXdt = curve56(x)
    cloc = [-10; -2];
    dXdt = polar_cw(x,cloc); 
end

% function dXdt = curve11(x)
%     cloc = [-2; -2];
%     dXdt = polar_cw(x,cloc); 
% end

function dXdt = curve81(x)
    cloc = [-4; -2];
    dXdt = polar_cw(x,cloc);
end

function dXdt = curve13(x)
    cloc = [1; -6];
    dXdt = polar_cw(x,cloc);
end

function dXdt = curve35(x)
    cloc = [-6; -8];
    dXdt = polar_cw(x,cloc);
end

function dXdt = trans352(x)
    dxdt = -1;
    dydt = -3*(x(2)+9);
    dXdt = [dxdt; dydt];  
end

function dXdt = curve353(x)
    cloc = [-8; -8];
    dXdt = polar_cw(x,cloc);
end

function dXdt = curve71(x)
    cloc = [-2; -4];
    dXdt = polar_ccw(x,cloc);
end

function dXdt = trans712(x)
    dxdt = -3*(x(1)+1);
    dydt = 0.5;  % use to be 1, slow down
    dXdt = [dxdt; dydt];  
end

function dXdt = curve78(x)
    cloc = [-6; -4];
    dXdt = polar_cw(x,cloc);
end

function dXdt = trans58(x)
    dxdt = 1;
    dydt = -6*(x(2)+3);
    dXdt = [dxdt; dydt];  
end

function dXdt = trans812(x)
    dxdt = 0.3;  % use to be 1, then 0.5
    dydt = -3*(x(2)+1);
    dXdt = [dxdt; dydt];  
end

function dXdt = polar_cw(x,cloc)
    r = sqrt((x(1)-cloc(1))^2 + (x(2)-cloc(2))^2);
    theta = atan2(x(2)-cloc(2),x(1)-cloc(1));
    drdt = 3*r*(1-r);
    dthetadt= -1;
    dxdt = -r*sin(theta)*dthetadt+cos(theta)*drdt;
    dydt = r*cos(theta)*dthetadt + sin(theta)*drdt;
    dXdt = [dxdt; dydt];  
end

function dXdt = polar_ccw(x,cloc)
    r = sqrt((x(1)-cloc(1))^2 + (x(2)-cloc(2))^2);
    theta = atan2(x(2)-cloc(2),x(1)-cloc(1));
    drdt = 3*r*(1-r);
    dthetadt= 1;
    dxdt = -r*sin(theta)*dthetadt+cos(theta)*drdt;
    dydt = r*cos(theta)*dthetadt + sin(theta)*drdt;
    dXdt = [dxdt; dydt];  
end

%%%%%%%%%%%%%%%%%%%%%%%
%%%%% wave function %%%
%%%%%%%%%%%%%%%%%%%%%%%

function wave = wave3(x,b,s1,s2,endlen)
    % b 0->1/2, bias term ranges from 0 to 0.5
    %s1 = 20;
    %s2 = 50;
    %endlen = 0.25;
    wave = 1/2*(tanh(s2*(x-endlen))+1)-tanh(s2*(x-(0.5+b)))-1+1/2*(tanh(s1*(x-(1-endlen)))+1)+1/2*(tanh(s1*(x-(1+endlen)))+1)+....
1/2*(tanh(s1*(x+endlen))+1)-tanh(s2*(x+(0.5-b)))-1+1/2*(tanh(s2*(x+(1-endlen)))+1)+1/2*(tanh(s2*(x+(1+endlen)))-1);
end
