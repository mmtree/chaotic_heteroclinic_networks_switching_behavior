%% create a grid of fps for C. elegans: forward, reversal, turn quiescence
clear all; close all; clc;
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',11)

%% specify specialty colors

pinkr = [0.9490    0.4000    0.4000];
purpleq = [0.2549    0.2667    0.8784];

%% Ideal matrix want to fit to (Target)

ideal_mat = [0 1/4 3/4 0;
    2/3 0 1/3 0;
    0 0.05 0 0.95;
    1 0 0 0];

fig = figure('position', [0, 0, 230, 200]); hold on;
xlab = [1 2 3 4];
ylab = [1 2 3 4];
colormap gray
cmap = colormap;
cmap = flip(cmap);
colormap(cmap)

imagesc(xlab,ylab,ideal_mat);

%pbaspect([1 1 1]);
set(gca,'xaxisLocation','top')
xticks([1 2 3 4]);
yticks([1 2 3 4]);
colorbar
%caxis([0 1]);
%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/Nich_Markov_mat_ideal','-dpdf','-r800')
%}

%%

% determine leave time
slow = 1/10;
med = 1/2;
fast = 1;

A1 = [med 0; % light blue: forward
    0 -3];

A2 = [med 0; % dark blue: quiescence
    0 -3];

A3 = [med 0; % yellow: reversal
    0 -3];


fp1_loc = [-3; -3];  % light blue
fp2_loc = [0; 0];    % dark blue
fp3_loc = [3; 3];    % yellow

b1 = -A1*fp1_loc;
b2 = -A2*fp2_loc;
b3 = -A3*fp3_loc;

%% plot flow lines

x0_vec = 14*(rand(2,400)-[0.5; 0.5])+[0;0];  % 2000 points
dt = 0.01;  % dt = 0.05 is best step size for getting heteroclinic curve correct
tspan = 0:dt:8;

fig = figure('position', [0, 0, 300, 260]); hold on;

for i = 1:length(x0_vec)
    x0 = x0_vec(:,i);
    x = NaN(length(tspan),2);
    x(1,:) = x0;
    for j = 2:length(tspan)
        dxdt = F_perturb(dt,x(j-1,:)',A1,A2,A3,b1,b2,b3);
        x(j,:) = x(j-1,:) + dt*dxdt';
    end
    plot(x(1,1),x(1,2),'color','c');
    %plot(x(:,1),x(:,2),'.','color','c');
    %plot(x(:,1),x(:,2),'.','color',[0.5 0.5 0.5]);
    plot(x(:,1),x(:,2),'color',[0.5 0.5 0.5]);
end

xlim([-7 7]);
ylim([-7 7]);


% plot single trajectory
tspan = 0:dt:8000;  % 8000
x0 = [-1; 0];
x = NaN(length(tspan),2);
x(1,:) = x0;
for j = 2:length(tspan)
    dxdt = F_perturb(dt,x(j-1,:)',A1,A2,A3,b1,b2,b3);
    x(j,:) = x(j-1,:) + dt*dxdt';
end

%%%%%%%%%%
plot(x(:,1),x(:,2),'.','color','k','MarkerSize',6);

% plot fp locations with color
scatter(fp1_loc(1),fp1_loc(2),30,'c','linewidth',1); %fp1: light blue, forward
scatter(fp2_loc(1),fp2_loc(2),30,'b','linewidth',1); %fp2: dark blue, quiescent
scatter(fp3_loc(1),fp3_loc(2),30,[0.9290, 0.6940, 0.1250],'linewidth',1); %fp3: yellow, reversal
axis off
%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/Nich_flowlines','-dpdf','-r800')
%}

%% compute movement between fixed points

rad = 0.6;
fp_vec = [];
for j = 1:length(x)
    d = Dist(x(j,:)');
    
    if d(1)< rad
        fp_vec = [fp_vec 1];
    elseif d(2) < rad
        fp_vec = [fp_vec 2];
    elseif d(3) < rad
        fp_vec = [fp_vec 3];
    elseif d(4) < rad   % use to be rad
        fp_vec = [fp_vec 4];
    else
        fp_vec = [fp_vec 0];
    end
    if mod(j,100000)==0
        disp(j)
    end
end

fp_vec_w_nostate = fp_vec;

fp_vec(fp_vec==0) = [];


%% compute leave times

fp_ts = fp_vec;
dwell_times = cell(4,1);  % ones cell for each state

st_val = fp_ts(1);
ct = 0;
for j = 1:length(fp_ts)
    if fp_ts(j) ==st_val
        ct = ct+1;
    else
        dwell_times{st_val} = [dwell_times{st_val} ct];
        st_val = fp_ts(j);
        ct=1; 
    end  
end

%% compute state timeseries and plot

mymap = [0 1 1;
    	0 0 1;
        0.9290, 0.6940, 0.1250;
        1 0 1;];
    
span = 5000:20800;
xvec = dt*span;
yvec = 1;
fig = figure('position', [0, 0, 400, 40]); hold on;
image(xvec,yvec,fp_ts(span))

%imagesc(fp_ts(span));
colormap(gca,mymap);
xlim([xvec(1) xvec(end)]);
axis off
%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/Nich_example_fp_ts','-dpdf','-r800')
%}

%% compute state timeseries and plot with nostate

fp_vec_w_nostate(fp_vec_w_nostate==0)=5;

mymap = [0 1 1;
    	0 0 1;
        0.9290, 0.6940, 0.1250;
        1 0 1
        0.5 0.5 0.5];
    
span = 5800:26000;
xvec = dt*span;
yvec = 1;
fig = figure('position', [0, 0, 400, 40]); hold on;
image(xvec,yvec,fp_vec_w_nostate(span))

%imagesc(fp_ts(span));
colormap(gca,mymap);
xlim([xvec(1) xvec(end)]);
axis off
%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/Nich_example_fp_ts_grey','-dpdf','-r800')
%}
%% compute next state after forward, bin by forward time

dwell_t_next_state = cell(2,1);  % cell 1: to state 2, cell 2: to state 3
fp_past = fp_ts(1);
ct = 0;
for j = 1:length(fp_ts)
    if fp_ts(j)==1
        ct = ct+1;
        fp_past = 1;
    elseif fp_ts(j)==2 && fp_past==1
        dwell_t_next_state{1} = [dwell_t_next_state{1} ct];  % add count to fp after = 2
        ct = 0;
        fp_past = 2;
    elseif fp_ts(j)==3 && fp_past==1
        dwell_t_next_state{2} = [dwell_t_next_state{2} ct];  % add count to fp after = 3
        ct = 0;
        fp_past = 3;
    end       
end

% compute bins for histrogram 
dwell_t_next_state2 = dwell_t_next_state{1};
dwell_t_next_state3 = dwell_t_next_state{2};

num_to2 = length(dwell_t_next_state2);
num_to3 = length(dwell_t_next_state3);

to2_time_bins = zeros(3,1); % three time bins
to3_time_bins = zeros(3,1); % three time bins

val = 1.3;   %1.3
cut1 = val/dt; %30
cut2 = val*10/dt; %300

for j = 1:num_to2
    if dwell_t_next_state2(j) <= cut1
        to2_time_bins(1) = to2_time_bins(1)+1;
    elseif (dwell_t_next_state2(j) > cut1) && (dwell_t_next_state2(j) <= cut2)
        to2_time_bins(2) = to2_time_bins(2)+1;
    elseif dwell_t_next_state2(j) > cut2
        to2_time_bins(3) = to2_time_bins(3)+1;
    end
end

for j = 1:num_to3
    if dwell_t_next_state3(j) <= cut1
        to3_time_bins(1) = to3_time_bins(1)+1;
    elseif (dwell_t_next_state3(j) > cut1) && (dwell_t_next_state3(j) <= cut2)
        to3_time_bins(2) = to3_time_bins(2)+1;
    elseif dwell_t_next_state3(j) > cut2
        to3_time_bins(3) = to3_time_bins(3)+1;
    end
end

data = [to2_time_bins to3_time_bins];
labels = [strcat("<", num2str(cut1)) strcat(num2str(cut1),"-",strcat(num2str(cut2))) strcat(">", num2str(cut2))];
labels_cat = categorical(labels);
labels_cat = reordercats(labels_cat,labels);

fig = figure('position', [0, 0, 120, 160]); hold on;
h=bar(labels_cat,data,'stacked');
h(1).FaceColor = purpleq;
h(2).FaceColor = pinkr;
xticklabels([]);
ylim([0 100]);


%%

% other optional colors: blue [0, 0.4470, 0.7410] and red [0.6350, 0.0780, 0.1840]
fig = figure('position', [0, 0, 150, 200]); hold on;
data_sum = sum(data,2);
plot(labels_cat,data(:,1)./data_sum,'-o','color',purpleq,'linewidth',2);
scatter(labels_cat,data(:,1)./data_sum,[],purpleq,'filled');
plot(labels_cat,data(:,2)./data_sum,'-o','color',pinkr,'linewidth',2);
scatter(labels_cat,data(:,2)./data_sum,[],pinkr,'filled');
ylim([0 1]);

%% path to Quiescence
ct_from_reversal = 0;
ct_from_forward = 0;
for j = 2:length(fp_vec)
    if fp_vec(j)==2
        if fp_vec(j-1)==1
            ct_from_forward = ct_from_forward+1;
        elseif fp_vec(j-1)==3
            ct_from_reversal = ct_from_reversal+1;
        end
    end
end

ct_from_forward = [ct_from_forward; 0];
ct_from_reversal = [ct_from_reversal; 0];

data2 = [ct_from_forward ct_from_reversal];
fig = figure('position', [0, 0, 200, 200]); hold on;
h=bar(data2,'stacked');
h(1).FaceColor = purpleq;
h(2).FaceColor = pinkr;

% set 3 display names for the 3 handles
set(h, {'DisplayName'}, {'Forward to Q','Reversal to Q'}')
% Legend will show names for each color
legend('location','best');

%% Histogram of dwell times for states
%{
range = [0 2000];
figure();
subplot(2,2,1);
histogram(dwell_times{1});
title('State 1 dwell time');
xlim(range);

subplot(2,2,2);
histogram(dwell_times{2});
title('State 2 dwell time');
xlim(range);

subplot(2,2,3);
histogram(dwell_times{3});
title('State 3 dwell time');
xlim(range);

subplot(2,2,4);
histogram(dwell_times{4});
title('State 4 dwell time');
xlim(range);
%}
%% Dwell time distribution of State 1 (forward)

fig = figure('position', [0, 0, 200, 150]); hold on;
histogram((3/130)*dwell_times{1},20,'FaceColor',[0.2 0.2 0.2],'EdgeColor',[0.2 0.2 0.2]);
xlim([0 60]);
%title('State 1 dwell time');
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
%print(fig,'figures/Nich_state1_dwell_time','-dpdf','-r800')

%% Dwell time distribution of State 3 (Reversal)

fig = figure('position', [0, 0, 200, 150]); hold on;
histogram((3/130)*dwell_times{3},10,'FaceColor',[0.2 0.2 0.2],'EdgeColor',[0.2 0.2 0.2]);
xlim([0 30]);
%title('State 3 dwell time');
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/Nich_state3_dwell_time','-dpdf','-r800')


%%
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


%%
Markov_mat = NaN(4,4);

ct1 = length(fp1_after);
ct2 = length(fp2_after);
ct3 = length(fp3_after);
ct4 = length(fp4_after);

Markov_mat(1,:) = [sum(fp1_after==1)/ct1 sum(fp1_after==2)/ct1 sum(fp1_after==3)/ct1 sum(fp1_after==4)/ct1];
Markov_mat(2,:) = [sum(fp2_after==1)/ct2 sum(fp2_after==2)/ct2 sum(fp2_after==3)/ct2 sum(fp2_after==4)/ct2];
Markov_mat(3,:) = [sum(fp3_after==1)/ct3 sum(fp3_after==2)/ct3 sum(fp3_after==3)/ct3 sum(fp3_after==4)/ct3];
Markov_mat(4,:) = [sum(fp4_after==1)/ct4 sum(fp4_after==2)/ct4 sum(fp4_after==3)/ct4 sum(fp4_after==4)/ct4];


%% plot Markov matrix
%fig = figure('position', [0, 0, 600, 400]); hold on;
%fig = figure('position', [0, 0, 250, 200]); hold on;

fig = figure('position', [0, 0, 230, 200]); hold on;
%figure();
xlab = [1 2 3 4];
ylab = [1 2 3 4];

colormap gray

cmap = colormap;
cmap = flip(cmap);
colormap(cmap)

imagesc(xlab,ylab,Markov_mat);

%pbaspect([1 1 1]);
set(gca,'xaxisLocation','top')
xticks([1 2 3 4]);
yticks([1 2 3 4]);
colorbar
%caxis([0 1]);
%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/Nich_Markov_mat','-dpdf','-r800')
%}

%% plot sawtooth graph
%{
xtest = -2:0.01:2;
ytest=[];
bias = 0.0;
for j = xtest
    ytest = [ytest sawtooth(j+bias)-bias];
end
figure();
plot(xtest,ytest);
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot wave function  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    b = 0.41; %0.4 % -1 <-> +1
    s1 = 5;
    s2 = 30;
    s3 = 70;
    L1 = 0.4;
    L2 = 0.05;
xtest = 0:0.001:1;
ytest=[];
for j = xtest
    ytest = [ytest wave3(j,b,s1,s2,s3,L1,L2)];
end

fig = figure('position', [0, 0, 100, 120]); hold on;
plot(xtest,0*ytest,'linewidth',2,'color',[0.75 0.75 0.75]);
plot(xtest,ytest,'linewidth',2,'color',[0.5 0.5 0.5]);
scatter(xtest(1:15:end),ytest(1:15:end),60,'.','k');
set(gca,'xticklabel',{[]},'yticklabel',{[]})
axis off
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
%print(fig,'figures/Nich_wave_plain','-dpdf','-r800');
%
fig = figure('position', [0, 0, 150, 100]); hold on;
histogram(ytest,30,'FaceColor',[0.2 0.2 0.2],'EdgeColor',[0.2 0.2 0.2]);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
%print(fig,'figures/Nich_wavedist_plain','-dpdf','-r800')

%% test

xtest = -1:0.001:1;
ytest=[];
for j = xtest
    ytest = [ytest wave3(j,b,s1,s2,s3,L1,L2)];
end

fig = figure('position', [0, 0, 100, 120]); hold on;
plot(xtest,0*ytest,'linewidth',2,'color',[0.75 0.75 0.75]);
plot(xtest,ytest,'linewidth',2,'color',[0.5 0.5 0.5]);
scatter(xtest(1:15:end),ytest(1:15:end),60,'.','k');
set(gca,'xticklabel',{[]},'yticklabel',{[]})
axis off

%% plot weights for functions
%{
figure();
x = -8:0.1:8;
y = -8:0.1:8;
[X,Y] = meshgrid(x,y);

Z = W1(X,Y) + W2(X,Y) + W3(X,Y) + W121(X,Y) + W122(X,Y) + W123(X,Y) +....
    W211(X,Y) + W212(X,Y) + W213(X,Y) + W231(X,Y) + W232(X,Y) + W233(X,Y) +....
       W321(X,Y) + W322(X,Y) + W323(X,Y) +....
       W311(X,Y) + W312(X,Y) + W313(X,Y) + W314(X,Y) + W315(X,Y) + W316(X,Y) + W317(X,Y)+....
       W131(X,Y) + W132(X,Y) + W133(X,Y) + W134(X,Y) + W135(X,Y) + W136(X,Y) + W137(X,Y);

surf(X,Y,Z);
xlabel('x');
ylabel('y');
%}

%% functions

function dxdt = F(x,A1,A2,A3,b1,b2,b3)
    alpha1 = W1(x(1),x(2));
    alpha2 = W2(x(1),x(2));
    alpha121 = W121(x(1),x(2));
    alpha122 = W122(x(1),x(2));
    alpha123 = W123(x(1),x(2));
    alpha211 = W211(x(1),x(2));
    alpha212 = W212(x(1),x(2));
    alpha213 = W213(x(1),x(2));
    alpha3 = W3(x(1),x(2));
    alpha231 = W231(x(1),x(2));
    alpha232 = W232(x(1),x(2));
    alpha233 = W233(x(1),x(2));
    alpha321 = W321(x(1),x(2));
    alpha322 = W322(x(1),x(2));
    alpha323 = W323(x(1),x(2));
    
    alpha311 = W311(x(1),x(2));
    alpha312 = W312(x(1),x(2));
    alpha313 = W313(x(1),x(2));
    alpha314 = W314(x(1),x(2));
    alpha315 = W315(x(1),x(2));
    alpha316 = W316(x(1),x(2));
    alpha317 = W317(x(1),x(2));
    
    alpha131 = W131(x(1),x(2));
    alpha132 = W132(x(1),x(2));
    alpha133 = W133(x(1),x(2));
    alpha134 = W134(x(1),x(2));
    alpha135 = W135(x(1),x(2));
    alpha136 = W136(x(1),x(2));
    alpha137 = W137(x(1),x(2));
    
    dxdt = alpha1*(A1*x+b1) + alpha2*(A2*x+b2) + alpha3*(A3*x+b3) +....
        alpha121*trans121(x) + alpha122*curve12(x) + alpha123*trans122(x)+....
        alpha211*trans211(x) + alpha212*curve21(x) + alpha213*trans212(x)+....
        alpha231*trans231(x) + alpha232*curve23(x) + alpha233*trans232(x)+....
        alpha321*trans321(x) + alpha322*curve32(x) + alpha323*trans322(x)+....
        alpha311*trans311(x) + alpha312*curve312(x) + alpha313*trans313(x)+....
        alpha314*curve314(x) + alpha315*trans315(x) + alpha316*curve316(x)+....
        alpha317*trans317(x)+....
        alpha131*trans131(x) + alpha132*curve132(x) + alpha133*trans133(x)+....
        alpha134*curve134(x) + alpha135*trans135(x) + alpha136*curve136(x)+....
        alpha137*trans137(x);
end


function dxdt = F_perturb(dt,x,A1,A2,A3,b1,b2,b3)
    alpha1 = W1(x(1),x(2));
    alpha2 = W2(x(1),x(2));
    alpha121 = W121(x(1),x(2));
    alpha122 = W122(x(1),x(2));
    alpha123 = W123(x(1),x(2));
    alpha211 = W211(x(1),x(2));
    alpha212 = W212(x(1),x(2));
    alpha213 = W213(x(1),x(2));
    alpha3 = W3(x(1),x(2));
    alpha231 = W231(x(1),x(2));
    alpha232 = W232(x(1),x(2));
    alpha233 = W233(x(1),x(2));
    alpha321 = W321(x(1),x(2));
    alpha322 = W322(x(1),x(2));
    alpha323 = W323(x(1),x(2));
    
    alpha311 = W311(x(1),x(2));
    alpha312 = W312(x(1),x(2));
    alpha313 = W313(x(1),x(2));
    alpha314 = W314(x(1),x(2));
    alpha315 = W315(x(1),x(2));
    alpha316 = W316(x(1),x(2));
    alpha317 = W317(x(1),x(2));
    
    alpha131 = W131(x(1),x(2));
    alpha132 = W132(x(1),x(2));
    alpha133 = W133(x(1),x(2));
    alpha134 = W134(x(1),x(2));
    alpha135 = W135(x(1),x(2));
    alpha136 = W136(x(1),x(2));
    alpha137 = W137(x(1),x(2));
    
    dxdt = alpha1*(A1*x+b1) + alpha2*(A2*x+b2) + alpha3*(A3*x+b3) +....
        alpha121*trans121(x) + alpha122*curve12(x) + alpha123*trans122(x)+....
        alpha211*trans211(x) + alpha212*curve21(x) + alpha213*trans212(x)+....
        alpha231*trans231(x) + alpha232*curve23(x) + alpha233*trans232(x)+....
        alpha321*trans321(x) + alpha322*curve32(x) + alpha323*trans322(x)+....
        alpha311*trans311(x) + alpha312*curve312(x) + alpha313*trans313(x)+....
        alpha314*curve314(x) + alpha315*trans315(x) + alpha316*curve316(x)+....
        alpha317*trans317(x)+....
        alpha131*trans131(x) + alpha132*curve132(x) + alpha133*trans133(x)+....
        alpha134*curve134(x) + alpha135*trans135(x) + alpha136*curve136(x)+....
        alpha137*trans137(x);
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%  Universal parameters  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% TO FP1 %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%% fit 1 - Fig5M %%%%%
    amp = 40; % 
    bias = -0.25;
    s1 = 55;
    s2 = 55;
    s3 = 12;
    L1 = 0.15;
    L2 = 0.32;
    %%%%%%%%%%%%%%%

    %%% sine kink fp2 to fp1 %%%
    xs1 = [-3; -2]; % starting point for sine wave
    dxdt_tmp = F(xs1,A1,A2,A3,b1,b2,b3);
    xe1 = xs1 + dt*dxdt_tmp; % ending point for sine wave
    if (x(1)>min(xs1(1),xe1(1))-1) && (x(1)<max(xs1(1),xe1(1))+1) && (x(2)>min(xs1(2),xe1(2))) && (x(2)<max(xs1(2),xe1(2)))
        L = abs(xs1(2)-xe1(2));
        extra_x = amp*(wave3((1/L)*(x(2)-xs1(2)),bias,s1,s2,s3,L1,L2)); % use wave instead
        extra_y = 0;
        dxdt = dxdt+[extra_x; extra_y];
    end
    
    %%% sine kink fp3 to fp1 %%%
    xs1 = [-3; -4]; % starting point for sine wave
    dxdt_tmp = F(xs1,A1,A2,A3,b1,b2,b3);
    xe1 = xs1 + dt*dxdt_tmp; % ending point for sine wave
    if (x(1)>min(xs1(1),xe1(1))-1) && (x(1)<max(xs1(1),xe1(1))+1) && (x(2)>min(xs1(2),xe1(2))) && (x(2)<max(xs1(2),xe1(2)))
        L = abs(xs1(2)-xe1(2));
        extra_x = amp*(wave3((1/L)*(x(2)-xs1(2)),bias,s1,s2,s3,L1,L2)); % use wave instead
        extra_y = 0;
        dxdt = dxdt+[extra_x; extra_y];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% TO FP2 %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    amp = 20; % 2
    bias=-0.15;
    s1 = 15;
    s2 = 10;
    s3 = 15;
    L1 = 0.2;
    L2 = 0.2;
    
    %%% sine kink fp1 to fp2 %%%
    xs1 = [0; -1]; % starting point for sine wave
    dxdt_tmp = F(xs1,A1,A2,A3,b1,b2,b3);
    xe1 = xs1 + dt*dxdt_tmp; % ending point for sine wave
    if (x(1)>min(xs1(1),xe1(1))-1) && (x(1)<max(xs1(1),xe1(1))+1) && (x(2)>min(xs1(2),xe1(2))) && (x(2)<max(xs1(2),xe1(2)))
        L = abs(xs1(2)-xe1(2));
        extra_x = amp*(wave3((1/L)*(x(2)-xs1(2)),bias,s1,s2,s3,L1,L2)); % use wave instead
        extra_y = 0;
        dxdt = dxdt+[extra_x; extra_y];
    end
    
    %%% sine kink fp3 to fp2 %%%
    xs1 = [0; 1]; % starting point for sine wave
    dxdt_tmp = F(xs1,A1,A2,A3,b1,b2,b3);
    xe1 = xs1 + dt*dxdt_tmp; % ending point for sine wave
    if (x(1)>min(xs1(1),xe1(1))-1) && (x(1)<max(xs1(1),xe1(1))+1) && (x(2)>min(xs1(2),xe1(2))) && (x(2)<max(xs1(2),xe1(2)))
        L = abs(xs1(2)-xe1(2));
        extra_x = amp*(wave3((1/L)*(x(2)-xs1(2)),bias,s1,s2,s3,L1,L2)); % use wave instead
        extra_y = 0;
        dxdt = dxdt+[extra_x; extra_y];
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% TO FP3 %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    amp = 20; % 2
    bias = 0.41; %0.4 % -1 <-> +1
    s1 = 5;
    s2 = 30;
    s3 = 70;
    L1 = 0.4;
    L2 = 0.05;
    
    %%% sine kink fp1 to fp3 %%%
    xs1 = [3; 4]; % starting point for sine wave
    dxdt_tmp = F(xs1,A1,A2,A3,b1,b2,b3);
    xe1 = xs1 + dt*dxdt_tmp; % ending point for sine wave
    if (x(1)>min(xs1(1),xe1(1))-1) && (x(1)<max(xs1(1),xe1(1))+1) && (x(2)>min(xs1(2),xe1(2))) && (x(2)<max(xs1(2),xe1(2)))
        L = abs(xs1(2)-xe1(2));
        extra_x = amp*(wave3((1/L)*(x(2)-xs1(2)),bias,s1,s2,s3,L1,L2)); % use wave instead
        extra_y = 0;
        dxdt = dxdt+[extra_x; extra_y];
    end
    
    %%% sine kink fp2 to fp3 %%%
    xs1 = [3; 2]; % starting point for sine wave
    dxdt_tmp = F(xs1,A1,A2,A3,b1,b2,b3);
    xe1 = xs1 + dt*dxdt_tmp; % ending point for sine wave
    if (x(1)>min(xs1(1),xe1(1))-1) && (x(1)<max(xs1(1),xe1(1))+1) && (x(2)>min(xs1(2),xe1(2))) && (x(2)<max(xs1(2),xe1(2)))
        L = abs(xs1(2)-xe1(2));
        extra_x = amp*(wave3((1/L)*(x(2)-xs1(2)),bias,s1,s2,s3,L1,L2)); % use wave instead
        extra_y = 0;
        dxdt = dxdt+[extra_x; extra_y];
    end    
end


%% Distance
function d = Dist(x)
    fp1 = [-3; -3];   % light blue
    fp2 = [0; 0];     % dark blue
    fp3 = [3; 3];     % yellow
    fp4 = [0; -6];     % magenta

    dfp1 = norm(x-fp1);
    dfp2 = norm(x-fp2);
    dfp3 = norm(x-fp3);
    %dfp4 = norm(x-fp4);
    dfp4 = min(abs(x(1)-6),abs(x(2)+6));  % turn region
    
    d = [dfp1; dfp2; dfp3; dfp4];
end

%% Weighting functions

% Weighting functions
function w1 = W1(x,y)  % weighting fp1
    s = 10;
    w1 = -(1/4)*(tanh(s*(y+2))-tanh(s*(y+4))).*(-tanh(s*(x+2))-tanh(-s*(x+4)));
end

function w2 = W2(x,y) % weighting fp2
    s = 10;
    w2 = (1/4)*(tanh(s*(y+1))-tanh(s*(y-1))).*(tanh(s*(x+1))-tanh(s*(x-1)));
end

function w121 = W121(x,y) % weighting for fp1 -> fp2
    s = 10;
    w121 = (1/4)*(tanh(s*(y+2))-tanh(s*(y+4))).*(-tanh(s*(x+2))-tanh(-s*(x+1)));
end

function w122 = W122(x,y) % weighting for fp1 -> fp2
    s = 10;
    w122 = -(1/4)*(tanh(s*(y+2))-tanh(s*(y+4))).*(-tanh(s*(x-1))-tanh(-s*(x+1)));
end

function w123 = W123(x,y) % weighting for fp1 -> fp2
    s = 10;
    w123 = -(1/4)*(tanh(s*(y+1))-tanh(s*(y+2))).*(-tanh(s*(x-1))-tanh(-s*(x+1)));
end

function w211 = W211(x,y) % weighting for fp2 -> fp1
    s = 10;
    w211 = (1/4)*(tanh(s*(y+1))-tanh(s*(y-1))).*(tanh(s*(x+2))-tanh(s*(x+1)));
end

function w212 = W212(x,y) % weighting for fp2 -> fp1
    s = 10;
    w212 = (1/4)*(tanh(s*(y+1))-tanh(s*(y-1))).*(tanh(s*(x+4))-tanh(s*(x+2)));
end

function w213 = W213(x,y) % weighting for fp2 -> fp1
    s = 10;
    w213 = (1/4)*(tanh(s*(y+2))-tanh(s*(y+1))).*(tanh(s*(x+4))-tanh(s*(x+2)));
end

function w3 = W3(x,y) % weighting fp3
    s = 10;
    w3 = -(1/4)*(tanh(s*(y-2))-tanh(s*(y-4))).*(-tanh(s*(x-2))-tanh(-s*(x-4)));
end

function w231 = W231(x,y) % weighting for fp2 -> fp3 
    s = 10;
    w231 = -(1/4)*(tanh(s*(y-1))-tanh(s*(y+1))).*(tanh(s*(x-1))-tanh(s*(x-2)));
end

function w232 = W232(x,y) % weighting for fp2 -> fp3 
    s = 10;
    w232 = -(1/4)*(tanh(s*(y-1))-tanh(s*(y+1))).*(tanh(s*(x-2))-tanh(s*(x-4)));
end

function w233 = W233(x,y) % weighting for fp2 -> fp3 
    s = 10;
    w233 = (1/4)*(tanh(s*(y-1))-tanh(s*(y-2))).*(tanh(s*(x-2))-tanh(s*(x-4)));
end

function w321 = W321(x,y) % weighting for fp2 -> fp3 
    s = 10;
    w321 = -(1/4)*(tanh(s*(y-2))-tanh(s*(y-4))).*(tanh(s*(x-2))-tanh(s*(x-1)));
end

function w322 = W322(x,y) % weighting for fp2 -> fp3 
    s = 10;
    w322 = -(1/4)*(tanh(s*(y-2))-tanh(s*(y-4))).*(tanh(s*(x-1))-tanh(s*(x+1)));
end

function w323 = W323(x,y) % weighting for fp2 -> fp3 
    s = 10;
    w323 = -(1/4)*(tanh(s*(y-1))-tanh(s*(y-2))).*(tanh(s*(x-1))-tanh(s*(x+1)));
end

function w311 = W311(x,y) % weighting for fp3 -> fp1
    s = 10;
    w311 = (1/4)*(tanh(s*(y-2))-tanh(s*(y-4))).*(tanh(s*(x-4))-tanh(s*(x-5)));
end

function w312 = W312(x,y) % weighting for fp3 -> fp1 
    s = 10;
    w312 = (1/4)*(tanh(s*(y-2))-tanh(s*(y-4))).*(tanh(s*(x-5))-tanh(s*(x-7)));
end

function w313 = W313(x,y) % weighting for fp3 -> fp1 
    s = 10;
    w313 = -(1/4)*(tanh(s*(y-2))-tanh(s*(y+5))).*(tanh(s*(x-5))-tanh(s*(x-7)));
end

function w314 = W314(x,y) % weighting for fp3 -> fp1 
    s = 10;
    w314 = -(1/4)*(tanh(s*(y+5))-tanh(s*(y+7))).*(tanh(s*(x-5))-tanh(s*(x-7)));
end

function w315 = W315(x,y) % weighting for fp3 -> fp1 
    s = 10;
    w315 = -(1/4)*(tanh(s*(y+5))-tanh(s*(y+7))).*(tanh(s*(x+2))-tanh(s*(x-5)));
end

function w316 = W316(x,y) % weighting for fp3 -> fp1 
    s = 10;
    w316 = -(1/4)*(tanh(s*(y+5))-tanh(s*(y+7))).*(tanh(s*(x+4))-tanh(s*(x+2)));
end

function w317 = W317(x,y) % weighting for fp3 -> fp1 
    s = 10;
    w317 = (1/4)*(tanh(s*(y+5))-tanh(s*(y+4))).*(tanh(s*(x+4))-tanh(s*(x+2)));
end

function w131 = W131(x,y) % weighting for fp1 -> fp3 
    s = 10;
    w131 = (1/4)*(tanh(s*(y+2))-tanh(s*(y+4))).*(tanh(s*(x+4))-tanh(s*(x+5)));
end

function w132 = W132(x,y) % weighting for fp1 -> fp3 
    s = 10;
    w132 = (1/4)*(tanh(s*(y+2))-tanh(s*(y+4))).*(tanh(s*(x+5))-tanh(s*(x+7)));
end

function w133 = W133(x,y) % weighting for fp1 -> fp3 
    s = 10;
    w133 = -(1/4)*(tanh(s*(y+2))-tanh(s*(y-5))).*(tanh(s*(x+5))-tanh(s*(x+7)));
end

function w134 = W134(x,y) % weighting for fp1 -> fp3 
    s = 10;
    w134 = -(1/4)*(tanh(s*(y-5))-tanh(s*(y-7))).*(tanh(s*(x+5))-tanh(s*(x+7)));
end

function w135 = W135(x,y) % weighting for fp1 -> fp3 
    s = 10;
    w135 = (1/4)*(tanh(s*(y-5))-tanh(s*(y-7))).*(tanh(s*(x+5))-tanh(s*(x-2)));
end

function w136 = W136(x,y) % weighting for fp1 -> fp3 
    s = 10;
    w136 = (1/4)*(tanh(s*(y-5))-tanh(s*(y-7))).*(tanh(s*(x-2))-tanh(s*(x-4)));
end

function w137 = W137(x,y) % weighting for fp1 -> fp3 
    s = 10;
    w137 = (1/4)*(tanh(s*(y-4))-tanh(s*(y-5))).*(tanh(s*(x-2))-tanh(s*(x-4)));
end


%% Connecting functions (heteroclinic orbits)

function dXdt = curve12(x)
    cloc = [-1; -2];
    dXdt = 4*polar_ccw(x,cloc);
end

function dXdt = trans121(x)
    dxdt =  8;
    dydt = -32*(x(2)+3);
    dXdt = [dxdt; dydt];  
end

function dXdt = trans122(x)
    dxdt =  -32*(x(1)+0);
    dydt = 8;
    dXdt = [dxdt; dydt];  
end

function dXdt = curve21(x)
    cloc = [-2; -1];
    dXdt = 4*polar_ccw(x,cloc);
end

function dXdt = trans211(x)
    dxdt =  -8;
    dydt = -32*(x(2)+0);
    dXdt = [dxdt; dydt];  
end

function dXdt = trans212(x)
    dxdt =  -32*(x(1)+3);
    dydt = -8;
    dXdt = [dxdt; dydt];  
end

function dXdt = curve23(x)
    cloc = [2; 1];
    dXdt = 4*polar_ccw(x,cloc);
end

function dXdt = trans231(x)
    dxdt =  8;
    dydt = -32*(x(2)+0);
    dXdt = [dxdt; dydt];  
end

function dXdt = trans232(x)
    dxdt =  -32*(x(1)-3);
    dydt = 8;
    dXdt = [dxdt; dydt];  
end

function dXdt = curve32(x)
    cloc = [1; 2];
    dXdt = 4*polar_ccw(x,cloc);
end

function dXdt = trans321(x)
    dxdt =  -8;
    dydt = -32*(x(2)-3);
    dXdt = [dxdt; dydt];  
end

function dXdt = trans322(x)
    dxdt =  -32*(x(1)-0);
    dydt = -8;
    dXdt = [dxdt; dydt];  
end

function dXdt = trans311(x)
    dxdt =  8;
    dydt = -16*(x(2)-3);
    dXdt = [dxdt; dydt]; 
end

function dXdt = curve312(x)
    cloc = [5; 2];
    dXdt = 2*polar_cw(x,cloc);
end

function dXdt = trans313(x)
    dxdt =  -16*(x(1)-6);
    dydt = -8;
    dXdt = [dxdt; dydt];  
end

function dXdt = curve314(x)
    cloc = [5; -5];
    dXdt = 2*polar_cw(x,cloc);
end

function dXdt = trans315(x)
    dxdt =  -8;
    dydt = -32*(x(2)+6);
    dXdt = [dxdt; dydt];  
end

function dXdt = curve316(x)
    cloc = [-2; -5];
    dXdt = polar_cw(x,cloc);
end

function dXdt = trans317(x)
    dxdt =  -16*(x(1)+3);
    dydt = 4;
    dXdt = [dxdt; dydt];
end

function dXdt = trans131(x)
    dxdt =  -40;     % make faster
    dydt = -40*(x(2)+3);
    dXdt = [dxdt; dydt];
end

function dXdt = curve132(x)
    cloc = [-5; -2];
    dXdt = 8*polar_cw(x,cloc);
end

function dXdt = trans133(x)
    dxdt =  -40*(x(1)+6);
    dydt = 40;    % make faster
    dXdt = [dxdt; dydt];  
end

function dXdt = curve134(x)
    cloc = [-5; 5];
    dXdt = 8*polar_cw(x,cloc);
end

function dXdt = trans135(x)
    dxdt =  40;    % make faster
    dydt = -40*(x(2)-6);
    dXdt = [dxdt; dydt];  
end

function dXdt = curve136(x)
    cloc = [2; 5];
    dXdt = 2*polar_cw(x,cloc);
end

function dXdt = trans137(x)
    dxdt =  -40*(x(1)-3);
    dydt = -20;      % make faster
    dXdt = [dxdt; dydt];  
end



%%
function dXdt = polar_cw(x,cloc)
    r = sqrt((x(1)-cloc(1))^2 + (x(2)-cloc(2))^2);
    theta = atan2(x(2)-cloc(2),x(1)-cloc(1));
    drdt = 10*r*(1-r);
    dthetadt= -3;
    dxdt = -r*sin(theta)*dthetadt+cos(theta)*drdt;
    dydt = r*cos(theta)*dthetadt + sin(theta)*drdt;
    dXdt = [dxdt; dydt];  
end

function dXdt = polar_ccw(x,cloc)
    r = sqrt((x(1)-cloc(1))^2 + (x(2)-cloc(2))^2);
    theta = atan2(x(2)-cloc(2),x(1)-cloc(1));
    drdt = 10*r*(1-r);
    dthetadt= 3;
    dxdt = -r*sin(theta)*dthetadt+cos(theta)*drdt;
    dydt = r*cos(theta)*dthetadt + sin(theta)*drdt;
    dXdt = [dxdt; dydt];  
end

%%
function saw = sawtooth(x)
    saw = zeros(1,length(x));
    for i = 1:length(x)
        if x(i) >=-1.5 && x(i)<-1.25
            saw(i) = -1.5-x(i);
        elseif x(i) >=-1.25 && x(i)<-3/4
            saw(i) = 1+x(i);
        elseif x(i) >=-3/4 && x(i)<-1/4
            saw(i) = -1/2-x(i);
        elseif x(i) >=-1/4 && x(i)<1/4
            saw(i) = x(i);
        elseif x(i) >=1/4 && x(i)<3/4
            saw(i) = 1/2-x(i);
        elseif x(i) >=3/4 && x(i)<=1.25
            saw(i) = -1+x(i);
        elseif x(i) >=1.25 && x(i)<=1.5
            saw(i) = 1.5-x(i);
        end
    end
    saw = 4*saw;
end

%%%%%%%%%%%%%%%%%%%%%%%
%%%%% wave function %%%
%%%%%%%%%%%%%%%%%%%%%%%
%{
function wave = wave(x,b)
    % b 0->1/2, bias term ranges from 0 to 0.5
    s = 10;
    %wave = tanh(s*x)-(tanh(s*(x-(1/2+b)))+1)+1*(tanh(s*(x-1))+1)-(tanh(s*(x+(1/2-b)))-1)+....
    %    1*(tanh(s*(x+1))-1)-(tanh(s*(x-3/2))+1)-(tanh(s*(x+3/2))-1);
    wave = tanh(s*x)-(tanh(s*(x-(1/2+b)))+1)+1*(tanh(s*(x-1))+1)-(tanh(s*(x+(1/2-b)))-1)+1*(tanh(s*(x+1))-1);
end

function wave = wave2(x,b)
    % b 0->1/2, bias term ranges from 0 to 0.5
    s1 = 20;
    s2 = 50;
    end_len = 0.25;
    wave = 1/2*(tanh(s2*(x-end_len))+1)-tanh(s2*(x-(0.5+b)))-1+1/2*(tanh(s1*(x-(1-end_len)))+1)+1/2*(tanh(s1*(x+end_len))+1)-tanh(s2*(x+(0.5-b)))-1+1/2*(tanh(s2*(x+(1-end_len)))+1);
end
%}

%%%%%%%%%%%%%%%%%%%%%%%
%%%%% wave function %%%
%%%%%%%%%%%%%%%%%%%%%%%

function wave = wave3(x,b,s1,s2,s3,L1,L2)
    if x>=0
        wave = (1/2)*(tanh(s1*(x-L1)) -tanh(s1*(-x-L1)))-tanh(s2*(x-(0.5+b)))-1+1/2*(tanh(s3*(x-(1-L2)))+1)+1/2*(tanh(s3*(x-(1+L2)))+1);
    else
        wave = (1/2)*(tanh(s1*(-x-L1)) -tanh(s1*(x-L1)))-tanh(s2*(-x-(0.5+b)))-1+1/2*(tanh(s3*(-x-(1-L2)))+1)+1/2*(tanh(s3*(-x-(1+L2)))+1);
    end
end

function wave = wave_test(x,b,s1,s2,s3,L1,L2)
    if x>=0
        wave = (1/2)*(tanh(s1*(x-L1)) -tanh(s1*(-x-L1)))-tanh(s2*(x-(0.5+b)))+1/2*(tanh(s3*(x-(1-L2))))+1/2*(tanh(s3*(x-(1+L2))));
    else
        wave = (1/2)*(tanh(s1*(-x-L1)) -tanh(s1*(x-L1)))-tanh(s2*(-x-(0.5+b)))+1/2*(tanh(s3*(-x-(1-L2))))+1/2*(tanh(s3*(-x-(1+L2))));
    end
end