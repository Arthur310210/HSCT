% example of standard linear chirps
clear all; close all;clc;
addpath('functions');

% Choose the experiment you want to try the following singals:
Expri = 3;

% Parameter settings:
Hz = 100 ; % sampling rate
L = 6 ; % time duration in second
t = (0:1/Hz:L)'; % Input time tick
t_idx = 100:500; % Observe time index
t_show = t(t_idx); % Observe time 
tDs = 1; % Hopping step
Low_Fre = 0; % Min frequecny
High_Fre = 50; % Max frequecny
sigma = 1; % Standard deviation for Gaussian window
Q = round(1.9143/(sqrt(sigma)*(1/Hz))); % Frequency points = 2Q+1
alpha_resol = 1/(2*Q+1);  % Frequency resolution
chirpReso = 0.2; % Chirp resolution
chirpMax = 16.55; % chirp max and resolution: -chirpMax:chirpReso:chirpMax

% Window function settings:
tt = (-Q:Q)*(1/Hz); %同 t % Frequency points = 2Q+1
% windows functions:
alpha1 = 1;
h1 = exp(-pi*alpha1*tt.^2'); % window g_0
Dh1 = dwindow(h1);
DDh1 = dwindow(Dh1);
DDDh1 = dwindow(DDh1);


% Input signal : 

if Expri==1
    % Exp4.1 (heat kernel, nice one.) :x2
    x1 = exp(2*pi*1i*(12*t + 12*exp(-0.5*(t-3).^2)));
    x2 = 0;      
    if1 = 12-6*(t-3).*exp(-0.5*(t-3).^2);
    if2 = 0;
elseif Expri==2
    %Exp4.2, 4.3 (Pure chirp) :x11,x12
    x1 = exp(2*pi*1i*(4*t.^2)) ;
    x2 = exp(2*pi*1i*(-pi*t.^2 + (24+6*pi).*t));
    if1 = 8*t;
    if2 = -2*pi*t+ (24+6*pi);
elseif Expri==3
    % Exp4.2 (sine) : x31,x32
    x1 = exp(2*pi*1i*(8*(t-2.2).^2));
    x2 = exp(2*pi*1i*(-2*cos((2/3)*pi*1*(t-2))+ 13*(t-2)));
    if1 = 16*(t-2.2);
    if2 = (4/3*pi)*sin((2/3)*pi*1*(t-2))+13;
elseif Expri==4
    % Exp4.4 : Three components case, 4-2 is used now. : x41, x42, x43
    %Exp4-2 (Fixed, 2D, 3D both have gap, nice.)
    x1 = exp(2*pi*1i*((5/6)*t.^3-(15/2)*t.^2+ 40*t+100)) + exp(2*pi*1i*(-(5/6)*t.^3+ (15/2)*t.^2+(35/2)*t));
    x2 = exp(2*pi*1i*(-3.5*t.^2 + (26+6*pi).*t));
    if1 = (5/2)*t.^2-(15)*t+ 40;
    if2 = -(5/2)*t.^2+ (15)*t+(35/2);
else
    disp('No experiment settings');
end
% Sum of signal:
x = x1 + x2;

% run CT,SCT,Proposed for different windows
[tfc1, tfrsq1, tfrsq1_proposed, tcrtic1, tfrtic1, ttrtic1] = ...
    TF_algorithm(x, Hz, Low_Fre, High_Fre, chirpMax, alpha_resol, tDs, chirpReso, h1, Dh1, DDh1, DDDh1);

return;

%% plot the upsampled signal and their IFs:
figure()
plot(t,if1,'LineWidth',2);
hold on
plot(t,if2,'LineWidth',2);
xlabel('time (s)');
ylabel('frequency (Hz)')
ylim([0 50])
set(gca,'fontsize',35);title('Instantaneous frequency');

%% Project CT, SCT, Proposed with the window g0 on the TF-plane:
[d1,r1,c1] = size(tfc1);
tfc1_2D = zeros(r1,c1); tfrsq1_2D = zeros(r1,c1); tfrsq1_proposed_2D = zeros(r1,c1);
for z =1:d1
    tfc1_2D = tfc1_2D + abs(squeeze(tfc1(z,:,:)));
    tfrsq1_2D = tfrsq1_2D + abs(squeeze(tfrsq1(z,:,:)));
    tfrsq1_proposed_2D = tfrsq1_proposed_2D + abs(squeeze(tfrsq1_proposed(z,:,:)));
end
C=1*1e3;
U = 0.8;
figure()
subplot(221)
image(ttrtic1',tfrtic1,(abs(tfc1_2D)./max(max(abs(tfc1_2D))))*C);
colormap(gray(256));
set(gca,'Ydir','normal');
set(gca,'Fontsize',18);
xlabel('Time (Sec)','Fontsize',20);
ylabel('Frequency (Hz)','Fontsize',20);
title('CT','Fontsize',25);
subplot(222)
image(ttrtic1',tfrtic1,(abs(tfrsq1_2D)./max(max(abs(tfrsq1_2D)))).^U*C);
colormap(gray(256));
set(gca,'Ydir','normal');
set(gca,'Fontsize',18);
xlabel('Time (Sec)','Fontsize',20);
ylabel('Frequency (Hz)','Fontsize',20);
title('SCT','Fontsize',25);
subplot(223)
image(ttrtic1',tfrtic1,(abs(tfrsq1_proposed_2D)./max(max(abs(tfrsq1_proposed_2D)))).^U*C);
colormap(gray(256));
set(gca,'Ydir','normal');
set(gca,'Fontsize',18);
xlabel('Time (Sec)','Fontsize',20);
ylabel('Frequency (Hz)','Fontsize',20);
title('Proposed','Fontsize',25);

tfc1_RE = Renyi_entropy(tfc1,2);
tfrsq1_RE = Renyi_entropy(tfrsq1,2);
tfrsq1_proposed_RE = Renyi_entropy(tfrsq1_proposed,2);
disp('Renyi entropy of CT is');
disp(tfc1_RE);
disp('Renyi entropy of SCT is');
disp(tfrsq1_RE);
disp('Renyi entropy of Proposed is');
disp(tfrsq1_proposed_RE);
%% Chirp-rate vs Magnitude
if Expri==2
    x10 = -6.28;
    x20 = +8;
    str = 'value $-7$';
    str2 = 'value $+12$';
    width = 2;
    t3 = knnsearch(ttrtic1', 3, 'k', 1);
    f24 = knnsearch(tfrtic1, 24, 'k',  1);   
elseif Expri==3
    x10 = -6.64;
    x20 = +16;
    str = 'value $-7$';
    str2 = 'value $+12$';
    width = 2;
    t3 = knnsearch(ttrtic1', 3.16, 'k', 1);
    f24 = knnsearch(tfrtic1, 15.7, 'k',  1);   
else
    x10 = -6.28;
    x20 = +8;
    str = 'value $-7$';
    str2 = 'value $+12$';
    width = 2;
    t3 = knnsearch(ttrtic1', 3, 'k', 1);
    f24 = knnsearch(tfrtic1, 24, 'k',  1);   
end

figure()
set(gcf,'color','w','units','normalized','position',[0 0.5 0.8 0.3])
subplot(131)  
plot(tcrtic1, abs(squeeze(tfc1(:, f24, t3))));
set(gca,'Fontsize',18);
xlabel('Chirp rate (\lambda)','Fontsize',20)
ylabel('Magnitude','Fontsize',20)
%title('3.16, 15.7','Fontsize',20)
grid on;
hold on
plot([x10,x10], [0,1.2], 'R', 'LineWidth', width,'DisplayName', str ,'LineStyle', ':');
plot([x20,x20], [0,1.2], 'r', 'LineWidth', width,'DisplayName', str2,'LineStyle', ':');
%legend('Location', 'southwest', 'FontSize', 10,'Interpreter', 'latex');


subplot(132)
%t3 = knnsearch(ttrtic1', 3, 'k', 1);
%f24 = knnsearch(tfrtic1, 24, 'k',  1);
plot(tcrtic1, abs(squeeze(tfrsq1(:, f24, t3))));
set(gca,'Fontsize',18);
xlabel('Chirp rate (\lambda)','Fontsize',20)
ylabel('Magnitude','Fontsize',20)
%title('SCT','Fontsize',20)
ylim([0,30])
grid on;
hold on
plot([x10,x10], [0,30], 'r', 'LineWidth', width,'DisplayName', str,'LineStyle', ':' );
plot([x20,x20], [0,30], 'r', 'LineWidth', width,'DisplayName', str,'LineStyle', ':' );

subplot(133)
%t3 = knnsearch(ttrtic1', 3, 'k', 1);
%f24 = knnsearch(tfrtic1, 24, 'k',  1);
plot(tcrtic1, abs(squeeze(tfrsq1_proposed(:, f24, t3))));
set(gca,'Fontsize',18);
xlabel('Chirp rate (\lambda)','Fontsize',20)
ylabel('Magnitude','Fontsize',20)
%title('Proposed','Fontsize',20)
ylim([0,30])
grid on;
hold on
plot([x10,x10], [0,30], 'r', 'LineWidth', width,'DisplayName', str,'LineStyle', ':' );
plot([x20,x20], [0,30], 'r', 'LineWidth', width,'DisplayName', str ,'LineStyle', ':');
%% 3d plot of CT with g_0
QN = 5 ;
D = tfc1(:,:,t_idx);
thresh = quantile(abs(D(:)),0.9999);
D(find(abs(D) < thresh * (10-QN+1)/10)) = thresh * (10-QN)/10 ;
 
for jj = 1: QN
    idx = find(abs(D) <= thresh * (10-jj+1)/10 & abs(D) > thresh * (10-jj)/10 );
    [I1,I2,I3] = ind2sub(size(D),idx);
    scatter3((I3+99)/Hz,tfrtic1(I2),tcrtic1(I1),20, [1 1 1]*(jj-1)/8, 'filled');
    hold on
end
 
view(30,30)
ylim([0 50])
zlim([-16.55 25])
set(gca,'Fontsize',23);
xlabel('time (s)','Fontsize',25);
ylabel('frequency (Hz)','Fontsize',25);
zlabel('chirp rate','Fontsize',25);
title('CT','Fontsize',25)
colormap(1-gray)
% colorbar

%% 3d plot of SCT with g_0
QN = 5 ;
D = tfrsq1(:,:,t_idx);
thresh = quantile(abs(D(:)),0.9999);
D(find(abs(D) < thresh * (10-QN+1)/10)) = thresh * (10-QN)/10 ;
 
for jj = 1: QN
    idx = find(abs(D) <= thresh * (10-jj+1)/10 & abs(D) > thresh * (10-jj)/10 );
    [I1,I2,I3] = ind2sub(size(D),idx);
    scatter3((I3+99)/Hz,tfrtic1(I2),tcrtic1(I1),20, [1 1 1]*(jj-1)/8, 'filled');
    hold on
end
 
view(30,30)
ylim([0 50])
zlim([-16.55 25])
set(gca,'Fontsize',23);
xlabel('time (s)','Fontsize',25);
ylabel('frequency (Hz)','Fontsize',25);
zlabel('chirp rate','Fontsize',25);
title('SCT','Fontsize',25)
colormap(1-gray)
% colorbar

%% 3d plot of Proposed with g_0
QN = 5 ;
D = tfrsq1_proposed(:,:,t_idx);
thresh = quantile(abs(D(:)),0.9999);
D(find(abs(D) < thresh * (10-QN+1)/10)) = thresh * (10-QN)/10 ;
 
for jj = 1: QN
    idx = find(abs(D) <= thresh * (10-jj+1)/10 & abs(D) > thresh * (10-jj)/10 );
    [I1,I2,I3] = ind2sub(size(D),idx);
    scatter3((I3+99)/Hz,tfrtic1(I2),tcrtic1(I1),20, [1 1 1]*(jj-1)/8, 'filled');
    hold on
end
 
view(30,30)
ylim([0 50])
zlim([-16.55 25])
set(gca,'Fontsize',23);
xlabel('time (s)','Fontsize',25);
ylabel('frequency (Hz)','Fontsize',25);
zlabel('chirp rate','Fontsize',25);
title('Proposed','Fontsize',25)
colormap(1-gray)
% colorbar



%% Project CT, SCT, Proposed with the window g0 on the TF-plane:
[d1,r1,c1] = size(tfc1);
tfc1_2D = zeros(r1,c1); tfrsq1_2D = zeros(r1,c1); tfrsq1_proposed_2D = zeros(r1,c1);
for z =1:d1
    tfc1_2D = tfc1_2D + squeeze(tfc1(z,:,:));
    tfrsq1_2D = tfrsq1_2D + squeeze(tfrsq1(z,:,:));
    tfrsq1_proposed_2D = tfrsq1_proposed_2D + squeeze(tfrsq1_proposed(z,:,:));
end
C=3*1e3;
U =0.8

% h = rotate3d;
% set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
% set(h, 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')')
% set(gcf, 'ResizeFcn', @align_axislabel)
% align_axislabel([], gca)

figure()
%subplot(221)
image(t,tfrtic1,abs(tfc1_2D)./max(max(abs(tfc1_2D)))*C);
colormap(gray(256));
set(gca,'Ydir','normal');
set(gca,'Fontsize',35);
xlabel('Time (Sec)','Fontsize',40);
ylabel('Frequency (Hz)','Fontsize',40);
%title('CT with the window g0 ','Fontsize',50);
%title(['e^{-\pit^{2}}',num2str(alpha1)])
%subplot(222)
figure()
image(t,tfrtic1,(abs(tfrsq1_2D)./max(max(abs(tfrsq1_2D)))).^U*C);
colormap(gray(256));
set(gca,'Ydir','normal');
set(gca,'Fontsize',35);
xlabel('Time (Sec)','Fontsize',40);
ylabel('Frequency (Hz)','Fontsize',40);
%title('SCT with the window g0','Fontsize',50);
%subplot(223)
figure()
image(t,tfrtic1,(abs(tfrsq1_proposed_2D)./max(max(abs(tfrsq1_proposed_2D)))).^U*C);
colormap(gray(256));
set(gca,'Ydir','normal');
set(gca,'Fontsize',35);
xlabel('Time (Sec)','Fontsize',40);
ylabel('Frequency (Hz)','Fontsize',40);
%title('Proposed method with the window g0','Fontsize',50);


%% projection of CT onto TF plane with g0
tfproj = zeros(size(tfc1,2),size(tfc1,3));
for i = 1:size(tfproj,1)
    for j = 1:size(tfproj,2)
        tfproj(i,j) = sum(abs(squeeze(tfc1(:,i,j))));
    end
end
figure()
imageSQ(t(t_idx), Hz*tfrtic1, abs(tfproj(:,t_idx)), 0.9999); axis xy; colormap(1-gray); 
xlabel('time (s)'); ylabel('frequency (Hz)'); title('CT on the TF-plane with g0'); 

%% projection of SCT onto TF plane with g0
tfproj = zeros(size(tfrsq1,2),size(tfrsq1,3));
for i = 1:size(tfproj,1)
    for j = 1:size(tfproj,2)
        tfproj(i,j) = sum(abs(squeeze(tfrsq1(:,i,j))));
    end
end
figure()
imageSQ(t(t_idx), Hz*tfrtic1, abs(tfproj(:,t_idx)), 0.9999); axis xy; colormap(1-gray); 
xlabel('time (s)'); ylabel('frequency (Hz)'); title('SCT on the TF-plane with g0');

%% projection of Proposed onto TF plane with g0
tfproj = zeros(size(tfrsq1_proposed,2),size(tfrsq1_proposed,3));
for i = 1:size(tfproj,1)
    for j = 1:size(tfproj,2)
        tfproj(i,j) = sum(abs(squeeze(tfrsq1_proposed(:,i,j))));
    end
end
figure()
imageSQ(t(t_idx), Hz*tfrtic1, abs(tfproj(:,t_idx)), 0.9999); axis xy; colormap(1-gray); 
xlabel('time (s)'); ylabel('frequency (Hz)'); title('Proposed on the TF-plane with g0'); 
