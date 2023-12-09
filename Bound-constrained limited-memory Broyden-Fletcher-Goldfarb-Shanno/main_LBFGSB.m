% Code package purpose: Accoplish the proposed method in the paper entitled
%                       "Near-surface Rayleigh Wave Dispersion Curve 
%                       Inversion Algorithms: A Comprehensive Comparison".
%
% paper status: major revision (journal: Surveys in Geophysics)  
%
% software version: MATLAB R2017a
%
% Acknowledgement: The forward modeling program used to generate 
%                  theoretical Rayleigh wave dispersion curves in this 
%                  code package was obtained from the  website 
%                  (https://github.com/eespr/MuLTI) provided by 
%                  Killingbeck et al. (2018);the L-BFGS-B codes available 
%               on the website (https://github.com/stephenbeckr/L-BFGS-B-C)
%                  by Stephen Becker were applied for this code package.
%
% Killingbeck et al. (2018): Killingbeck, S. F., Livermore, P. W., 
%                            Booth, A. D., & West, L. J. (2018). Multimodal 
%                            layered transdimensional inversion of seismic 
%                            dispersion curves with depth constraints. 
%                            Geochemistry, Geophysics, Geosystems, 19(12), 
%                            4957-4971.
%
% Date: 2023/11/26
%
% Developed by: Xiao-Hui Yang, Currently working at 
%               Chengdu University of Information Technology
%
% Email: yangxh@cuit.edu.cn / xiao-hui.yang@hotmail.com
%
% Note: the inversion performance of the L-BFGS-B method for Rayleigh wave 
%       inversion can refer to the paper entitiled "Near-surface Rayleigh 
%       Wave Dispersion Curve Inversion Algorithms: A Comprehensive 
%       Comparison"; the users can cite this paper for scientific research.
% 

clear;
clc;
close all;

myFontSize = 20;
myMarkerSize = 20;

%% Reset the seed by clock for random number generation
rand('seed',sum(100*clock))
randn('seed',sum(100*clock))

%% Raw data (a numerical example, measured dispersion curves)
curve_00 = xlsread('numerical_fun.xls'); % fundamental mode
curve_01 = xlsread('numerical_1st.xls'); % 1st higher mode
curve_02 = xlsread('numerical_2nd.xls'); % 2nd higher mode
curve_03 = xlsread('numerical_3rd.xls'); % 3rd higher mode

% fundamental mode
f_00_original = curve_00(:,1)'; f_00_original = f_00_original(:)';
dispersion_00_original = curve_00(:,2); 
dispersion_00_original = dispersion_00_original(:)';
% first higher mode
f_01_original = curve_01(:,1)'; f_01_original = f_01_original(:)';
dispersion_01_original = curve_01(:,2); 
dispersion_01_original = dispersion_01_original(:)';
% second higher mode
f_02_original = curve_02(:,1); f_02_original = f_02_original(:)';
dispersion_02_original = curve_02(:,2); 
dispersion_02_original = dispersion_02_original(:)';
% third higher mode
f_03_original = curve_03(:,1); f_03_original = f_03_original(:)';
dispersion_03_original = curve_03(:,2); 
dispersion_03_original = dispersion_03_original(:)';

% interpolation process - frequency and dispersion values for invertion use
f_00_min = 6; f_00_max = 48;
f_01_min = 31.5; f_01_max = 69.5;
f_02_min = 29.5; f_02_max = 49;
f_03_min = 64; f_03_max = 100;
df = 0.5;
f_00 = f_00_min:df:f_00_max;
f_01 = f_01_min:df:f_01_max;
f_02 = f_02_min:df:f_02_max;
f_03 = f_03_min:df:f_03_max;
dispersion_00 = interp1(f_00_original,dispersion_00_original,f_00);
dispersion_01 = interp1(f_01_original,dispersion_01_original,f_01);
dispersion_02 = interp1(f_02_original,dispersion_02_original,f_02);
dispersion_03 = interp1(f_03_original,dispersion_03_original,f_03);

dispersion_all_cell = cell(1,4);
dispersion_all_cell{1} = dispersion_00;
dispersion_all_cell{2} = dispersion_01;
dispersion_all_cell{3} = dispersion_02;
dispersion_all_cell{4} = dispersion_03;


f = 6:df:100;
index_vec_all = cell(1,4);
index_vec_all{1} = [find(f == f_00_min) find(f == f_00_max)];
index_vec_all{2} = [find(f == f_01_min) find(f == f_01_max)];
index_vec_all{3} = [find(f == f_02_min) find(f == f_02_max)];
index_vec_all{4} = [find(f == f_03_min) find(f == f_03_max)];

modes_num_vec = [1 2 3 4]; % available modes of dispersion curves
index_vec = index_vec_all(modes_num_vec);

dispersions_R_true = []; % integrate different modes of dispersion curves
for jj = 1:1:length(modes_num_vec)
    temp = modes_num_vec(jj);
    dispersions_R_true = [dispersions_R_true dispersion_all_cell{temp}];
end

%% Actual model parameters (a numerical example)
h_true = [2 4 5 5];
h_true2 = [h_true 0];
Vs_true = [400 200 300 500 650];
Vp_true = [700 300 500 900 1100]; % primary wave velocity (all layers)
den_true = [1.9 1.7 1.8 2.0 2.1]; % density (all layers)

layers_num = length(Vp_true);

Vs_true_profile = [h_true Vs_true];

Vs_profile_lower = [1.3 2.6 3.25 3.25 260 130 195 325 500]; % search space
Vs_profile_upper = [2.7 5.4 6.75 6.75 540 350 405 675 780]; % search space

%% L-BFGS-B inversion
Vp = Vp_true;
den = den_true;

num_unknown = length(Vs_profile_lower);
Y_raw = zeros(1,num_unknown);
for j = 1:num_unknown
    Y_raw(j) = Vs_profile_lower(j) + (Vs_profile_upper(j)-...
        Vs_profile_lower(j))*rand(1,1);
end

theta_0 = Y_raw;


my_fun = @(theta_BFGS)BFGSobjfun(theta_BFGS, layers_num, ...
    dispersions_R_true, Vp, den, f,modes_num_vec,index_vec_all,...
    Vs_profile_lower,Vs_profile_upper);
my_grad = @(theta_BFGS)my_gradFun(theta_BFGS,layers_num, ...
    dispersions_R_true, Vp, den, f,modes_num_vec,index_vec_all);


% Solve NNLS with L-BFGS-B
L_fun = @(theta_BFGS)fminunc_wrapper( theta_BFGS, my_fun, my_grad); 

tstart=tic;
% Request very high accuracy for this test:
opts    = struct( 'factr', 1e4, 'pgtol', 1e-8, 'm', 10);
opts.printEvery     = 5;
opts.x0 = theta_0(:); % initial value
opts.pgtol = 1e-6; % tolerance
opts.maxIts = 350; % max iterations
opts.maxTotalIts = 500000;
if length(Vs_profile_upper) > 10000
    opts.m  = 50;
end
% Run the algorithm:
global time_global_index Obj_global_index
time_global_index = [];
Obj_global_index = [];

tic
[Y_hat, ~, info] = lbfgsb(L_fun, Vs_profile_lower(:), ...
    Vs_profile_upper(:), opts );
Y_hat = Y_hat(:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_index = time_global_index(:);
Obj_temp_index = sqrt(Obj_global_index);

clear time_global_index Obj_global_index

Obj_index = zeros(size(Obj_temp_index,1),size(Obj_temp_index,2));

for i = 1:length(Obj_temp_index)
    Obj_index(i) = min(Obj_temp_index(1:i));
end
all_index = [time_index Obj_index];

dispersions_R_inverted = calDispersions_2(Y_hat,Vp_true,den_true,f,...
    modes_num_vec,index_vec_all);

%% draw plots
% draw inverted dispersion curves
figure(1)
plot(curve_00(:,1),curve_00(:,2),'k.','MarkerSize',myMarkerSize);
hold on
plot(curve_01(:,1),curve_01(:,2),'k.','MarkerSize',myMarkerSize);
plot(curve_02(:,1),curve_02(:,2),'k.','MarkerSize',myMarkerSize);
plot(curve_03(:,1),curve_03(:,2),'k.','MarkerSize',myMarkerSize);
drawDispersionsCompareField_5(dispersions_R_inverted,f,modes_num_vec,...
    index_vec,myFontSize) % inverted
axis([0 100 200 600]);
set(gca,'XTick',0:20:100);
set(gca,'YTick',100:100:700);
set(gca,'FontName','Times New Roman','FontSize',myFontSize);

% draw objective function curve
figure(2)
plot(all_index(:,1),all_index(:,2),'b','Linewidth',1.5);
xlabel('Time [s]','FontSize',myFontSize);
ylabel('Obj','FontSize',myFontSize);
set(gca,'FontName','Times New Roman','FontSize',myFontSize);

% draw Vs profile
figure(3)
drawProfile_fun(Y_hat,h_true,Vs_true)
xlabel('Shear-wave velocity [m/s]','FontSize',myFontSize);
ylabel('Depth [m]','FontSize',myFontSize);
axis([100 800 0 25]);
set(gca,'FontName','Times New Roman','FontSize',myFontSize);
set(gca,'XTick',0:200:800);
set(gca,'YTick',0:5:25);
% set(figure(3),'Position',[680,100,560,880]); %[left, bottom, width, height]
