% Code package purpose: Accoplish the proposed method in the paper entitled
%                       "Near-surface Rayleigh Wave Dispersion Curve 
%                       Inversion Algorithms: A Comprehensive Comparison".
%
% paper status: major revision (journal: Surveys in Geophysics)  
%
% software version: MATLAB R2017a
%
% Acknowledgement:
%      The Genetic Algorithm Toolbox for MATLAB was developed at the
%      Department of Automatic Control and Systems Engineering of The
%      University of Sheffield, UK, in order to make GA's accessible to the
%      control engineer within the framework of a existing computer-aided
%      control system design package. The toolbox was written with the
%      support of a UK SERC grant, and the final version (v1.2) was
%      completed in 1994.

%      The Toolbox was originally developed for MATLAB v4.2. It has also
%      been successfully used with subsequent versions up to and including
%      MATLAB 7.

%      For a more detailed introduction to the capabilities and use of the 
%      GA Toolbox, please refer to the two introductory papers and the
%      Toolbox User's Guide, all of which are available at the GA Toolbox
%      homepage at http://codem.group.shef.ac.uk/index.php/ga-toolbox.

%      The GA Toolbox is copyright the original authors and The University 
%      of Sheffield, and is published here under the GNU General Public 
%      License. (See http://www.fsf.org/licenses/licenses.html)
% Acknowledgement: The forward modeling program used to generate 
%                  theoretical Rayleigh wave dispersion curves in this 
%                  code package was obtained from the  website 
%                  (https://github.com/eespr/MuLTI) provided by 
%                  Killingbeck et al. (2018)
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
% Note: the inversion performance of the GA algorithm for Rayleigh wave 
%       inversion can refer to the paper entitiled "Near-surface Rayleigh 
%       Wave Dispersion Curve Inversion Algorithms: A Comprehensive 
%       Comparison"; the users can cite this paper for scientific research.
% 

clear;
clc;
close all;

%% add relative path
addpath (genpath('GA_ToolBox')) % add relative path

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

validation_x = dispersions_R_true;

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

%% GA inversion
Vp = Vp_true;
den = den_true;

LIND = 60; % Length of individual vars.
NVAR = 2*layers_num-1; % No. of decision variables
NIND = 100; % No. of individuals
GGAP = 0.9; % Generation gap
XOV = 0.7; % Crossover rate
MUTR = 0.0175; % Mutation rate
MAXGEN = 100; % No. of generations

FieldD = [LIND*ones(1,NVAR); Vs_profile_lower; Vs_profile_upper; ...
    ones(1,NVAR); zeros(1,NVAR); zeros(1,NVAR); zeros(1,NVAR)];


time_index = zeros(MAXGEN,1);
Obj_index = zeros(MAXGEN,1);

tic
disp('Performing inversion ...')
% Initialise population
Chrom = crtbp(NIND, LIND*NVAR); % Create binary population
variable = bs2rv(Chrom, FieldD);
% Evaluate objective fn.
ObjV = GAobjfun(variable, layers_num, dispersions_R_true, Vp, den, ...
    f,modes_num_vec,index_vec_all); 
Gen = 0; % Counter
% Begin generational loop
while Gen < MAXGEN + 1
%     Gen
    % Assign fitness values to entire population
    FitnV = ranking(ObjV);
    % Select individuals for breeding
    SelCh = select('sus', Chrom, FitnV, GGAP);
    % Recombine individuals (crossover)
    SelCh = recombin('xovsp', SelCh, XOV);
    % Apply mutation
    SelCh = mut(SelCh, MUTR);
    % Evaluate offspring, call objective function
    ObjVSel = GAobjfun(bs2rv(SelCh, FieldD), layers_num, ...
        dispersions_R_true, Vp, den, f,modes_num_vec,index_vec_all);
    % Reinsert offspring into population
    [Chrom ObjV]=reins(Chrom, SelCh, 1, 1, ObjV, ObjVSel);
    % time and accuracy
    [Y I] = min(ObjV);
    variable = bs2rv(Chrom, FieldD);
    optimal_value = variable(I,:);
    Obj_index(Gen+1) = Y;
    time_index(Gen+1) = toc;
    % Increment counter
    Gen = Gen+1;
end
all_index = [time_index Obj_index];
disp('Inversion finished ...')

Phen = bs2rv(Chrom, FieldD);

[Y I] = min(ObjV);
variable = bs2rv(Chrom, FieldD);
optimal_value = variable(I,:);

Y_hat = optimal_value;

toc

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


%% remove relative path
rmpath (genpath('GA_ToolBox')) % remove relative path