%%%% Mulitlayer network analysis of musicians and nonmusicians. This script
%%%% uses SVD normalisation before performing group comparisons 
%%%% This is the code for reproducing the results from our recent preprint https://www.biorxiv.org/content/10.1101/2024.01.02.573886v1 where we perform a multilayer network analysis of resting-state MEG data of musicians and non-musicians.
clear all
close all
clc
%% KM's paths
% addpath(genpath('/Users/lpxknm/Documents/OMEGA/')) %%%% m_AEC.mat = musicians; nm_AEC.mat = nonmusicians
% addpath(genpath('/Users/lpxknm/Documents/MATLAB/my_matlab/'))
addpath(genpath('/Users/kanad/Documents/MATLAB/'))
%% PT's paths
% addpath(genpath('G:\linux\matlab\multilayer_compare'))

%% JM's paths
% addpath(genpath('/Users/jilmeier/Documents/WORK/multilayer_project/codes_JM_021117/AEC_data')) %%%% m_AEC.mat = musicians; nm_AEC.mat = nonmusicians
% cd('/Users/jilmeier/Documents/WORK/multilayer_project/codes_JM_021117/AEC_data') %%%% m_AEC.mat = musicians; nm_AEC.mat = nonmusicians
% 
% addpath('/Users/jilmeier/Documents/WORK/ALS_data_Matlab/BCT/2017_01_15_BCT')
     
% addpath('/Users/jilmeier/Documents/WORK/multilayer_project/codes_JM_021117/GraphVar-1.03/src/ext/GenLouvain2.0/privat')

%% Get data
tic
m_AEC = load('/Users/kanad/Documents/MATLAB/OMEGA_March19/m_AEC_231117.mat') % musicians
nm_AEC = load('/Users/kanad/Documents/MATLAB/OMEGA_March19/nm_AEC_231117.mat') % nonmusicians
% prof_mu = [2:7,10,11,14:27];
% bands_low = [1.0000    4.0000    6.5000    8.5000   10.5000   12.5000
% 16.5000   20.5000];
% bands_high = [3.5000    6.0000    8.0000   10.0000   12.0000   16.0000
% 20.0000   28.0000];
% 

%% Data for individual participants
%%%%%%%%%%%%%%%%%%%%%%%% musicians %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_aec_t = squeeze(m_AEC.AEC(:,:,3,:)); %%% 6.5-8Hz
m_aec_a = squeeze(m_AEC.AEC(:,:,4,:)); %%% 8.5-10Hz
m_aec_b = squeeze(m_AEC.AEC(:,:,5,:)); %%% 10.5-12Hz
m_aec_g = squeeze(m_AEC.AEC(:,:,6,:)); %%% 12.5-16Hz
m_aec_g2 = squeeze(m_AEC.AEC(:,:,7,:)); %%% 16.5-20Hz
%%%%%%%%%%%%%%%%%%%%%%% nonmusicians %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nm_aec_t = squeeze(nm_AEC.AEC(:,:,3,:)); %%% 6.5-8Hz
nm_aec_a = squeeze(nm_AEC.AEC(:,:,4,:)); %%% 8.5-10Hz
nm_aec_b = squeeze(nm_AEC.AEC(:,:,5,:)); %%% 10.5-12Hz
nm_aec_g = squeeze(nm_AEC.AEC(:,:,6,:)); %%% 12.5-16Hz
nm_aec_g2 = squeeze(nm_AEC.AEC(:,:,7,:)); %%% 16.5-20Hz
%% Layer specific assignment (duplex - alpha and beta band) and other parameters
%%%%%% Weight normalisation using the maximum link weight in the network so
%%%%%% that the values of the Adj matrix range from 0 1 (useful for single
%%%%%% layer analysis). The normalised Adj matrix is currently NOT used in
%%%%%% the multilayer analysis as SVD is used.
%%%%%% weight_conversion is part of the Brain Connectivity Toolbox
%%% 
FC_m = squeeze(mean(m_AEC.AEC,4)); %%%% 
FC_nm = squeeze(mean(nm_AEC.AEC,4));
N = size(FC_m,1); %% number of nodes
E = (N^2-N)/2; 
M = 2; %%% number of layers
% av = 0.1:0.05:1.1;
th = 0:0.05:pi/2; 
% Layer specific assigment and make sure the matrices are symmetric!
no_it = size(m_AEC.AEC,4); % Corresponds to number of participants

for i=1:no_it
    m_layer1(:,:,i)=abs(zeros(N,N)+tril(m_aec_t(:,:,i), -1)+tril(m_aec_t(:,:,i), -1)');
    m_layer2(:,:,i)=abs(zeros(N,N)+tril(m_aec_a(:,:,i), -1)+tril(m_aec_a(:,:,i), -1)');
    m_layer3(:,:,i)=abs(zeros(N,N)+tril(m_aec_b(:,:,i), -1)+tril(m_aec_b(:,:,i), -1)');
    m_layer4(:,:,i)=abs(zeros(N,N)+tril(m_aec_g(:,:,i), -1)+tril(m_aec_g(:,:,i), -1)');
    m_layer5(:,:,i)=abs(zeros(N,N)+tril(m_aec_g2(:,:,i), -1)+tril(m_aec_g2(:,:,i), -1)');
    
    nm_layer1(:,:,i)=abs(zeros(N,N)+tril(nm_aec_t(:,:,i), -1)+tril(nm_aec_t(:,:,i), -1)');
    nm_layer2(:,:,i)=abs(zeros(N,N)+tril(nm_aec_a(:,:,i), -1)+tril(nm_aec_a(:,:,i), -1)');
    nm_layer3(:,:,i)=abs(zeros(N,N)+tril(nm_aec_b(:,:,i), -1)+tril(nm_aec_b(:,:,i), -1)');
    nm_layer4(:,:,i)=abs(zeros(N,N)+tril(nm_aec_g(:,:,i), -1)+tril(nm_aec_g(:,:,i), -1)');
    nm_layer5(:,:,i)=abs(zeros(N,N)+tril(nm_aec_g2(:,:,i), -1)+tril(nm_aec_g2(:,:,i), -1)');
end
%%

avg_m_layer1 = mean(m_layer1,3);
avg_m_layer2 = mean(m_layer2,3);
avg_m_layer3 = mean(m_layer3,3);
avg_m_layer4 = mean(m_layer4,3);
avg_m_layer5 = mean(m_layer5,3);
%%%%%%%%%%%%%%%%%%%%%
avg_nm_layer1 = mean(nm_layer1,3);
avg_nm_layer2 = mean(nm_layer2,3);
avg_nm_layer3 = mean(nm_layer3,3);
avg_nm_layer4 = mean(nm_layer4,3);
avg_nm_layer5 = mean(nm_layer5,3);
%% Pre-allocate - Graph measures of interest
% 1) degree correlation
dc_m = zeros(1,no_it);
dc_nm = zeros(1,no_it);
% 2) participation coefficent
pc_m = zeros(N,no_it);
pc_nm = zeros(N,no_it);
% 3) community structure
Qii_m = zeros(1,no_it);
Qii_nm = zeros(1,no_it);

% mat_av_m = zeros(N,N,2);
% mat_av_nm = zeros(N,N,2);
mat_av_norm_m = zeros(N,N,2);
mat_av_norm_nm = zeros(N,N,2);
X = zeros(N,numel(th),numel(th));

%%
for it = 1:no_it
    
%%%%%%%%%%%%%%%%%%%%%%%%%% change average connectivity %%%%%%%%%%%%%%%%%%%%
%%%% Musicians - Check layer assignment carefully
    mat_av_m(:,:,1) = m_layer1(:,:,it);
    mat_av_m(:,:,2) = m_layer2(:,:,it);
    mat_av_m(:,:,3) = m_layer3(:,:,it);
    mat_av_m(:,:,4) = m_layer4(:,:,it);
    mat_av_m(:,:,5) = m_layer5(:,:,it);
%%%% Nonmusicians - Check layer assignment carefully  
    mat_av_nm(:,:,1) = nm_layer1(:,:,it);
    mat_av_nm(:,:,2) = nm_layer2(:,:,it);
    mat_av_nm(:,:,3) = nm_layer3(:,:,it);
    mat_av_nm(:,:,4) = nm_layer4(:,:,it);
    mat_av_nm(:,:,5) = nm_layer5(:,:,it);
%%%%%%%%%%%%%%%%%%%%%%%% SVD normalisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Musicians    
%     c_m = mean(mean(mean(mat_av_m))); % Mean of the whole multilayer,
%     including diagonal
% Calculate intra-layer average connectivity per layer, write all entries
% in a vector, except for diagonal
    lv_av_m=zeros(2*E,2);
    u=1;
    for i=1:N
        for j=1:N
            if i~=j
                lv_av_m(u,1)=mat_av_m(i,j,1);
                lv_av_m(u,2)=mat_av_m(i,j,2);
                lv_av_m(u,3)=mat_av_m(i,j,3);
                lv_av_m(u,4)=mat_av_m(i,j,4);
                lv_av_m(u,5)=mat_av_m(i,j,5);
                u=u+1;
            end
        end
    end
    % lv is a vector with all link weights of layer one in its first column and all link weights of layer 2 in its second column
    cm_1=mean(lv_av_m(:,1)); %average connectivity of layer one
    cm_2=mean(lv_av_m(:,2)); %average connectivity of layer two 
    cm_3=mean(lv_av_m(:,3));
    cm_4=mean(lv_av_m(:,4));
    cm_5=mean(lv_av_m(:,5));
    cm_new=mean([cm_1, cm_2, cm_3, cm_4 , cm_5]);
% cm_new = 1;
%%%%%%%%%%%%%%%%%%%% Five layer block adjacency matrix %%%%%%%%%%%%%%%%%%%%
   Adj_m = [mat_av_m(:,:,1) eye(N)*cm_new eye(N)*cm_new eye(N)*cm_new eye(N)*cm_new;...
           eye(N)*cm_new mat_av_m(:,:,2) eye(N)*cm_new  eye(N)*cm_new eye(N)*cm_new;... 
           eye(N)*cm_new eye(N)*cm_new mat_av_m(:,:,3) eye(N)*cm_new eye(N)*cm_new;...
           eye(N)*cm_new eye(N)*cm_new eye(N)*cm_new mat_av_m(:,:,4) eye(N)*cm_new;...
           eye(N)*cm_new eye(N)*cm_new eye(N)*cm_new eye(N)*cm_new mat_av_m(:,:,5)];
%%%%%%%%%%%%%%%%%%%% Four layer block adjacency matrix %%%%%%%%%%%%%%%%%%%%
%    Adj_m = [mat_av_m(:,:,1) eye(N)*cm_new eye(N)*cm_new eye(N)*cm_new;...
%            eye(N)*cm_new mat_av_m(:,:,2) eye(N)*cm_new  eye(N)*cm_new;... 
%            eye(N)*cm_new eye(N)*cm_new mat_av_m(:,:,3) eye(N)*cm_new;...
%            eye(N)*cm_new eye(N)*cm_new eye(N)*cm_new mat_av_m(:,:,4)];
%%%%%%%%%%%%%%%%%%%% 3 layer block adjacency matrix %%%%%%%%%%%%%%%%%%%%%%%
%    Adj_m = [mat_av_m(:,:,1) eye(N)*cm_new eye(N)*cm_new;...
%            eye(N)*cm_new mat_av_m(:,:,2) eye(N)*cm_new;... 
%            eye(N)*cm_new eye(N)*cm_new mat_av_m(:,:,3)];
%%%%%%%%%%%%%%%%%%% 2 layer block adjacency matrix %%%%%%%%%%%%%%%%%%%%%%%%
%    Adj_m = [mat_av_m(:,:,1) eye(N)*cm_new;...
%            eye(N)*cm_new mat_av_m(:,:,2)];

    Adj_norm_m = eignorm(Adj_m); %%%% SVD normalisation
    % Adj_norm_m is not completely symmetric. It needs to be made
    % symmetric.
    Adj_norm_m=zeros(5*N,5*N)+tril(Adj_norm_m, -1)+tril(Adj_norm_m, -1)';  % 78*3 
    
    mat_av_norm_m(:,:,1) = Adj_norm_m(1:N,1:N);
    mat_av_norm_m(:,:,2) = Adj_norm_m(N+1:N*2,N+1:N*2);
    mat_av_norm_m(:,:,3) = Adj_norm_m(N*2+1:N*3,N*2+1:N*3);
    mat_av_norm_m(:,:,4) = Adj_norm_m(N*3+1:N*4,N*3+1:N*4);
    mat_av_norm_m(:,:,5) = Adj_norm_m(N*4+1:N*5,N*4+1:N*5);
%%%%%%%%%%%%%%% Save individual data for averaging %%%%%%%%%%%%%%%%%%%%%%%%
    svd_av_norm_m(:,:,1,it) = Adj_norm_m(1:N,1:N);
    svd_av_norm_m(:,:,2,it) = Adj_norm_m(N+1:N*2,N+1:N*2);
    svd_av_norm_m(:,:,3,it) = Adj_norm_m(N*2+1:N*3,N*2+1:N*3);
    svd_av_norm_m(:,:,4,it) = Adj_norm_m(N*3+1:N*4,N*3+1:N*4);
    svd_av_norm_m(:,:,5,it) = Adj_norm_m(N*4+1:N*5,N*4+1:N*5);
%%%%%%%%%%%%%%%%%%%%%%%%%% nonmusicians   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% cnm_new = mean(mean(mean(mat_av_nm)));
lv_av_nm=zeros(2*E,2);
    u=1;
    for i=1:N
        for j=1:N
            if i~=j
                lv_av_nm(u,1)=mat_av_nm(i,j,1);
                lv_av_nm(u,2)=mat_av_nm(i,j,2);
                lv_av_nm(u,3)=mat_av_nm(i,j,3);
                lv_av_nm(u,4)=mat_av_nm(i,j,4);
                lv_av_nm(u,5)=mat_av_nm(i,j,5);
                u=u+1;
            end
        end
    end
% lv is a vector with all link weights of layer one in its first column and all link weights of layer 2 in its second column
    cnm_1=mean(lv_av_nm(:,1)); %average connectivity of layer one
    cnm_2=mean(lv_av_nm(:,2)); %average connectivity of layer two
    cnm_3=mean(lv_av_nm(:,3)); 
    cnm_4=mean(lv_av_nm(:,4));
    cnm_5=mean(lv_av_nm(:,5));
    cnm_new=mean([cnm_1, cnm_2, cnm_3, cnm_4, cnm_5]);
% cnm_new = 1 ;
%%%%%%%%%%%%%%%%%%%% Five layer block adjacency matrix %%%%%%%%%%%%%%%%%%%%
    Adj_nm = [mat_av_nm(:,:,1) eye(N)*cnm_new eye(N)*cnm_new eye(N)*cnm_new eye(N)*cnm_new;...
             eye(N)*cnm_new mat_av_nm(:,:,2) eye(N)*cnm_new eye(N)*cnm_new eye(N)*cnm_new;...
             eye(N)*cnm_new eye(N)*cnm_new mat_av_nm(:,:,3) eye(N)*cnm_new eye(N)*cnm_new;...
             eye(N)*cnm_new eye(N)*cnm_new eye(N)*cnm_new mat_av_nm(:,:,4) eye(N)*cnm_new;...
             eye(N)*cnm_new eye(N)*cnm_new eye(N)*cnm_new eye(N)*cnm_new mat_av_nm(:,:,5)];
%%%%%%%%%%%%%%%%%%%% Four layer block adjacency matrix %%%%%%%%%%%%%%%%%%%%
%     Adj_nm = [mat_av_nm(:,:,1) eye(N)*cnm_new eye(N)*cnm_new eye(N)*cnm_new;...
%              eye(N)*cnm_new mat_av_nm(:,:,2) eye(N)*cnm_new eye(N)*cnm_new;...
%              eye(N)*cnm_new eye(N)*cnm_new mat_av_nm(:,:,3) eye(N)*cnm_new;...
%              eye(N)*cnm_new eye(N)*cnm_new eye(N)*cnm_new mat_av_nm(:,:,4)];
%%%%%%%%%%%%%%%%% 3 layer block adjacency matrix  %%%%%%%%%%%%%%%%%%%%%%%%%
%          Adj_nm = [mat_av_nm(:,:,1) eye(N)*cnm_new eye(N)*cnm_new;...
%              eye(N)*cnm_new mat_av_nm(:,:,2) eye(N)*cnm_new;...
%              eye(N)*cnm_new eye(N)*cnm_new mat_av_nm(:,:,3)];
%%%%%%%%%%%%%%%%%%% 2 layer block adjacency matrix %%%%%%%%%%%%%%%%%%%%%%%%
%          Adj_nm = [mat_av_nm(:,:,1) eye(N)*cnm_new;...
%              eye(N)*cnm_new mat_av_nm(:,:,2)];
         
    Adj_norm_nm = eignorm(Adj_nm);
    Adj_norm_nm=zeros(5*N,5*N)+tril(Adj_norm_nm, -1)+tril(Adj_norm_nm, -1)'; % 78*3 
    
    mat_av_norm_nm(:,:,1) = Adj_norm_nm(1:N,1:N);
    mat_av_norm_nm(:,:,2) = Adj_norm_nm(N+1:N*2,N+1:N*2);
    mat_av_norm_nm(:,:,3) = Adj_norm_nm(N*2+1:N*3,N*2+1:N*3);
    mat_av_norm_nm(:,:,4) = Adj_norm_nm(N*3+1:N*4,N*3+1:N*4);
    mat_av_norm_nm(:,:,5) = Adj_norm_nm(N*4+1:N*5,N*4+1:N*5);
%%%%%%%%%%%%%%%%% Save individual participants data for averaging %%%%%%%%%
    svd_av_norm_nm(:,:,1,it) = Adj_norm_nm(1:N,1:N);
    svd_av_norm_nm(:,:,2,it) = Adj_norm_nm(N+1:N*2,N+1:N*2);
    svd_av_norm_nm(:,:,3,it) = Adj_norm_nm(N*2+1:N*3,N*2+1:N*3);
    svd_av_norm_nm(:,:,4,it) = Adj_norm_nm(N*3+1:N*4,N*3+1:N*4);
    svd_av_norm_nm(:,:,5,it) = Adj_norm_nm(N*4+1:N*5,N*4+1:N*5);
% %%%%%%%%%%%%%%%%%%%%% degree correlation - musicians %%%%%%%%%%%%%%%%%%%%%%
%     dc_t_m = degreecorr(mat_av_norm_m);
% %     dc_m(it) = dc_t_m(1,2);
%     dc_m_12(it) = dc_t_m(1,2);
%     dc_m_13(it) = dc_t_m(1,3);
%     dc_m_14(it) = dc_t_m(1,4);
%     dc_m_15(it) = dc_t_m(1,5);
%     
%     dc_m_23(it) = dc_t_m(2,3);
%     dc_m_24(it) = dc_t_m(2,4);
%     dc_m_25(it) = dc_t_m(2,5);
%     
%     dc_m_34(it) = dc_t_m(3,4);
%     dc_m_35(it) = dc_t_m(3,5);
%     
%     dc_m_45(it) = dc_t_m(4,5);
% %%%%%%%%%%%%%%%%%%%% degree correlation - nonmusicians %%%%%%%%%%%%%%%%%%%%
%     dc_t_nm = degreecorr(mat_av_norm_nm);
% %     dc_nm(it) = dc_t_nm(1,2);
%     dc_nm_12(it) = dc_t_nm(1,2);
%     dc_nm_13(it) = dc_t_nm(1,3);
%     dc_nm_14(it) = dc_t_nm(1,4);
%     dc_nm_15(it) = dc_t_nm(1,5);
%     
%     dc_nm_23(it) = dc_t_nm(2,3);
%     dc_nm_24(it) = dc_t_nm(2,4);
%     dc_nm_25(it) = dc_t_nm(2,5);
%     
%     dc_nm_34(it) = dc_t_nm(3,4);
%     dc_nm_35(it) = dc_t_nm(3,5);
%     
%     dc_nm_45(it) = dc_t_nm(4,5);
%%%%%%%%%%%%%%%%% Layer participation coefficient %%%%%%%%%%%%%%%%%%%%%%%%%    
%      pc_m(:,it) = partc_layers(mat_av_norm_m); 
%%%%%%%%%%%%%%%%%%%%%%% Entropy - musicians %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     entrop_m(:,it) = entropylayer(mat_av_norm_m);
%%%%%%%%%%%%%%%%%%%% Layer participation coefficient %%%%%%%%%%%%%%%%%%%%%%    
%      pc_nm(:,it) = partc_layers(mat_av_norm_nm); 
%%%%%%%%%%%%%%%%%%%%%%% Entropy - nonmusicians %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     entrop_nm(:,it) = entropylayer(mat_av_norm_nm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Efficiency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             for i = 1:size(mat_av_norm_m,3) % Musicians % Ignore this for now
%                 Eglob_m(i,it) = efficiency_wei(mat_av_norm_m(:,:,i));
%                 Eloc_m(:,it) = efficiency_wei(mat_av_norm_m(:,:,i),1);
%                 Eglob_nm(i,it) = efficiency_wei(mat_av_norm_nm(:,:,i));
%                 Eloc_nm(:,it) = efficiency_wei(mat_av_norm_nm(:,:,i),1);
%             end
%%%%%%%%%%%%%%%%%%%%%%%% Community - musicians %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A_m{1} = mat_av_norm_m(:,:,1);
    A_m{2} = mat_av_norm_m(:,:,2);
    A_m{3} = mat_av_norm_m(:,:,3);
    A_m{4} = mat_av_norm_m(:,:,4);
    A_m{5} = mat_av_norm_m(:,:,5);
    
    gamma = 1;
    omega = cm_new;
    
    T_m=length(A_m);
    B_m=spalloc(N*T_m,N*T_m,N*N*T_m+2*N*T_m); %% 
    twomu_m=0;
    for s_m=1:T_m % For each layer
        k_m=sum(A_m{s_m});
        twom_m=sum(k_m);    % Sums up all the matric elements of one layer
        twomu_m=twomu_m+twom_m; % Add this to twomu_m
        indx_m=(1:N)+(s_m-1)*N; 
        B_m(indx_m,indx_m)=A_m{s_m}-gamma*k_m'*k_m/twom_m;
    end
%     twomu_m=twomu_m+2*omega*N*(T_m-1);
     twomu_m=twomu_m+omega*N*T_m*(T_m-1);
    B_m = B_m + omega*spdiags(ones(N*T_m,2),[-N,N],N*T_m,N*T_m);
    [D_m,Q_m] = genlouvain(B_m);
    Qii_m(:,it) = Q_m/twomu_m; %%%%%%%%%%%%%%%%%%
    Dii_m{it}(:,:) = reshape(D_m,N,T_m); % Save nodal assignment per participant into a cell
%%%%%%%%%%%%%%%%%%%%%%%% Community - nonmusicians %%%%%%%%%%%%%%%%%%%%%%%%%   
    A_nm{1} = mat_av_norm_nm(:,:,1);
    A_nm{2} = mat_av_norm_nm(:,:,2);
    A_nm{3} = mat_av_norm_nm(:,:,3);
    A_nm{4} = mat_av_norm_nm(:,:,4);
    A_nm{5} = mat_av_norm_nm(:,:,5);
    
%     gamma = 1;
      omega = cnm_new;

    T_nm=length(A_nm);
    B_nm=spalloc(N*T_nm,N*T_nm,N*N*T_nm+2*N*T_nm);
    twomu_nm=0;
    for s_nm=1:T_nm
        k_nm=sum(A_nm{s_nm});
        twom_nm=sum(k_nm);
        twomu_nm=twomu_nm+twom_nm;
        indx_nm=(1:N)+(s_nm-1)*N; % Bug fixed!
        B_nm(indx_nm,indx_nm)=A_nm{s_nm}-gamma*k_nm'*k_nm/twom_nm;
    end
%     twomu_nm=twomu_nm+2*omega*N*(T_nm-1);
    twomu_nm=twomu_nm+omega*N*(T_nm-1);
    B_nm = B_nm + omega*spdiags(ones(N*T_nm,2),[-N,N],N*T_nm,N*T_nm);
    [D_nm,Q_nm] = genlouvain(B_nm);
    Qii_nm(:,it) = Q_nm/twomu_nm; %%%%%%%%%%%%%%%%%%
    Dii_nm{it}(:,:) = reshape(D_nm,N,T_nm);%%%% Save nodal assignment per participant. 
    %%%%% Not sure how to test Dii between groups. 

    disp(sprintf('Iteration: %f \n',it))
    toc   
   
end
%% Statistics 
% pc_nm = squeeze(pc_nm);
% pc_m = squeeze(pc_m);

%%%%% test pc difference
%%%%%% pc_m and pc_nm are 78 x 31. Perform 78 tests and correct using FDR
% for m=1:78
%     [p_pc(m) n_h_pc(m)] = ranksum(pc_m(m,:),pc_nm(m,:)); %%%% 
%     [p_entrop(m) n_h_entrop(m)] = ranksum(entrop_m(m,:),entrop_nm(m,:));
% end
% [h_pc, crit_p_pc, adj_ci_cvrg_pc, adj_p_pc] = fdr_bh(p_pc,0.05,'pdep','yes') 
% [h_entrop, crit_p_entrop, adj_ci_cvrg_entrop, adj_p_entrop] = fdr_bh(p_entrop,0.05,'pdep','yes') 
% % Some plots for descriptive analysis
% avg_pc_m = mean(pc_m,2);
% avg_pc_nm = mean(pc_nm,2);
% avg_pc = [avg_pc_m,avg_pc_nm];

% figure
% boxplot(avg_pc)
% violin(avg_pc)

% figure
% plot(avg_pc_m)
% hold on
% plot((avg_pc_nm),'r')

%% Qii
%%%%% test community structrue difference
%%%% Qii_m and Qii_nm are 1x31 vectors, perform one test
[P_Qii H_Qii] = ranksum(Qii_m,Qii_nm);


Qii_temp = [Qii_m',Qii_nm'];%%%%% boxplot(temp) will plot two boxes, one per colum size(temp) = 31 x 2
figure
violin(temp)

figure
h = boxplot(temp);
set(h,{'linew'},{4});

%%%%%%%%%%%%%%%%%%%%%%%%%%% histogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Qii distribution - histogram
% figure
% h = histfit(Qii_m,50)
% set(h(1),'facecolor','g'); set(h(1),'facealpha',1); set(gca,'YLim',[0 10])
% hold on
% h2 = histfit(Qii_nm,50)
% set(h2(1),'facecolor','b');set(h2(1),'facealpha',0.50); set(gca,'YLim',[0 10])
% axis tight
% title('Community structure in Musicians and Nonmusicians')
%% DC
%%%% dc_m and dc_nm are 1 x 31 vectors, perform one test
% [P_dc H_dc] = ranksum(dc_m,dc_nm);

[P_dc12 H_dc12] = ranksum(dc_m_12,dc_nm_12);
[P_dc13 H_dc13] = ranksum(dc_m_13,dc_nm_13);
[P_dc14 H_dc14] = ranksum(dc_m_14,dc_nm_14); % Significant
[P_dc15 H_dc15] = ranksum(dc_m_15,dc_nm_15);

[P_dc23 H_dc23] = ranksum(dc_m_23,dc_nm_23);
[P_dc24 H_dc24] = ranksum(dc_m_24,dc_nm_24); % Significant
[P_dc25 H_dc25] = ranksum(dc_m_25,dc_nm_25);

[P_dc34 H_dc34] = ranksum(dc_m_34,dc_nm_34);
[P_dc35 H_dc35] = ranksum(dc_m_35,dc_nm_35);

[P_dc45 H_dc45] = ranksum(dc_m_45,dc_nm_45);

[h_dc, crit_p_dc, adj_ci_cvrg_dc, adj_p_dc]=fdr_bh(P_dc14,q,'pdep','yes')

all_dc_p = [P_dc12, P_dc13, P_dc14, P_dc15, P_dc23, P_dc24, P_dc25, P_dc34, P_dc35, P_dc45]; 
[h_dc, crit_p_dc, adj_ci_cvrg_dc, adj_p_dc]=fdr_bh(all_dc_p,0.05,'pdep','yes')
%%

tiledlayout(1,2);
nexttile
notBoxPlot(dc14,'jitter',0.5,'style','sdline')
xticklabels({'Musicians','Non-musicians'})
ylabel('Degree correlations')
title('DC between layers 1&4')

box
nexttile
notBoxPlot(dc24,'jitter',0.5,'style','sdline')
xticklabels({'Musicians','Non-musicians'})
ylabel('Degree correlations')
title('DC between layers 2&4')
box


figure
notBoxPlot(Qii_temp,'jitter',0.5,'style','sdline')
box
xticklabels({'Musicians','Non-musicians'})
ylabel('Modularity (Q)')
title('Modularity between groups')

%% Global efficiency
% for l = 1:size(Eglob_m,2)
% [P_Eglob(l), H_Eglob(l)] = ranksum(Eglob_m(:,l),Eglob_nm(:,l));    
% end
% [h_Eglob, crit_p_Eglob, adj_ci_cvrg_Eglob, adj_p_Eglob] = fdr_bh(P_Eglob,0.05,'pdep','yes')
% [P_Eglob, H_Eglob] = ranksum(Eglob_m,Eglob_nm);
% Eglob = [Eglob_m',Eglob_nm'];
% violin(Eglob)

%%%%% local efficiency
% for j = 1:78
%     [P_Eloc(j) n_H_Eloc(j)] = ranksum(Eloc_m(j,:),Eloc_nm(j,:));
% end
% [h_Eloc, crit_p_Eloc, adj_ci_cvrg_Eloc, adj_p_Eloc] = fdr_bh(P_Eloc,0.05,'pdep','yes')


%% NBS analysis - single layer - non singificant
% m_svd1 = squeeze(svd_av_norm_m(:,:,1,:));
% avg_m_svd1 = mean(m_svd1,3);
% nm_svd1 = squeeze(svd_av_norm_nm(:,:,1,:));
% avg_nm_svd1 = mean(nm_svd1,3);
% 
% m_svd2 = squeeze(svd_av_norm_m(:,:,2,:));
% avg_m_svd2 = mean(m_svd2,3);
% nm_svd2 = squeeze(svd_av_norm_nm(:,:,2,:));
% avg_nm_svd2 = mean(nm_svd2,3);
% 
% m_svd3 = squeeze(svd_av_norm_m(:,:,3,:));
% avg_m_svd3 = mean(m_svd3,3);
% nm_svd3 = squeeze(svd_av_norm_nm(:,:,3,:));
% avg_nm_svd3 = mean(nm_svd3,3);
% 
% m_svd4 = squeeze(svd_av_norm_m(:,:,4,:));
% avg_m_svd4 = mean(m_svd4,3);
% nm_svd4 = squeeze(svd_av_norm_nm(:,:,4,:));
% avg_nm_svd4 = mean(nm_svd4,3);
% 
% m_svd5 = squeeze(svd_av_norm_m(:,:,5,:));
% avg_m_svd5 = mean(m_svd5,3);
% nm_svd5 = squeeze(svd_av_norm_nm(:,:,5,:));
% avg_nm_svd5 = mean(nm_svd5,3);
%% Correlation
% %%%%%% 6.5-8Hz
% t_m_svd1 = squeeze(mean(mean(m_svd1)));
% %%%%% 8.5-10Hz
% a1_m_svd2 = squeeze(mean(mean(m_svd2))); %%%% 0.1
% %%%%% 10.5-12Hz
% a2_m_svd3 = squeeze(mean(mean(m_svd3)));
% %%%% 12.5-16Hz
% b1_m_svd4 = squeeze(mean(mean(m_svd4)));
% %%%% 16.5-20Hz
% b2_m_svd5 = squeeze(mean(mean(m_svd5)));

%%
% figure 
% subplot(321)
% imagesc(avg_m_svd1)
% colormap jet
% 
% subplot(322)
% imagesc(avg_m_svd2)
% colormap jet
% 
% subplot(323)
% imagesc(avg_m_svd3)
% colormap jet
% 
% 
% subplot(324)
% imagesc(avg_m_svd4)
% colormap jet
% 
% subplot(325)
% imagesc(avg_m_svd5)
% colormap jet


%%
yrs_exp = [9
25
20
20
10
9
10
12
10
20
18
19
12
14
20
10
18
8
11
16
15
16
10
10
16
12
15
5
5
6
10]; 

figure
scatter(yrs_exp, Qii_m)
lsline

% ltheta = cat(3,m_svd1,nm_svd1); % 4-6Hz
% htheta = cat(3,m_svd2,nm_svd2);  % 6.5-8Hz
% lalpha = cat(3,m_svd3,nm_svd3); % 8.5-10Hz
% 
% raw_ltheta = cat(3,m_layer1,nm_layer1);
% raw_htheta = cat(3,m_layer2,nm_layer2);
% raw_alpha = cat(3,m_layer3,nm_layer3);


% %% Single layer community 
% gamma1 = 1;
% for o = 1:size(m_svd1,3)
% [M_m1(:,o), Q_m1(:,o)] = community_louvain(m_svd1(:,:,o),gamma1,[],'modularity');
% [M_m2(:,o), Q_m2(:,o)] = community_louvain(m_svd2(:,:,o),gamma1,[],'modularity');
% [M_m3(:,o), Q_m3(:,o)] = community_louvain(m_svd3(:,:,o),gamma1,[],'modularity');
% [M_m4(:,o), Q_m4(:,o)] = community_louvain(m_svd4(:,:,o),gamma1,[],'modularity');
% [M_m5(:,o), Q_m5(:,o)] = community_louvain(m_svd5(:,:,o),gamma1,[],'modularity');
% end
% 
% for q = 1:size(nm_svd1,3)
% [M_nm1(:,q), Q_nm1(:,q)] = community_louvain(nm_svd1(:,:,q),gamma1,[],'modularity');
% [M_nm2(:,q), Q_nm2(:,q)] = community_louvain(nm_svd2(:,:,q),gamma1,[],'modularity');
% [M_nm3(:,q), Q_nm3(:,q)] = community_louvain(nm_svd3(:,:,q),gamma1,[],'modularity');
% [M_nm4(:,q), Q_nm4(:,q)] = community_louvain(nm_svd4(:,:,q),gamma1,[],'modularity');
% [M_nm5(:,q), Q_nm5(:,q)] = community_louvain(nm_svd5(:,:,q),gamma1,[],'modularity');
% end

%% statistical testing - 
% [P_l1, H_l1] = ranksum(Q_m1,Q_nm1); % Significant
% [P_l2, H_l2] = ranksum(Q_m2,Q_nm2); % 
% [P_l3, H_l3] = ranksum(Q_m3,Q_nm3);
% [P_l4, H_l4] = ranksum(Q_m4,Q_nm4);
% [P_l5, H_l5] = ranksum(Q_m5,Q_nm5);

% Q_l1 = [Q_m1',Q_nm1'];
% Q_l2 = [Q_m2',Q_nm2'];
% Q_l3 = [Q_m3',Q_nm3'];

% figure
% subplot(131)
% boxplot(Q_l1)
% % violin(Q_l1)
% 
% subplot(132)
% boxplot(Q_l2)
% % violin(Q_l2)
% 
% subplot(133)
% boxplot(Q_l3)
% % violin(Q_l3)
%%
%{ 
%%%%%%%%%%%%%%%%%%%%%%% Single layer analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Single layer efficiency - Global
Eglob_m1 = zeros(1,it); Eloc_m1 = zeros(N,it);
Eglob_m2 = zeros(1,it); Eloc_m2 = zeros(N,it);
Eglob_m3 = zeros(1,it); Eloc_m3 = zeros(N, it);
%%%%% Single layer efficiency - Local
Eglob_nm1 = zeros(1,it); Eloc_nm1 = zeros(N,it);
Eglob_nm2 = zeros(1,it); Eloc_nm2 = zeros(N,it);
Eglob_nm3 = zeros(1,it); Eloc_nm3 = zeros(N, it);

tic
for t = 1:size(m_svd1,3)
   Eglob_m1(t) = efficiency_wei(m_svd1(:,:,t));
   Eloc_m1(:,t) = efficiency_wei(m_svd1(:,:,t),1); 
   Eglob_m2(t) = efficiency_wei(m_svd2(:,:,t));
   Eloc_m2(:,t) = efficiency_wei(m_svd2(:,:,t),1); 
   Eglob_m3(t) = efficiency_wei(m_svd3(:,:,t));
   Eloc_m3(:,t) = efficiency_wei(m_svd3(:,:,t),1); 
   Eglob_m4(t) = efficiency_wei(m_svd4(:,:,t));
   Eloc_m4(:,t) = efficiency_wei(m_svd4(:,:,t),1); 
   Eglob_m5(t) = efficiency_wei(m_svd5(:,:,t));
   Eloc_m5(:,t) = efficiency_wei(m_svd5(:,:,t),1); 
   disp(sprintf('Musician: %f \n',t))
   toc
end

tic
for u = 1:size(nm_svd1,3)
   Eglob_nm1(u) = efficiency_wei(nm_svd1(:,:,u));
   Eloc_nm1(:,u) = efficiency_wei(nm_svd1(:,:,u),1); 
   Eglob_nm2(u) = efficiency_wei(nm_svd2(:,:,u));
   Eloc_nm2(:,u) = efficiency_wei(nm_svd2(:,:,u),1); 
   Eglob_nm3(u) = efficiency_wei(nm_svd3(:,:,u));
   Eloc_nm3(:,u) = efficiency_wei(nm_svd3(:,:,u),1); 
   Eglob_nm4(u) = efficiency_wei(nm_svd4(:,:,u));
   Eloc_nm4(:,u) = efficiency_wei(nm_svd4(:,:,u),1); 
   Eglob_nm5(u) = efficiency_wei(nm_svd5(:,:,u));
   Eloc_nm5(:,u) = efficiency_wei(nm_svd5(:,:,u),1); 
   disp(sprintf('Non-musician: %f \n',u))
   toc
end

%%%%%%%%%%%%%% Single layer global efficiency - not significant %%%%%%%%%%%
[P_Eglobm1, H_Eglobm1] = ranksum(Eglob_m1,Eglob_nm1);[P_Eglobm2, H_Eglobm2] = ranksum(Eglob_m2,Eglob_nm2);
[P_Eglobm3, H_Eglobm3] = ranksum(Eglob_m3,Eglob_nm3);[P_Eglobm4, H_Eglobm4] = ranksum(Eglob_m4,Eglob_nm4);
[P_Eglobm5, H_Eglobm5] = ranksum(Eglob_m5,Eglob_nm5);
%%%%%%%%%%%%%%% Single layer local efficiency 
for l = 1:78
    [P_Elocl1(l) n_H_Elocl1(l)] = ranksum(Eloc_m1(l,:),Eloc_nm1(l,:));
    [P_Elocl2(l) n_H_Elocl2(l)] = ranksum(Eloc_m2(l,:),Eloc_nm2(l,:));
    [P_Elocl3(l) n_H_Elocl3(l)] = ranksum(Eloc_m3(l,:),Eloc_nm3(l,:));
    [P_Elocl4(l) n_H_Elocl4(l)] = ranksum(Eloc_m4(l,:),Eloc_nm4(l,:));
    [P_Elocl5(l) n_H_Elocl5(l)] = ranksum(Eloc_m5(l,:),Eloc_nm5(l,:));
end


[h_Eloc_l1, crit_p_l1, adj_ci_cvrg_l1, adj_p_l1] = fdr_bh(P_Elocl1,0.05,'pdep','yes')
[h_Eloc_l2, crit_p_l2, adj_ci_cvrg_l2, adj_p_l2] = fdr_bh(P_Elocl2,0.05,'pdep','yes')
[h_Eloc_l3, crit_p_l3, adj_ci_cvrg_l3, adj_p_l3] = fdr_bh(P_Elocl3,0.05,'pdep','yes')
[h_Eloc_l4, crit_p_l4, adj_ci_cvrg_l4, adj_p_l4] = fdr_bh(P_Elocl4,0.05,'pdep','yes')
[h_Eloc_l5, crit_p_l5, adj_ci_cvrg_l5, adj_p_l5] = fdr_bh(P_Elocl5,0.05,'pdep','yes')
%%
%}

savename = '/Users/kanad/Documents/MATLAB/OMEGA_March19/5layers_results_300319.mat'
save([savename])



dim=5*78; %%% layers X nodes
agree_m=zeros(dim,dim); % Pre allocate

%it should be a block adjacency matrix as big as the 3-layer one, so 3x78
for i=1:31
    Dii_m_new{i}=reshape(Dii_m{i},dim,1);
    Dii_nm_new{i}=reshape(Dii_nm{i},dim,1);
end


for AAL_reg1 = 1:dim % Number of AAL regions (need a symmetric matrix at the end)
    %agree_m_line = zeros(78,1); % For one layer
    for AAL_reg2=1:dim
        
        for i=1:31 % Number of participants

            %for lyr = 1:3 % Number of layers

          Comm_m1=Dii_m_new{1,i}(AAL_reg1); % Community in musicians
          Comm_m2=Dii_m_new{1,i}(AAL_reg2); 
          
          %test_m=bsxfun(@eq,Dii_m{1,i},Comm_m); % create test_m for every layer
          if Comm_m1==Comm_m2
              
            agree_m(AAL_reg1, AAL_reg2)=agree_m(AAL_reg1, AAL_reg2)+1;
         
          end
          
            %end

        end 
    
        %agree_m(AAL_reg,:) = agree_m_line; % Add them in an agreement matrix

    end
end




%collapse it


N=78;
agree_m_mat(:,:,1) = agree_m(1:N,1:N);
agree_m_mat(:,:,2) = agree_m(N+1:N*2,N+1:N*2);
agree_m_mat(:,:,3) = agree_m(N*2+1:N*3,N*2+1:N*3);
agree_m_mat(:,:,4) = agree_m(N*3+1:N*4,N*3+1:N*4);
agree_m_mat(:,:,5) = agree_m(N*4+1:N*5,N*4+1:N*5);
    
%now sum them up

agree_m_col=sum(agree_m_mat,3);
%put diagonal to zero  
for i=1:N
    agree_m_col(i,i)=0;
end

 
%now also make it proportional: 

%maximum is 31*5=155

agree_m_col=agree_m_col./155;

imagesc(agree_m_col)   

%agree_nm=agree_nm./31;

%and put the diagonal back to zero;

% for i=1:dim
%         agree_nm(i,i)=0;
% end
%imagesc(agree_nm)

%for no further confusion
agree_m=agree_m_col;



%now also make it proportional: 

% agree_m=agree_m./31;
% 
% %and put the diagonal back to zero;
% 
% for i=1:dim
%     
%         agree_m(i,i)=0;
%   
% end
% 
% 
% imagesc(agree_m)

%%%%%%%%%%%%%%%%%non musicians the same%%%%%%%%%%%%%%%%%%%%%%
dim=5*78;
agree_nm=zeros(dim,dim); % Pre allocate

%it should be a block adjacency matrix as big as the 3-layer one, so 3x78

for AAL_reg1 = 1:dim % Number of AAL regions (need a symmetric matrix at the end)
    for AAL_reg2=1:dim
        
        for i=1:31 % Number of participants

          Comm_m1=Dii_nm_new{1,i}(AAL_reg1); % Community in musicians
          Comm_m2=Dii_nm_new{1,i}(AAL_reg2); 
          if Comm_m1==Comm_m2
            agree_nm(AAL_reg1, AAL_reg2)=agree_nm(AAL_reg1, AAL_reg2)+1;
          end
        end 

    end
end

%collapse it


N=78;
agree_nm_mat(:,:,1) = agree_nm(1:N,1:N);
agree_nm_mat(:,:,2) = agree_nm(N+1:N*2,N+1:N*2);
agree_nm_mat(:,:,3) = agree_nm(N*2+1:N*3,N*2+1:N*3);
agree_nm_mat(:,:,4) = agree_nm(N*3+1:N*4,N*3+1:N*4);
agree_nm_mat(:,:,5) = agree_nm(N*4+1:N*5,N*4+1:N*5);
    
%now sum them up

agree_nm_col=sum(agree_nm_mat,3);
%put diagonal to zero  
for i=1:N
    agree_nm_col(i,i)=0;
end

 
%now also make it proportional: 

%maximum is 31*5=155

agree_nm_col=agree_nm_col./155;

imagesc(agree_nm_col)   

%agree_nm=agree_nm./31;

%and put the diagonal back to zero;

% for i=1:dim
%         agree_nm(i,i)=0;
% end
%imagesc(agree_nm)

%for no further confusion
agree_nm=agree_nm_col;

%%

figure (1)
colormap jet
subplot(121)
imagesc(agree_m)
colorbar
caxis([0 1])
subplot(122)
imagesc(agree_nm)
colorbar
caxis([0 1])


%run the genlouvain algorithm on these block-matrices
%treat it as 234 separate nodes
%we need the quality function B_m
%NN=78;
%%

dim=78;

gamma = 1;
    
    u_tot_m=sum(sum(agree_m))/2;
    B_m=zeros(dim,dim); %%   
    degrs=sum(agree_m,2);
    B_m=agree_m-gamma*(degrs*degrs')/(2*u_tot_m);

    [D_m_agree,Q_m_agree] = genlouvain(B_m);
    Q_m_agree_norm = Q_m_agree/(2*u_tot_m); %%%%%%%%%%%%%%%%%%

%now for non-musicians the same
gamma = 1;
    u_tot_nm=sum(sum(agree_nm))/2;
    
    B_nm=zeros(dim,dim); %% 
    
    %we need the degrees per node
    
    degrs_nm=sum(agree_nm,2);
    
    B_nm=agree_nm-gamma*(degrs_nm*degrs_nm')/(2*u_tot_nm);
    [D_nm_agree,Q_nm_agree] = genlouvain(B_nm);
    Q_nm_agree_norm = Q_nm_agree/(2*u_tot_nm); %%%%%%%%%%%%%%%%%%

% we have three versus 2 communities here




%%
% % % D_m_agree1 = D_m_agree(1:N,:); D_m_agree2 = D_m_agree(N+1:N*2,:); D_m_agree3 = D_m_agree(N*2+1:N*3,:); D_m_agree4 = D_m_agree(N*3+1:N*4,:);
% % % D_m_agree5 = D_m_agree(N*4+1:N*5,:);
% % % 
% % % D_nm_agree1 = D_nm_agree(1:N,:); D_nm_agree2 = D_nm_agree(N+1:N*2,:); D_nm_agree3 = D_nm_agree(N*2+1:N*3,:); D_nm_agree4 = D_nm_agree(N*3+1:N*4,:);
% % % D_nm_agree5 = D_nm_agree(N*4+1:N*5,:); 
% % % %%
% [zRand_1,SR_1,SAR_1,VI_1] = zrand(D_m_agree1,D_nm_agree1)
% 
% [zRand_2,SR_2,SAR_2,VI_2] = zrand(D_m_agree2,D_nm_agree2)
% 
% [zRand_3,SR_3,SAR_3,VI_3] = zrand(D_m_agree3,D_nm_agree3)
% 
% [zRand_4,SR_4,SAR_4,VI_4] = zrand(D_m_agree4,D_nm_agree4)
% 
% % [zRand_5,SR_5,SAR_5,VI_5] = zrand(D_m_agree5,D_nm_agree5)
% 
[zRand, SR, SAR, VI] = zrand(D_m_agree,D_nm_agree)
%zRand = 95.13; SR = 0.77; SAR = 0.54; VI = 0.95
%%
% close all
% 
% for mu = 1:size(m_svd1,3)
%    figure (mu)
%    colormap jet
%    imagesc(m_svd1(:,:,mu))
%    caxis([0 1])
% end
 
% % % 
% for nomu = 1:size(nm_svd5,3)
%    figure (nomu)
%    colormap jet
%    imagesc(nm_svd5(:,:,nomu))
%    caxis([0 1])
% end


