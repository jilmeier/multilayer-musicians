%%%% Multilayer network analysis of musicians and nonmusicians. This script
%%%% uses SVD normalisation before performing group comparisons 
clear all
close all
clc

%% JM's paths

% Add your own paths to GenLouvain-master and network_metrics folder

%addpath('/Users/jilmeier/Documents/musician/code/OMEGA_Feb24/GenLouvain-master')
%addpath('/Users/jilmeier/Documents/musician/code/OMEGA_Feb24/network_metrics')

%load('/Users/jilmeier/Documents/musician/code/OMEGA_Feb24/workspace_for_SVD_results.mat')


%% Get data

%m_AEC = load('/Users/kanad/Documents/MATLAB/OMEGA_March19/m_AEC_231117.mat') % musicians
%nm_AEC = load('/Users/kanad/Documents/MATLAB/OMEGA_March19/nm_AEC_231117.mat') % nonmusicians
% prof_mu = [2:7,10,11,14:27];
% bands_low = [1.0000    4.0000    6.5000    8.5000   10.5000   12.5000
% 16.5000   20.5000];
% bands_high = [3.5000    6.0000    8.0000   10.0000   12.0000   16.0000
% 20.0000   28.0000];
% 

%% Data for individual participants
%%%%%%%%%%%%%%%%%%%%%%% musicians %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_aec_t = squeeze(m_AEC.AEC(:,:,3,:)); %%% 6.5-8Hz
m_aec_a = squeeze(m_AEC.AEC(:,:,4,:)); %%% 8.5-10Hz
m_aec_b = squeeze(m_AEC.AEC(:,:,5,:)); %%% 10.5-12Hz
m_aec_g = squeeze(m_AEC.AEC(:,:,6,:)); %%% 12.5-16Hz
m_aec_g2 = squeeze(m_AEC.AEC(:,:,7,:)); %%% 16.5-20Hz
%%%%%%%%%%%%%%%%%%%%%% nonmusicians %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% community structure
Qii_m = zeros(1,no_it);
Qii_nm = zeros(1,no_it);


mat_av_norm_m = zeros(N,N,2);
mat_av_norm_nm = zeros(N,N,2);

omega_idx=0;
uiui=1;
nsamples=100;

Qii_m_it_100=[]; 
Dii_m_it_100=[];

%% analyze community structure over all individuals for all interlayer connectivity values, here omegas

for omega=0:0.01:1
    omega_idx=omega_idx+1;
    for it = 1:no_it %perform the following analysis for all 31x2 participants
        %     
        % %%%%%%%%%%%%%%%%%%%%%%%%%% change average connectivity %%%%%%%%%%%%%%%%%%%%
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
        % %%%%%%%%%%%%%%%%%%%%%%%% SVD normalisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            cm_new=omega;
        % cm_new = 1;
        % %%%%%%%%%%%%%%%%%%%% Five layer block adjacency matrix %%%%%%%%%%%%%%%%%%%%
           Adj_m = [mat_av_m(:,:,1) eye(N)*cm_new eye(N)*cm_new eye(N)*cm_new eye(N)*cm_new;...
                   eye(N)*cm_new mat_av_m(:,:,2) eye(N)*cm_new  eye(N)*cm_new eye(N)*cm_new;... 
                   eye(N)*cm_new eye(N)*cm_new mat_av_m(:,:,3) eye(N)*cm_new eye(N)*cm_new;...
                   eye(N)*cm_new eye(N)*cm_new eye(N)*cm_new mat_av_m(:,:,4) eye(N)*cm_new;...
                   eye(N)*cm_new eye(N)*cm_new eye(N)*cm_new eye(N)*cm_new mat_av_m(:,:,5)];

            Adj_norm_m = eignorm(Adj_m); %%%% SVD normalisation
            %if we want no SVD normalization: 
            % Adj_norm_m = Adj_m;
            
            % Adj_norm_m is not completely symmetric. It needs to be made
            % symmetric.
            Adj_norm_m=zeros(5*N,5*N)+tril(Adj_norm_m, -1)+tril(Adj_norm_m, -1)';  % 78*3 
            %extract the new interlayer connectivity
            inter_m=Adj_norm_m(1,79);
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
            cnm_new=omega;

        %%%%%%%%%%%%%%%%%%%% Five layer block adjacency matrix %%%%%%%%%%%%%%%%%%%%
            Adj_nm = [mat_av_nm(:,:,1) eye(N)*cnm_new eye(N)*cnm_new eye(N)*cnm_new eye(N)*cnm_new;...
                     eye(N)*cnm_new mat_av_nm(:,:,2) eye(N)*cnm_new eye(N)*cnm_new eye(N)*cnm_new;...
                     eye(N)*cnm_new eye(N)*cnm_new mat_av_nm(:,:,3) eye(N)*cnm_new eye(N)*cnm_new;...
                     eye(N)*cnm_new eye(N)*cnm_new eye(N)*cnm_new mat_av_nm(:,:,4) eye(N)*cnm_new;...
                     eye(N)*cnm_new eye(N)*cnm_new eye(N)*cnm_new eye(N)*cnm_new mat_av_nm(:,:,5)];

            Adj_norm_nm = eignorm(Adj_nm);
            %if we want no SVD normalization: 
            % Adj_norm_nm = Adj_nm;
            Adj_norm_nm=zeros(5*N,5*N)+tril(Adj_norm_nm, -1)+tril(Adj_norm_nm, -1)'; % 78*3 
            inter_nm=Adj_norm_nm(1,79);
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

        %%%%%%%%%%%%%%%%%%%%%%%% Community - musicians %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            A_m{1} = mat_av_norm_m(:,:,1);
            A_m{2} = mat_av_norm_m(:,:,2);
            A_m{3} = mat_av_norm_m(:,:,3);
            A_m{4} = mat_av_norm_m(:,:,4);
            A_m{5} = mat_av_norm_m(:,:,5);

            gamma = 1; %fixed value of gamma

            T_m=length(A_m); %5
            B_m=spalloc(N*T_m,N*T_m,N*N*T_m+2*N*T_m); %% 
            twomu_m=0;
            for s_m=1:T_m % For each layer
                k_m=sum(A_m{s_m}); %sums up weights per region
                twom_m=sum(k_m);    % Sums up all the matrix elements of one layer
                twomu_m=twomu_m+twom_m; % Add this to twomu_m
                indx_m=(1:N)+(s_m-1)*N; 
                B_m(indx_m,indx_m)=A_m{s_m}-gamma*k_m'*k_m/twom_m;
            end

            twomu_m=twomu_m+omega*N*T_m*(T_m-1);
            B_m = B_m + omega*spdiags(ones(N*T_m,2),[-N,N],N*T_m,N*T_m);
            %repeat the Louvain algorithm to obtain robust results

            for i=1:nsamples
                [D_m_it_100(:,i) ,Q_m_it_100(i)] = genlouvain(B_m, 10000,0);

            end
            Qii_m_it_100(:,it) = Q_m_it_100./twomu_m; 
            Dii_m_it_100{it}(:,:) = D_m_it_100;



            % %%%%%%%%%%%%%%%%%%%%%%%% Community - nonmusicians %%%%%%%%%%%%%%%%%%%%%%%%%   
            A_nm{1} = mat_av_norm_nm(:,:,1);
            A_nm{2} = mat_av_norm_nm(:,:,2);
            A_nm{3} = mat_av_norm_nm(:,:,3);
            A_nm{4} = mat_av_norm_nm(:,:,4);
            A_nm{5} = mat_av_norm_nm(:,:,5);

            gamma = 1;

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
            twomu_nm=twomu_nm+omega*N*(T_nm-1);
            B_nm = B_nm + omega*spdiags(ones(N*T_nm,2),[-N,N],N*T_nm,N*T_nm);

            for i=1:nsamples
                [D_nm_it_100(:,i) ,Q_nm_it_100(i)] = genlouvain(B_nm, 10000,0);
            end
            Qii_nm_it_100(:,it) = Q_nm_it_100./twomu_nm; 
            Dii_nm_it_100{it}(:,:) = D_nm_it_100;


            %disp(sprintf('Iteration: %f \n',it))


    end

%% NBS analysis - single layer - non singificant
m_svd1 = squeeze(svd_av_norm_m(:,:,1,:));
avg_m_svd1 = mean(m_svd1,3);
nm_svd1 = squeeze(svd_av_norm_nm(:,:,1,:));
avg_nm_svd1 = mean(nm_svd1,3);

m_svd2 = squeeze(svd_av_norm_m(:,:,2,:));
avg_m_svd2 = mean(m_svd2,3);
nm_svd2 = squeeze(svd_av_norm_nm(:,:,2,:));
avg_nm_svd2 = mean(nm_svd2,3);

m_svd3 = squeeze(svd_av_norm_m(:,:,3,:));
avg_m_svd3 = mean(m_svd3,3);
nm_svd3 = squeeze(svd_av_norm_nm(:,:,3,:));
avg_nm_svd3 = mean(nm_svd3,3);

m_svd4 = squeeze(svd_av_norm_m(:,:,4,:));
avg_m_svd4 = mean(m_svd4,3);
nm_svd4 = squeeze(svd_av_norm_nm(:,:,4,:));
avg_nm_svd4 = mean(nm_svd4,3);

m_svd5 = squeeze(svd_av_norm_m(:,:,5,:));
avg_m_svd5 = mean(m_svd5,3);
nm_svd5 = squeeze(svd_av_norm_nm(:,:,5,:));
avg_nm_svd5 = mean(nm_svd5,3);


%test this range of t values 0.5:0.1:5.0
%degrees of freedom 31+31-2=df
df=31*2-2;
alpha=0.05;
t_values = 0.5:0.1:5.0;

% Calculate p-values
p_values = 2 * (1 - tcdf(abs(t_values), df));

%perform a t test
P=[];
for i=1:5 %over all layers
    
    [~, P{i}] = ttest2(squeeze(svd_av_norm_m(:,:,i,:)), squeeze(svd_av_norm_nm(:,:,i,:)), 'Dim', 3);
    Deci_default=P{i}<0.05;
    %do the bonferroni correction: divide 0.05 by the number of performed test

    %derive the number of performed tests
    N=78;
    t=N*(N-1)/2; %N(N-1)/2

    signif_bon=0.05/t;

    %now let us see which p values survive this correction

    Deci=P{i}<signif_bon;

    imagesc(Deci);

    %yes we have 4 connections left

    [ind,k]=find(Deci==1);

    sig_con{i}=[ind k];
    %none of the connections survive Bonferroni correction
    %now perform NBS
    indx_alpha=0;
    for alpha=1:length(p_values)%0.01:0.01:0.05
        indx_alpha=indx_alpha+1;
        disp(['layer ' num2str(i) ' p limit ' num2str(p_values(alpha))]);
        [p_value_lcc(i,indx_alpha), comp_sizes_T4HC, comps_ori_T4HC, max_sz_T4HC, P_bigcomp_full_T4HC]=nbs_corrected_ttest(P{i},p_values(alpha), squeeze(svd_av_norm_m(:,:,i,:)), squeeze(svd_av_norm_nm(:,:,i,:)), 1000);
        %alpha is the cutoff threshold for the loaded p value matrix
    end
end

%find(p_value_lcc<0.05)
%imagesc(p_value_lcc<0.05)
%how big is the significant component?

%now for the last two tests, let us see what the connected component is

%extract it
%so we need a list of affected connections

%we have the last P_bigcomp_full matrix as a matrix which only includes the
%links of the largest conn component


aff_conn_T4HC=sum(P_bigcomp_full_T4HC{1},2);



%% Matrix plots for Figures 2 and 3
figure 
subplot(321)
imagesc(avg_m_svd1)
colormap jet

subplot(322)
imagesc(avg_m_svd2)
colormap jet

subplot(323)
imagesc(avg_m_svd3)
colormap jet


subplot(324)
imagesc(avg_m_svd4)
colormap jet

subplot(325)
imagesc(avg_m_svd5)
colormap jet
%---------------
figure 
subplot(321)
imagesc(avg_nm_svd1)
colormap jet

subplot(322)
imagesc(avg_nm_svd2)
colormap jet

subplot(323)
imagesc(avg_nm_svd3)
colormap jet


subplot(324)
imagesc(avg_nm_svd4)
colormap jet

subplot(325)
imagesc(avg_nm_svd5)
colormap jet

%%
Q_m_omega{omega_idx}=Qii_m_it_100;
Q_nm_omega{omega_idx}=Qii_nm_it_100;

Qii_temp_it = [mean(Q_m_omega{omega_idx})',mean(Q_nm_omega{omega_idx})'];
figure
notBoxPlot(Qii_temp_it,'jitter',0.5,'style','sdline')
box
xticklabels({'Musicians','Non-musicians'})
ylabel('Modularity (Q)')
title('Modularity between groups')
savefig(['modularity_m' num2str(omega) '.fig'])

%% Build an agreement matrix over all participants

dim=5*78; 

agree_m_100=zeros(dim,dim); %
for AAL_reg1 = 1:dim % Number of AAL regions (need a symmetric matrix at the end)
    for AAL_reg2=1:dim %dim=390, thus for every node in every layer
        for i=1:31 % Number of participants
            for j=1:nsamples %go over the 100 realizations of genlouvain
              Comm_m1=Dii_m_it_100{1,i}(AAL_reg1,j); % Community in musicians
              Comm_m2=Dii_m_it_100{1,i}(AAL_reg2,j);              
              if Comm_m1==Comm_m2
                agree_m_100(AAL_reg1, AAL_reg2)=agree_m_100(AAL_reg1, AAL_reg2)+1;
              end
            end
        end 
    end
end

N=78;
agree_m_mat_100(:,:,1) = agree_m_100(1:N,1:N);
agree_m_mat_100(:,:,2) = agree_m_100(N+1:N*2,N+1:N*2);
agree_m_mat_100(:,:,3) = agree_m_100(N*2+1:N*3,N*2+1:N*3);
agree_m_mat_100(:,:,4) = agree_m_100(N*3+1:N*4,N*3+1:N*4);
agree_m_mat_100(:,:,5) = agree_m_100(N*4+1:N*5,N*4+1:N*5);
    
%now sum them up over all layers

agree_m_col_100=sum(agree_m_mat_100,3);
%put diagonal to zero  
for i=1:N
    agree_m_col_100(i,i)=0;
end

 
%now also make it proportional: 
%maximum is 31*5*1000 = 155 000
zz=31*5*nsamples;
agree_m_col_100=agree_m_col_100./zz;
agree_m_100=agree_m_col_100;

%%%%%%%%%%%%%%%%% non-musicians the same %%%%%%%%%%%%%%%%%%%%%%
dim=5*78;
agree_nm_100=zeros(dim,dim); %
for AAL_reg1 = 1:dim % Number of AAL regions (need a symmetric matrix at the end)
    for AAL_reg2=1:dim %dim=390, thus for every node in every layer
        for i=1:31 % Number of participants
            for j=1:nsamples %go over the 100 realizations of genlouvain
              Comm_m1=Dii_nm_it_100{1,i}(AAL_reg1,j); % Community in non-musicians
              Comm_m2=Dii_nm_it_100{1,i}(AAL_reg2,j); 
              if Comm_m1==Comm_m2
                agree_nm_100(AAL_reg1, AAL_reg2)=agree_nm_100(AAL_reg1, AAL_reg2)+1;
              end     
            end
        end 
    end
end
N=78;
agree_nm_mat_100(:,:,1) = agree_nm_100(1:N,1:N);
agree_nm_mat_100(:,:,2) = agree_nm_100(N+1:N*2,N+1:N*2);
agree_nm_mat_100(:,:,3) = agree_nm_100(N*2+1:N*3,N*2+1:N*3);
agree_nm_mat_100(:,:,4) = agree_nm_100(N*3+1:N*4,N*3+1:N*4);
agree_nm_mat_100(:,:,5) = agree_nm_100(N*4+1:N*5,N*4+1:N*5);
%now sum them up over all layers
agree_nm_col_100=sum(agree_nm_mat_100,3);
%put diagonal to zero  
for i=1:N
    agree_nm_col_100(i,i)=0;
end

%now also make it proportional: 
%maximum is 31*5*1000 = 155 000
zz=31*5*nsamples;
agree_nm_col_100=agree_nm_col_100./zz;
agree_nm_100=agree_nm_col_100;
%% Perform modularity algorithm on collapsed agreement matrices


dim=78;

gamma = 1;
    
u_tot_m=sum(sum(agree_m_100))/2;

B_m=zeros(dim,dim); %%   

degrs=sum(agree_m_100,2);

B_m=agree_m_100-gamma*(degrs*degrs')/(2*u_tot_m);


%nsamples=1000;
for i=1:nsamples
    [D_m_agree_100_it(:,i) ,Q_m_agree_100_it(i)] = genlouvain(B_m, 10000,0);
end
Qii_m_agree_100_it = Q_m_agree_100_it./(2*u_tot_m); 

dim2=78;
agree_m_100_2=zeros(dim2,dim2); %
for AAL_reg1 = 1:dim2 % Number of AAL regions (need a symmetric matrix at the end)
    %agree_m_line = zeros(78,1); % For one layer
    for AAL_reg2=1:dim2 %dim=390, thus for every node in every layer
        
        %for i=1:31 % Number of participants
            for j=1:nsamples %go over the 100 realizations of genlouvain
            %for lyr = 1:3 % Number of layers

              Comm_m1=D_m_agree_100_it(AAL_reg1,j); % Community in musicians
              Comm_m2=D_m_agree_100_it(AAL_reg2,j); 

              %test_m=bsxfun(@eq,Dii_m{1,i},Comm_m); % create test_m for every layer
              if Comm_m1==Comm_m2

                agree_m_100_2(AAL_reg1, AAL_reg2)=agree_m_100_2(AAL_reg1, AAL_reg2)+1;

              end
          
            end

        %end 
    
        %agree_m(AAL_reg,:) = agree_m_line; % Add them in an agreement matrix

    end
end

    
%now also make it proportional: 
%maximum is 1000
agree_m_100_2=agree_m_100_2./nsamples;
%imagesc(agree_m_100_2)

%and now run community detection again on it
dim=78;
gamma = 1;
u_tot_m_2=sum(sum(agree_m_100_2))/2;

B_m=zeros(dim,dim); %%   

degrs=sum(agree_m_100_2,2);

B_m=agree_m_100_2-gamma*(degrs*degrs')/(2*u_tot_m_2);
for i=1:100
    [D_m_agree_100_it_2(:,i) ,Q_m_agree_100_it_2(i)] = genlouvain(B_m, 10000,0);
end
Qii_m_agree_100_it_2 = Q_m_agree_100_it_2./(2*u_tot_m_2); 

%here we saw that we only have one max
[U,~, idx]=unique(D_m_agree_100_it_2', 'rows');
cnt = histc(idx,unique(idx));

if cnt == nsamples
    no_mods_m(omega_idx)=max(D_m_agree_100_it_2(:,1));
else
    disp(omega);
    problem_m(uiui)=omega;
    uiui=uiui+1;
end

%%
%same for non-musicians
dim=78;
gamma = 1;
u_tot_nm=sum(sum(agree_nm_100))/2;

B_nm=zeros(dim,dim); %% 

%we need the degrees per node

degrs_nm=sum(agree_nm_100,2);

B_nm=agree_nm_100-gamma*(degrs_nm*degrs_nm')/(2*u_tot_nm);
%nsamples=1000;
for i=1:nsamples
    [D_nm_agree_100_it(:,i) ,Q_nm_agree_100_it(i)] = genlouvain(B_nm, 10000,0);
end
Qii_nm_agree_100_it = Q_nm_agree_100_it./(2*u_tot_nm); %%%%%%%%%%%%%%%%%%



%---------- for the 100 iterations of the algorithm
dim2=78;
agree_nm_100_2=zeros(dim2,dim2); %
for AAL_reg1 = 1:dim2 % Number of AAL regions (need a symmetric matrix at the end)
    %agree_m_line = zeros(78,1); % For one layer
    for AAL_reg2=1:dim2 %dim=390, thus for every node in every layer
        
        %for i=1:31 % Number of participants
            for j=1:nsamples %go over the 100 realizations of genlouvain
            %for lyr = 1:3 % Number of layers

              Comm_m1=D_nm_agree_100_it(AAL_reg1,j); % Community in musicians
              Comm_m2=D_nm_agree_100_it(AAL_reg2,j); 

              %test_m=bsxfun(@eq,Dii_m{1,i},Comm_m); % create test_m for every layer
              if Comm_m1==Comm_m2

                agree_nm_100_2(AAL_reg1, AAL_reg2)=agree_nm_100_2(AAL_reg1, AAL_reg2)+1;

              end
          
            end

        %end 
    
        %agree_m(AAL_reg,:) = agree_m_line; % Add them in an agreement matrix

    end
end

    
%now also make it proportional: 
%maximum is 100
agree_nm_100_2=agree_nm_100_2./nsamples;
%imagesc(agree_nm_100_2)

%and now run community detection again on it
dim=78;
gamma = 1;
u_tot_nm_2=sum(sum(agree_nm_100_2))/2;

B_nm=zeros(dim,dim); %%   

degrs=sum(agree_nm_100_2,2);

B_nm=agree_nm_100_2-gamma*(degrs*degrs')/(2*u_tot_nm_2);
for i=1:nsamples
    [D_nm_agree_100_it_2(:,i) ,Q_nm_agree_100_it_2(i)] = genlouvain(B_nm, 10000,0);
end
Qii_nm_agree_100_it_2 = Q_nm_agree_100_it_2./(2*u_tot_nm_2); 


%imagesc(D_nm_agree_100_it_2)
%we need to make sure these are just stripes
[U,~, idx]=unique(D_nm_agree_100_it_2', 'rows');
cnt = histc(idx,unique(idx));

if cnt == nsamples
    no_mods_nm(omega_idx)=max(D_nm_agree_100_it_2(:,1));
else
    disp(omega);
    problem_nm(uiui)=omega;
    uiui=uiui+1;
    %we have a case in here and then we need for omega=1m repeat the
    %process one more time
    dim2=78;
    agree_nm_100_3=zeros(dim2,dim2); %
    for AAL_reg1 = 1:dim2 % Number of AAL regions (need a symmetric matrix at the end)
        %agree_m_line = zeros(78,1); % For one layer
        for AAL_reg2=1:dim2 %dim=390, thus for every node in every layer

            %for i=1:31 % Number of participants
                for j=1:nsamples %go over the 100 realizations of genlouvain
                %for lyr = 1:3 % Number of layers

                  Comm_m1=D_nm_agree_100_it_2(AAL_reg1,j); % Community in musicians
                  Comm_m2=D_nm_agree_100_it_2(AAL_reg2,j); 

                  %test_m=bsxfun(@eq,Dii_m{1,i},Comm_m); % create test_m for every layer
                  if Comm_m1==Comm_m2

                    agree_nm_100_3(AAL_reg1, AAL_reg2)=agree_nm_100_3(AAL_reg1, AAL_reg2)+1;

                  end

                end

            %end 

            %agree_m(AAL_reg,:) = agree_m_line; % Add them in an agreement matrix

        end
    end


    %now also make it proportional: 
    %maximum is 100
    agree_nm_100_2=agree_nm_100_3./nsamples;
    %imagesc(agree_nm_100_2)

    %and now run community detection again on it
    dim=78;
    gamma = 1;
    u_tot_nm_2=sum(sum(agree_nm_100_3))/2;

    B_nm=zeros(dim,dim); %%   

    degrs=sum(agree_nm_100_3,2);

    B_nm=agree_nm_100_3-gamma*(degrs*degrs')/(2*u_tot_nm_2);
    for i=1:nsamples
        [D_nm_agree_100_it_3(:,i) ,Q_nm_agree_100_it_3(i)] = genlouvain(B_nm, 10000,0);
    end
    Qii_nm_agree_100_it_3 = Q_nm_agree_100_it_3./(2*u_tot_nm_2); %%%%%%%%%%%%%%%%%%


    %imagesc(D_nm_agree_100_it_3)
    %we need to make sure these are just stripes
    [U,~, idx]=unique(D_nm_agree_100_it_3', 'rows');
    cnt = histc(idx,unique(idx));

    if cnt == nsamples
        no_mods_nm(omega_idx)=max(D_nm_agree_100_it_3(:,1));
        D_nm_agree_100_it_2=D_nm_agree_100_it_3;
    else
        disp(omega);
        problem_nm_again(uiui)=omega;
        uiui=uiui+1;
        dim2=78;
        agree_nm_100_4=zeros(dim2,dim2); %
        for AAL_reg1 = 1:dim2 % Number of AAL regions (need a symmetric matrix at the end)
            %agree_m_line = zeros(78,1); % For one layer
            for AAL_reg2=1:dim2 %dim=390, thus for every node in every layer

                %for i=1:31 % Number of participants
                    for j=1:nsamples %go over the 100 realizations of genlouvain
                    %for lyr = 1:3 % Number of layers

                      Comm_m1=D_nm_agree_100_it_3(AAL_reg1,j); % Community in musicians
                      Comm_m2=D_nm_agree_100_it_3(AAL_reg2,j); 

                      %test_m=bsxfun(@eq,Dii_m{1,i},Comm_m); % create test_m for every layer
                      if Comm_m1==Comm_m2

                        agree_nm_100_4(AAL_reg1, AAL_reg2)=agree_nm_100_4(AAL_reg1, AAL_reg2)+1;

                      end

                    end

                %end 

                %agree_m(AAL_reg,:) = agree_m_line; % Add them in an agreement matrix

            end
        end


        %now also make it proportional: 
        %maximum is 100
        agree_nm_100_2=agree_nm_100_4./nsamples;
        %imagesc(agree_nm_100_2)

        %and now run community detection again on it
        dim=78;
        gamma = 1;
        u_tot_nm_2=sum(sum(agree_nm_100_4))/2;

        B_nm=zeros(dim,dim); %%   

        degrs=sum(agree_nm_100_4,2);

        B_nm=agree_nm_100_4-gamma*(degrs*degrs')/(2*u_tot_nm_2);
        for i=1:nsamples
            [D_nm_agree_100_it_4(:,i) ,Q_nm_agree_100_it_4(i)] = genlouvain(B_nm, 10000,0);
        end
        Qii_nm_agree_100_it_4 = Q_nm_agree_100_it_4./(2*u_tot_nm_2); %%%%%%%%%%%%%%%%%%


        %imagesc(D_nm_agree_100_it_4)
        %we need to make sure these are just stripes
        [U,~, idx]=unique(D_nm_agree_100_it_4', 'rows');
        cnt = histc(idx,unique(idx));

        if cnt == nsamples
            no_mods_nm(omega_idx)=max(D_nm_agree_100_it_4(:,1));
            D_nm_agree_100_it_2=D_nm_agree_100_it_4;
        else
            disp(omega);
            problem_nm_again_again(uiui)=omega;
            uiui=uiui+1;

        end
        
    end
    
end

%% Calculate the rand index on the group-level module assignments

addpath('/Users/jilmeier/Documents/musician')

%rand index over all partitions for musicians and non-musicians
[zRand, SR, SAR(omega_idx), VI] = zrand(D_m_agree_100_it_2(:,1),D_nm_agree_100_it_2(:,1));

%save the community assignments
D_nm_omega(:,omega_idx)=D_nm_agree_100_it_2(:,1);
D_m_omega(:,omega_idx)=D_m_agree_100_it_2(:,1);

%savename = ['5layers_results_17052024_omega'  num2str(omega)  '.mat']
%save([savename])
disp(omega);
end

%savename = ['5layers_results_17052024_omega_all.mat']
%save([savename])

%plots with omega variable
% load this workspace in case you do not want to rerun the above loop
%load('/Users/jilmeier/Documents/musician/code/OMEGA_Feb24/upload to git/upload to git for revision/do not upload - workspaces needed to reproduce results/workspace_for_SVD_results.mat')
omega_var=0:0.01:1;

%% Figure 4 plot over all interlayer coupling values
%add standard deviations
for om=1:length(omega_var)
    Q_m_std(om)=std(Q_m_omega{om}(:));
    Q_nm_std(om)=std(Q_nm_omega{om}(:));
end

figure;
errorbar(omega_var, Q_m_mean, Q_m_std, "LineWidth",2)
hold on;
errorbar(omega_var, Q_nm_mean, Q_nm_std, "LineWidth",2)
xlabel('interlayer coupling c')
ylabel('Modularity (Q)')
legend('Musicians', 'Non-Musicians', 'FontSize', 20, 'Location', 'northwest')
ax=gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
savefig(['modularity_omega.fig'])

%% Average Rand index over all interlayer coupling values

mean(SAR)

%ans =

%    0.6876

%% Number of modules
plot(omega_var,no_mods_nm)
hold on;
plot(omega_var,no_mods_m)


%% Statistics for modularity over all interlayer coupling values c
for om=1:length(Q_m_omega)
    [P_Qii_it(omega_idx) H_Qii_it STATS] = ranksum(mean(Q_m_omega{om}),mean(Q_nm_omega{om}));
    zstats(om)=STATS.zval;
end
%plot(omega_var,P_Qii_it)
%always significantly different
%imagesc(P_Qii_it<0.05)
%correct for multiple comparisons
[h_dc, crit_p_dc, adj_ci_cvrg_dc, corrected_p_values_P_Qii_it] = fdr_bh(P_Qii_it(:),0.05,'dep','yes');


%% Relation between modularity and age
%but now take all omega values together and recalculate Spearman
ages_m=[
30
50
42
28
27
19
31
30
28
28
26
26
23
20
27
26
24
20
20
41
27
21
25
27
20
33
20
34
33
43
27    
];

ages_nm=[
 26
31
32
31
32
33
22
25
28
39
24
21
28
36
21
26
24
32
30
33
35
27
41
28
19
22
21
21
23
23
18];


for i=1:length(Q_m_omega)
   Q_m_all(:,:,i)= Q_m_omega{i};
   Q_nm_all(:,:,i)= Q_nm_omega{i};
end
Q_combi_all=[Q_m_all Q_nm_all];
YEARS=repmat([ages_m' ages_nm'],nsamples,1, 101);
[R,P]=corr(YEARS(:), Q_combi_all(:), 'type', 'Spearman');

%now correct for age variable in the years of training vs. modularity
%relationship

%% relation between modularity and years of training
%load('workspace_for_SVD_results.mat')
%omega_var=0:0.01:1;
% Years of musical experience
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

YEARS=repmat(yrs_exp',nsamples,1);
big_vector_Q_m=[];
big_vector_years_m=[];
for i = 1:numel(Q_m_omega)
    % Convert each matrix into a vector and concatenate
    big_vector_Q_m = vertcat(big_vector_Q_m, Q_m_omega{i}(:));
    big_vector_years_m=vertcat(big_vector_years_m, YEARS(:));
end

[R_Exp,P_Exp]=corr(big_vector_years_m, big_vector_Q_m, 'type', 'Spearman');

YEARS=repmat(ages_m',nsamples,1);
big_vector_Q_m=[];
big_vector_years_m=[];
for i = 1:numel(Q_m_omega)
    % Convert each matrix into a vector and concatenate
    big_vector_Q_m = vertcat(big_vector_Q_m, Q_m_omega{i}(:));
    big_vector_years_m=vertcat(big_vector_years_m, YEARS(:));
end

[R_AGE,P_AGE]=corr(big_vector_years_m, big_vector_Q_m, 'type', 'Spearman');

rho_yage=R_AGE;
rho_partial = (R_Exp - rho_yage * rho_yage) / sqrt((1 - rho_yage * rho_yage) * (1 - rho_yage * rho_yage));
%this is the corrected spearman correlation rho_partial
% Calculate the number of observations (n)
n = length(big_vector_years_m);

% Calculate the standard error of the partial correlation coefficient
se_rho_partial = sqrt((1 - rho_partial^2) / (n - 3));

% Calculate the t-statistic
t_stat = rho_partial / se_rho_partial;

% Calculate the degrees of freedom
df = n - 3;

% Calculate the two-tailed p-value
p_value = 2 * (1 - tcdf(abs(t_stat), df));


%% Relation between modularity and age at which participants started taking lessons

%age - years of training
for i=1:length(Q_m_omega)
   Q_m_all(:,:,i)= Q_m_omega{i};
   Q_nm_all(:,:,i)= Q_nm_omega{i};
end
YEARS=repmat([ages_m-yrs_exp]',nsamples,1, 101);
[R,P]=corr(YEARS(:), Q_m_all(:), 'type', 'Spearman');


%% predictive model: to predict whether somebody is a musicians based on
% multilayer modularity

%make an average over all omega and plot it
average_over_cells = [];  % Initialize with zeros or NaNs
for i = 1:numel(Q_m_omega)
    average_over_cells(:,i) = mean(Q_m_omega{i});
end
average_Q_m=mean(average_over_cells,2);

average_over_cells = [];  % Initialize with zeros or NaNs
for i = 1:numel(Q_m_omega)
    average_over_cells(:,i) = mean(Q_nm_omega{i});
end
average_Q_nm=mean(average_over_cells,2);


avg_modularity = [average_Q_m; average_Q_nm]'; % average modularity values
is_musician = [ones(31,1); zeros(31,1)]'; % 1 indicates musician, 0 indicates non-musician

% Combine data into a table
data = table(avg_modularity', is_musician', 'VariableNames', {'AvgModularity', 'IsMusician'});

% Fit logistic regression model
model = fitglm(data, 'IsMusician ~ AvgModularity', 'Distribution', 'binomial', 'Link', 'logit');

% Display model summary
disp(model);

% Predict probabilities
predicted_probs = predict(model, data);

% Convert probabilities to binary predictions using 0.5 threshold
predicted_classes = predicted_probs >= 0.5;
predicted_classes = double(predicted_classes);
% Confusion matrix
confusion_mat = confusionmat(data.IsMusician, predicted_classes);
disp('Confusion Matrix:');
disp(confusion_mat);

% Calculate accuracy
accuracy = sum(predicted_classes == data.IsMusician) / length(data.IsMusician);
disp(['Accuracy: ', num2str(accuracy)]);

% Plot ROC curve
[X, Y, T, AUC] = perfcurve(data.IsMusician, predicted_probs, 1);
figure;
plot(X, Y);
xlabel('False positive rate');
ylabel('True positive rate');
title(['ROC Curve (AUC = ', num2str(AUC), ')']);

% Generate a range of values for AvgModularity
x_range = linspace(min(avg_modularity), max(avg_modularity), 100)';

% Create a table for prediction
x_range_table = table(x_range, 'VariableNames', {'AvgModularity'});

% Predict probabilities using the logistic regression model
y_range_probs = predict(model, x_range_table);

% Plot the data points
figure;
hold on;
gscatter(avg_modularity, is_musician, is_musician, 'br', 'xo');
plot(x_range, y_range_probs, 'k-', 'LineWidth', 2);
xlabel('Average Modularity');
ylabel('Probability of Being a Musician');
title('Logistic Regression - Probability of Being a Musician vs. Average Modularity');
legend('Non-Musician', 'Musician', 'Logistic Regression', 'Location', 'Best');
hold off;


%% predictive model over different omega

for jj=1:length(omega_var)
    avg_modularity=[mean(Q_m_omega{jj}) mean(Q_nm_omega{jj})];
    
    is_musician = [ones(31,1); zeros(31,1)]'; % 1 indicates musician, 0 indicates non-musician

    % Combine data into a table
    data = table(avg_modularity', is_musician', 'VariableNames', {'AvgModularity', 'IsMusician'});

    % Fit logistic regression model
    model = fitglm(data, 'IsMusician ~ AvgModularity', 'Distribution', 'binomial', 'Link', 'logit');

    % Display model summary
    disp(model);

    % Predict probabilities
    predicted_probs = predict(model, data);

    % Convert probabilities to binary predictions using 0.5 threshold
    predicted_classes = predicted_probs >= 0.5;
    predicted_classes = double(predicted_classes);
    % Confusion matrix
    confusion_mat = confusionmat(data.IsMusician, predicted_classes);
    disp('Confusion Matrix:');
    disp(confusion_mat);

    % Calculate accuracy
    accuracy(jj) = sum(predicted_classes == data.IsMusician) / length(data.IsMusician);
    disp(['Accuracy: ', num2str(accuracy(jj))]);

    

end

plot(omega_var,accuracy)
xlabel('omega','FontSize', 20)
ylabel('Accuracy of predictive model', 'FontSize', 20)

%% 4-fold cross-validation

avg_modularity = [average_Q_m; average_Q_nm]'; % average modularity values
is_musician = [ones(31,1); zeros(31,1)]'; % 1 indicates musician, 0 indicates non-musician

% Combine data into a table
data = table(avg_modularity', is_musician', 'VariableNames', {'AvgModularity', 'IsMusician'});
n=62;
% Randomly shuffle the data
rng('default'); % For reproducibility
shuffled_indices = randperm(n);
data = data(shuffled_indices, :);

% Define the number of folds
k = 4;

% Split data into k folds
cv = cvpartition(n, 'KFold', k);

% Initialize variables to store results
accuracy = zeros(k, 1);
predicted_probs_all = cell(k, 1);
predicted_classes_all = cell(k, 1);
actual_classes_all = cell(k, 1);

for i = 1:k
    % Get training and test indices for the current fold
    train_idx = cv.training(i);
    test_idx = cv.test(i);
    
    % Prepare training and test data
    train_data = data(train_idx, :);
    test_data = data(test_idx, :);
    
    % Train the logistic regression model
    model = fitglm(train_data, 'IsMusician ~ AvgModularity', 'Distribution', 'binomial', 'Link', 'logit');
    
    % Predict probabilities for the test set
    predicted_probs = predict(model, test_data);
    
    % Convert probabilities to binary predictions using 0.5 threshold
    predicted_classes = predicted_probs >= 0.5;
    
    % Store results
    predicted_probs_all{i} = predicted_probs;
    predicted_classes_all{i} = predicted_classes;
    actual_classes_all{i} = test_data.IsMusician;
    
    % Calculate accuracy for the current fold
    accuracy(i) = sum(predicted_classes == test_data.IsMusician) / length(test_data.IsMusician);
end

% Average accuracy across all folds
mean_accuracy = mean(accuracy);
disp(['Mean Accuracy: ', num2str(mean_accuracy)]);

% Aggregate all predicted probabilities and actual classes
all_predicted_probs = cell2mat(predicted_probs_all);
all_actual_classes = cell2mat(actual_classes_all);

% Plot ROC curve
[X, Y, T, AUC] = perfcurve(all_actual_classes, all_predicted_probs, 1);
figure;
plot(X, Y);
xlabel('False positive rate');
ylabel('True positive rate');
title(['ROC Curve (AUC = ', num2str(AUC), ')']);


%% 4-fold cross-validation for each omega

for jj=1:length(omega_var)
    avg_modularity=[mean(Q_m_omega{jj}) mean(Q_nm_omega{jj})];
    
    is_musician = [ones(31,1); zeros(31,1)]'; % 1 indicates musician, 0 indicates non-musician

    % Combine data into a table
    data = table(avg_modularity', is_musician', 'VariableNames', {'AvgModularity', 'IsMusician'});

    % Randomly shuffle the data
    rng('default'); % For reproducibility
    shuffled_indices = randperm(n);
    data = data(shuffled_indices, :);

    % Define the number of folds
    k = 4;

    % Split data into k folds
    cv = cvpartition(n, 'KFold', k);

    % Initialize variables to store results
    accuracy = zeros(k, 1);
    predicted_probs_all = cell(k, 1);
    predicted_classes_all = cell(k, 1);
    actual_classes_all = cell(k, 1);

    for i = 1:k
        % Get training and test indices for the current fold
        train_idx = cv.training(i);
        test_idx = cv.test(i);

        % Prepare training and test data
        train_data = data(train_idx, :);
        test_data = data(test_idx, :);

        % Train the logistic regression model
        model = fitglm(train_data, 'IsMusician ~ AvgModularity', 'Distribution', 'binomial', 'Link', 'logit');

        % Predict probabilities for the test set
        predicted_probs = predict(model, test_data);

        % Convert probabilities to binary predictions using 0.5 threshold
        predicted_classes = predicted_probs >= 0.5;

        % Store results
        predicted_probs_all{i} = predicted_probs;
        predicted_classes_all{i} = predicted_classes;
        actual_classes_all{i} = test_data.IsMusician;

        % Calculate accuracy for the current fold
        accuracy(i) = sum(predicted_classes == test_data.IsMusician) / length(test_data.IsMusician);
    end

    % Average accuracy across all folds
    mean_accuracy(jj) = mean(accuracy);
    disp(['Mean Accuracy: ', num2str(mean_accuracy(jj))]);

    % Aggregate all predicted probabilities and actual classes
    all_predicted_probs = cell2mat(predicted_probs_all);
    all_actual_classes = cell2mat(actual_classes_all);

    % Plot ROC curve
    [X, Y, T, AUC] = perfcurve(all_actual_classes, all_predicted_probs, 1);
%     figure;
%     plot(X, Y);
%     xlabel('False positive rate');
%     ylabel('True positive rate');
%     title(['ROC Curve (AUC = ', num2str(AUC), ')']);

end

plot(omega_var,mean_accuracy)
xlabel('interlayer coupling c','FontSize', 20)
ylabel('Accuracy of predictive model', 'FontSize', 20)





%% prediction of training years

% Example data
avg_modularity = [average_Q_m]'; % average modularity values
is_musician = [ones(31,1)]'; % 1 indicates musician, 0 indicates non-musician
training_years = yrs_exp'; % training years (only for musicians)

% Combine data into a table
data = table(avg_modularity', is_musician', training_years', 'VariableNames', {'AvgModularity', 'IsMusician', 'TrainingYears'});

% Filter data to only include musicians
musician_data = data(data.IsMusician == 1, :);
musician_data = rmmissing(musician_data); % Remove rows with missing values

% Fit linear regression model
model = fitlm(musician_data, 'TrainingYears ~ AvgModularity');

% Display model summary
disp(model);

% Predicted training years
predicted_training_years = predict(model, musician_data);

% Plot the results
figure;
scatter(musician_data.AvgModularity, musician_data.TrainingYears, 'b', 'filled');
hold on;
plot(musician_data.AvgModularity, predicted_training_years, 'r', 'LineWidth', 2);
xlabel('Average Modularity');
ylabel('Training Years');
title('Linear Regression - Training Years vs. Average Modularity');
legend('Observed Data', 'Fitted Line', 'Location', 'Best');
hold off;

% Calculate performance metrics
residuals = musician_data.TrainingYears - predicted_training_years;
rmse = sqrt(mean(residuals.^2));
disp(['Root Mean Squared Error: ', num2str(rmse)]);


%% Figure 5: Plot on template brain if module assignment is higher than chance
% color of the corresponding module if frequency of participation in that
% module over all omega is larger than 0.33
%check for the main module with regard to frequency
%load('workspace_for_SVD_results.mat')

numCols = size(D_nm_omega, 2);

% Threshold for more than 33.333% frequency
threshold = numCols * (1/3);

% Initialize an array to store the result for each row
result = NaN(size(D_nm_omega, 1), 1);

% Loop over each row
for i = 1:size(D_nm_omega, 1)
    % Get the current row
    row = D_nm_omega(i, :);
    
    % Count occurrences of each unique element in the row
    uniqueVals = unique(row);
    counts = histc(row, uniqueVals); % Use histc to count occurrences of each unique value
    
    % Find values that occur more frequently than the threshold
    moreFrequentValues = uniqueVals(counts > threshold);
    
    % Store the result if there is a value meeting the criteria
    if ~isempty(moreFrequentValues)
        result(i) = moreFrequentValues(1); % Take the first frequent value found
    end
end

%imagesc(mean(D_m_omega==1,2))
addpath('/Users/jilmeier/Documents/Virtual_DBS/code/BrainVisual')
data_m=result;

% 
colourbar_threshold=[]; % can be used to adjust the colour range (experimental)
mesh_type = 'spm_canonical'; % assume that input contains 78 AAL ROIs
nr_views=6; % #views of the cortical surface in the figures
colour_range=[min(data_m),max(data_m)]; % for display: colour_range will be based on the data; alternatively, you can provide a maximum and minimum value
%}

% get AAL labels
 
[aalID, aalind,fullnames,everyID,allnames] = aal_get_numbers( 'Precentral_L' );
        tmplabels = char(allnames);
        cfg.allnames=tmplabels;
        
% Use only the most superfial areas
indices_in_same_order_as_in_Brainwave = select_ROIs_from_full_AAL(cfg);
labels = tmplabels(indices_in_same_order_as_in_Brainwave,:); %78 labels

% plot

[colourbar_handle, patch_handles] = PaintBrodmannAreas_new2_clean_musician(labels, data_m, length(data_m),length(data_m),nr_views, colour_range, colourbar_threshold, mesh_type);
set(gcf,'Tag','ShowBrainFigure');
%tit
savefig(['modules_omega_nonmusicians_frequency_above1third.fig'])

% do the same for musicians

% color of the corresponding module if frequency of participation in that
% module over all omega is larger than 0.33

%D_nm_omega_mean=mean(D_nm_omega==1,2); %average times over all omega how often this region is 
%in the community 1
%check for the main module with regard to frequency
%D_nm_omega_mean=D_nm_omega_mean>0.3333;

numCols = size(D_m_omega, 2);

% Threshold for more than 33.333% frequency
threshold = numCols * (1/2);

% Initialize an array to store the result for each row
result = NaN(size(D_m_omega, 1), 1);

% Loop over each row
for i = 1:size(D_m_omega, 1)
    % Get the current row
    row = D_m_omega(i, :);
    
    % Count occurrences of each unique element in the row
    uniqueVals = unique(row);
    counts = histc(row, uniqueVals); % Use histc to count occurrences of each unique value
    
    % Find values that occur more frequently than the threshold
    moreFrequentValues = uniqueVals(counts > threshold);
    
    % Store the result if there is a value meeting the criteria
    if ~isempty(moreFrequentValues)
        result(i) = moreFrequentValues(1); % Take the first frequent value found
    end
end

% Display result
disp('Values more frequent than chance (50%) in each row:');
disp(result);

%imagesc(mean(D_m_omega==1,2))
addpath('/Users/jilmeier/Documents/Virtual_DBS/code/BrainVisual')
data_m=result;


% 
colourbar_threshold=[]; % can be used to adjust the colour range (experimental)
mesh_type = 'spm_canonical'; % assume that input contains 78 AAL ROIs
nr_views=6; % #views of the cortical surface in the figures
colour_range=[min(data_m),max(data_m)]; % for display: colour_range will be based on the data; alternatively, you can provide a maximum and minimum value
%}

% get AAL labels
 
[aalID, aalind,fullnames,everyID,allnames] = aal_get_numbers( 'Precentral_L' );
        tmplabels = char(allnames);
        cfg.allnames=tmplabels;
        
% Use only the most superfial areas
indices_in_same_order_as_in_Brainwave = select_ROIs_from_full_AAL(cfg);
labels = tmplabels(indices_in_same_order_as_in_Brainwave,:); %78 labels

% plot

[colourbar_handle, patch_handles] = PaintBrodmannAreas_new2_clean_musician(labels, data_m, length(data_m),length(data_m),nr_views, colour_range, colourbar_threshold, mesh_type);
set(gcf,'Tag','ShowBrainFigure');
%tit
savefig(['modules_omega_musicians_frequency_above50percent.fig'])



