%% comparison between SVD and no SVD normalization results

load('workspace_for_noSVD.mat')
%save it as noSVD:
D_nm_noSVD=D_nm_agree_omega(:,1);

load('workspace_for_SVD_results.mat')

D_nm_SVD=D_nm_agree_omega(:,1);


[zRand, SR, SAR_SVD, VI] = zrand(D_nm_SVD,D_nm_noSVD);



imagesc(D_nm_SVD-D_nm_noSVD)

%plot the difference on brain

addpath('/Users/jilmeier/Documents/Virtual_DBS/code/BrainVisual')
data_m=D_nm_SVD-D_nm_noSVD;


% 
colourbar_threshold=[]; % can be used to adjust the colour range (experimental)
mesh_type = 'spm_canonical'; % assume that input contains 78 AAL ROIs
nr_views=6; % #views of the cortical surface in the figures
colour_range=[min(data_m),max(data_m)]; % for display: colour_range will be based on the data; alternatively, you can provide a maximum and minimum value
%}

%% get AAL labels
 
[aalID, aalind,fullnames,everyID,allnames] = aal_get_numbers( 'Precentral_L' );
        tmplabels = char(allnames);
        cfg.allnames=tmplabels;
        
% Use only the most superfial areas
indices_in_same_order_as_in_Brainwave = select_ROIs_from_full_AAL(cfg);
labels = tmplabels(indices_in_same_order_as_in_Brainwave,:); %78 labels

%% plot

addpath('/Users/jilmeier/Documents/Virtual_DBS/code/BrainVisual')

[colourbar_handle, patch_handles] = PaintBrodmannAreas_new2_clean_musician(labels, data_m, length(data_m),length(data_m),nr_views, colour_range, colourbar_threshold, mesh_type);
set(gcf,'Tag','ShowBrainFigure');
%tit
