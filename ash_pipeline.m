%% create lesion and tdcs simulation

clear; clc;
%% cd in the folder where the input data are and set path to folders
cd path2input_data

[path2mri, ~, ~] = fileparts(pwd);
path2mri = [path2mri, filesep, 'data', filesep];
sysvar = ';';

figfolder = [path2mri,'figures'];
if ~exist(figfolder, 'dir')
    system(sprintf('mkdir %s', figfolder));
end

LesionConductivity = 0.126:0.1:1.654;
idx_subj = [34,35,41,42,44,46:48,50,51,53:57];
%% addpath to toolboxes

addpath('path2ash_project')
addpath('path2SimNIBS/matlab') 
addpath('path2fieldtrip/')
ft_defaults; 
addpath path2fieldtrip/external/iso2mesh/


%% create General Head Model with SimNIBS

idx_linda = [34,35,41,42,44,46:48,50,51,53:57];
for i=1:length(idx_linda)
    subjID = ['sub-0', num2str(idx_linda(i))];
    fprintf(sprintf('SIMNIBS - subject: %s\n', subjID));
        
    foldername = [path2mri, subjID];
    if ~exist(foldername, 'dir')
        system(sprintf('mkdir %s', foldername));
        system(sprintf('mv %s %s', ...
            [path2mri, subjID, '_mri_linda.nii'], ...
            [path2mri, subjID, filesep, subjID, '_mri_linda.nii']));
    end
    c = sprintf('cd %s; headreco all --cat %s %s_mri_linda.nii -d no-conform', foldername, subjID, subjID);
    if ~exist([foldername, filesep, subjID, '.msh'], 'file')
        system(c);
    end
end

%% use LINDA to segment the lesion
for subj=1:length(idx_subj)
    subjID = ['sub-0', num2str(idx_subj(subj))];
    path2msh_folder = [path2mri, subjID, filesep];
    path2msh_file = [path2msh_folder, subjID,'_mri_linda.nii'];
    fileID = fopen([path2msh_folder,'run_linda_',subjID,'.R'],'w');
    fprintf(fileID, '#!/usr/bin/env Rscript\n');
    fprintf(fileID, 'library(LINDA)\n');
    fprintf(fileID, ['output = linda_predict(''',path2msh_file,''')\n']);
    fclose(fileID);

    system(['chmod +x ',path2msh_folder,'run_linda_',subjID,'.R']);
    
    fh_file = ['run_linda_',subjID,'.R'];
    cmd = sprintf('cd %s ; Rscript %s', path2msh_folder, fh_file);
    system(cmd)

end

%% create lesion mesh
for subj=1:length(idx_subj)
    subjID = ['sub-0', num2str(idx_subj(subj))];
    path2linda = [path2mri, subjID,'/linda/Prediction3_native.nii.gz'];
    path2msh_folder = [path2mri, subjID, filesep];
    path2msh_file = [path2msh_folder, subjID,'.msh'];
    tic
    [lesion_hm, idx_lesion_tet] = create_lesion_mesh(path2msh_file,path2linda,subjID,path2mri);
    toc
    mesh_save_gmsh4(lesion_hm, [path2msh_folder,subjID,'_lesion.msh' ]);
    toc
end
%% visualize
% for subj=1:length(idx_subj)
%     subjID = ['sub-0', num2str(idx_subj(subj))];
%     path2linda = [path2mri, subjID,'/linda/Prediction3_native.nii.gz'];
%     path2msh_folder = [path2mri, subjID,'/'];
%     path2msh_file = [path2msh_folder, subjID,'.msh'];
%     %%
%     with_lesion = mesh_load_gmsh4([path2msh_folder,subjID,'_lesion.msh' ]);
%     mesh_show_surface(with_lesion)
%     hold on,
%     mesh_show_surface(with_lesion,'region_idx',1011)
% 
% end

%% run tdcs simulations
for subj=3:length(idx_subj)
    subjID = ['sub-0', num2str(idx_subj(subj))];
    path2msh_folder = [path2mri, subjID,'/'];
    path2msh_file = [path2msh_folder, subjID,'.msh'];
    run_tdcs_simulation(path2msh_folder, path2msh_file,subj, idx_subj);   
end

%% produce figures

run_analysis
 
%% measure lesion and plot figure volumes/max rel diff

for subj=1:length(idx_subj)
    subjID = ['sub-0', num2str(idx_subj(subj))];
    path2linda = [path2mri, subjID,'/linda/Prediction3_native.nii.gz'];
    path2msh_folder = [path2mri, subjID,'/'];
    path2msh_file = [path2msh_folder, subjID,'.msh'];
    geo = mesh_load_gmsh4([path2msh_folder,'sub-0', num2str(idx_subj(subj)),'_lesion.msh']);
    info_lesion.data(subj,1)= sum(geo.tetrahedron_regions==11);

    info_lesion.data(subj,2) = sum(elemvolume(geo.nodes,geo.tetrahedra(geo.tetrahedron_regions==11,:)));

    mask = ft_read_mri(path2linda);
    info_lesion.data(subj,3) = sum(sum(sum(mask.anatomy)));
end

save info_lesion info_lesion 
%%
load([path2mri, 'info_lesion.mat'])
load([path2mri, 'max_diff.mat'])

ii = info_lesion.data(:,2);
ii(13)=[];
max_diff(13) = [];
[a,b] = sort(ii);
a(end+1) = 1.83e+05;    %lesion volume subj401
c = max_diff(b);
c(end+1) = 63;          %max err for subj401
%%

figure, plot(a/1000,c,'.', 'MarkerSize',40,'linewidth', 3,'color', [0 0.4470 0.7410]),
set(gca,'Fontsize',28)
xlabel('Lesion volume (cm^3)');
ylabel('Max relative difference in ipsi target (%)');
set(gcf, 'Position', get(0, 'Screensize'));

saveas(gcf,[path2mri,'/figures/vol_err_noline.png'])

%% stats on lesion volume

iqr_lesion = iqr(a); %8.9856e+04
min_lesion = min(a); %2.5962e+03 
max_lesion = max(a); %183000
median_lesion = median(a); %3.7655e+04




