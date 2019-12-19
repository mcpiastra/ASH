%% TOTAL script
% Creates a Volume Conduction Model for chronic stroke patients
% Compares the developed volume conduction model with outcomes of SimNIBS 

% Requirements pipeline:
%   SimNIBS 2.1.2 (including: SPM12 and CAT12)
%   LINDA version 0.5.0
%   Fieldtrip 
%   R version 3.5.3 (2019-03-11)
%   Add SimNIBS and Fieltrip to the Matlab path
%   Replace the \simnibs_2.1.2\matlab\functions\standard_cond.m matlab function of SimNIBS 

% USAGE:
%   LesionConductivity: Lesion conductivity values 
%   cond_lesion:  Lesion conductivity value without decimals 
%   To run segmentations uncomment line 31 (commented due to time
%   consuming process)
%   To run simulations uncomment line 135 and line 151 (commented due to time
%   consuming process)

clear all;
close all;
clc;
%% Segmententation whole head
% SimNIBS
c = sprintf('cd Data/401 ; headreco all --cat 401 401_T1_L.nii -d no-conform');   
                                        % all = all reconstructions steps, 
                                        % 401 = name output folder, 
                                    	% 401_T1_L.nii = MRI data 
                                      	% -d no-conform = to keep the axis of the MRI (Saturnino et al., 2018)
%system(c);                              % Call SimNIBS from the terminal

% Load Segmented Data
% Segmentation mesh by SimNIBS
m1_seg = mesh_load_gmsh4([pwd, filesep, 'Data/401/401_5.msh']); % load mesh segmented by SimNIBS into scalp, skull, GM, WM and CSF 
m2_seg = mesh_load_gmsh4([pwd, filesep, 'Data/401/401_5.msh']); % load mesh segmented by SimNIBS into scalp, skull, GM, WM, CSF and lesion 
%% Segmentation lesion
% LINDA
b = fprintf('Rscript library(LINDA); output=linda_predict(cd Data/401/401_T1_L.nii)');  
system(b);

% Load Lesion mask segmented by LINDA:
lesion = ft_read_mri('Data/401/linda/Prediction3_native.nii.gz'); % load lesion mask
%% Creating a mesh 
tic;
% 1. Find coordinates of lesion mask
[r,c,v] = ind2sub(size(lesion.anatomy),find(lesion.anatomy == 1)); % Coordinates lesion voxels
lesion_vox = [r c v];                                               % voxels of the lesion
sc = lesion_vox*lesion.transform(1:3,1:3);
lesion_coordinates = sc + repmat(lesion.transform(1:3,4)', size(sc,1),1); % coordintes of lesion center
wx = 1; wy = 0.9375; wz = wy; % voxel dimensions

% 2. Find centers of elements:
centers_triangles = mesh_get_triangle_centers(m2_seg); % Centers triangles
centers_tetrahedron = mesh_get_tetrahedron_centers(m2_seg); % Centers tetrahedrons 

% Filter on distance between voxel and element
% tetrahedrons
[idx_vox_tet,D] = knnsearch(lesion_coordinates,centers_tetrahedron); %index voxel closest to element, D = distance between element and closest voxel
A_tet = [idx_vox_tet, D];
A_tet(:,3) = 1:size(A_tet,1); %element index
A_tet(A_tet(:,2)>2,:) = []; % Only keep indeces with a distance shorter than 2mm

% Triangles
[idx_vox_tri,D] = knnsearch(lesion_coordinates,centers_triangles); %index voxel closest to element, D = distance between element and closest voxel
A_tri = [idx_vox_tri, D];
A_tri(:,3) = 1:size(A_tri,1); %element index
A_tri(A_tri(:,2)>2,:) = []; % Only keep indeces with a distance shorter than 2mm



% 3. Find elements within lesion voxels
% tetrahedrons
for i = 1:length(A_tet);
    if centers_tetrahedron((A_tet(i,3)),1) > lesion_coordinates(A_tet(i,1),1) - wx/2 && centers_tetrahedron((A_tet(i,3)),1) < lesion_coordinates(A_tet(i,1),1) + wx/2 && ...
            centers_tetrahedron((A_tet(i,3)),2) > lesion_coordinates(A_tet(i,1),2) - wy/2 && centers_tetrahedron((A_tet(i,3)),2) < lesion_coordinates(A_tet(i,1),2) + wy/2 && ...
            centers_tetrahedron((A_tet(i,3)),3) > lesion_coordinates(A_tet(i,1),3) - wz/2 && centers_tetrahedron((A_tet(i,3)),3) < lesion_coordinates(A_tet(i,1),3) + wz/2 
        K_tet(i,1) = A_tet(i,3); % gives indices of lesion elements
    end
end
K_tet(K_tet==0) = []; % delete zero's

% Triangles
for i = 1:length(A_tri);
    if centers_triangles((A_tri(i,3)),1) > lesion_coordinates(A_tri(i,1),1) - wx/2 && centers_triangles((A_tri(i,3)),1) < lesion_coordinates(A_tri(i,1),1) + wx/2 && ...
            centers_triangles((A_tri(i,3)),2) > lesion_coordinates(A_tri(i,1),2) - wy/2 && centers_triangles((A_tri(i,3)),2) < lesion_coordinates(A_tri(i,1),2) + wy/2 && ...
            centers_triangles((A_tri(i,3)),3) > lesion_coordinates(A_tri(i,1),3) - wz/2 && centers_triangles((A_tri(i,3)),3) < lesion_coordinates(A_tri(i,1),3) + wz/2 
        K_tri(i,1) = A_tri(i,3); % gives indices of lesion elements
    end
end
K_tri(K_tri==0) = []; % delete zero's

% 4. Replace corresponding labels with lesion label (10)
m2_seg.triangle_regions(K_tri)=1010; 
m2_seg.tetrahedron_regions(K_tet)=10; 

% Save modified mesh of Model 2
mesh_save_gmsh4(m2_seg, 'Data/401/401_6');
t=toc;
%% Electrode configuration
% Electrode configuration 1 (stimulate left MC)
configuration(1).e1 = 'C3'; % Location eletrode 1
configuration(1).e2 = 'Fp2';% Location eletrode 2

% Electrode configuration 2 (stimulate right MC)
configuration(2).e1='C4';% Location eletrode 1
configuration(2).e2='Fp1';% Location eletrode 2

% General stimulation settings
S = sim_struct('SESSION');
S.poslist{1} = sim_struct('TDCSLIST');
S.poslist{1}.currents = [0.002, -0.002];                % Current flow through each channel (mA)
    
% First Electrode 
S.poslist{1}.electrode(1).channelnr = 1;                % Connect the electrode to the first channel
S.poslist{1}.electrode(1).shape = 'ellipse';            % Elliptical shape
S.poslist{1}.electrode(1).dimensions = [10, 10];        % Dimension in mm
S.poslist{1}.electrode(1).thickness = 3;                % 3 mm thickness
    
% Second Electrode
S.poslist{1}.electrode(2).channelnr = 2;                % Connect the electrode to the second channel
S.poslist{1}.electrode(2).shape = 'ellipse';            % Elliptical shape
S.poslist{1}.electrode(2).dimensions = [10, 10];        % Dimension in mm
S.poslist{1}.electrode(2).thickness = 3;                % 3 mm thickness

%% tDCS simulations 
 for j = 1:length(configuration) % For every configuration tDCS simulation
    S.poslist{1}.electrode(1).centre = configuration(j).e1;  % Location electrode 1             
    S.poslist{1}.electrode(2).centre = configuration(j).e2;  % Location electrode 2
    folder_name = ['Data/401/', num2str(configuration(j).e1), '_', num2str(configuration(j).e2),'_simulation'];
    mkdir(folder_name); % Make folder for output
    
    %% Calculations 
    % Model 1
    S.fnamehead = 'Data/401/401_5.msh';                              % Mesh to simulate stimulation
    S.pathfem = ['Data/401/', num2str(configuration(j).e1), '_', num2str(configuration(j).e2),'_simulation/simulation_5']; % Set path for the simulation output
         
%       run_simnibs(S)                                        % Run the simulation

    % Load simulation data
    temp = mesh_load_gmsh4(['Data/401/', num2str(configuration(j).e1), '_', num2str(configuration(j).e2),'_simulation/simulation_5/401_5_TDCS_1_scalar.msh']);           % load simulation of Model 1
    temp.max = []; temp.max_index = []; temp.max_perc= [];
    m1.configuration(j) = temp; clear temp; 
    
    % Model 2
    S.fnamehead = 'Data/401/401_6.msh';                                         % Mesh to simulate stimulation
    cond_lesion = 126:100:1654;
    LesionConductivity = 0.126:0.1:1.654;
    for i = 1:length(cond_lesion)
        S.pathfem =sprintf('Data/401/%s_%s_simulation/simulation_6_%d',configuration(j).e1, configuration(j).e2, cond_lesion(i));        % Folder for the simulation output   %.2f', cond_lesion(i));                            
        S.poslist{1,1}.cond(10).value = LesionConductivity(i);                                 % Lesion conductivity
        S.subpath = 'Data/401/m2m_401_5';
    
%           run_simnibs(S)                                % Run the simulation
    end

    % Load simulation data
    for i = 1:length(cond_lesion)
        m2(j).configuration(i).mesh_name = sprintf('m2_%d', cond_lesion(i));
        file_name = sprintf('Data/401/%s_%s_simulation/simulation_6_%d/401_6_TDCS_1_scalar.msh',configuration(j).e1, configuration(j).e2, cond_lesion(i));
        m2(j).configuration(i).output = mesh_load_gmsh4(file_name);
    end
    
%% Analysis - General results
% Range E in GM
for i = 1:length(cond_lesion)
    max_E(j).configuration(i)=max(m2(j).configuration(i).output.element_data{2, 1}.tetdata(m2(j).configuration(i).output.tetrahedron_regions==2));
  	min_E(j).configuration(i)=min(m2(j).configuration(i).output.element_data{2, 1}.tetdata(m2(j).configuration(i).output.tetrahedron_regions==2));
end
max_E_tot(j)=max(max_E(j).configuration)
min_E_tot(j)=min(min_E(j).configuration);
%% Result Analysis - Whole GM Volume 
% Max E in GM of Model 1
 m1.configuration(j).max = max(m1.configuration(j).element_data{2, 1}.tetdata(m1.configuration(j).tetrahedron_regions==2)); % max E 
 m1.configuration(j).max_index = find(m1.configuration(j).element_data{2, 1}.tetdata==m1.configuration(j).max); % Index of max E
 % Max E in GM of Model 2
 for i = 1:length(LesionConductivity)
    m2(j).configuration(i).max = max(m2(j).configuration(i).output.element_data{2, 1}.tetdata(m2(j).configuration(i).output.tetrahedron_regions==2)); % max  E 
    m2(j).configuration(i).max_index = find(m2(j).configuration(i).output.element_data{2, 1}.tetdata==m2(j).configuration(i).max); % Index of max  E
    m2(j).configuration(i).max_perc = (m2(j).configuration(i).max - m1.configuration(j).max)./m1.configuration(j).max*100;
 end
   
     %% Analysis - Target area  
    %% Find Motor Cortex (MC) 
%     mesh_show_surface(m2(j).configuration(1).output); end % Shows only gray matter
%     
% %    Click on motor cortex
%     dcm_obj(j) = datacursormode(gcf);
%     k = 0;
%     while k == 0
%         set(dcm_obj(j),'DisplayStyle','window','SnapToDataVertex','on','Enable','on');
%         k = waitforbuttonpress;
%     end
%     c_info.configuration(j) = getCursorInfo(dcm_obj(j));
%     MC_center(j).configuration = c_info.configuration(j).Position;    % coordinates of center of sphere of MC
%%   Pointed by Medical Doctor 
    MC_center(1).configuration= [-19.095 21.61 33.0];   % Ipsilesional motor cortex  
    MC_center(2).configuration= [50.02 3.593 19.85];    % Contralesional motor cortex
    
    % Make sphere representing the target area MC
    r =10;        % radius in mm of sphere
    
    % Find index triangles of sphere 
    D_MC_tri = pdist2(centers_triangles,MC_center(j).configuration,'euclidean'); %Distance between all tetrahedron centers and center of MC
    D_MC_tri(:,2) = 1:size(D_MC_tri,1); % Add tetrahedron index
    D_MC_tri(D_MC_tri(:,1)>r,:) = []; % Only keep indeces of tetrahedrons with a distance shorter than 10mm, so it is in the sphere
    for ii=1:length(D_MC_tri)
        if m2(j).configuration(1).output.triangle_regions(D_MC_tri(ii,2))==1002 % Only de GM within the sphere
           MC_tri_index(j).configuration(ii,1) = D_MC_tri(ii,2);
        end
    end
    MC_tri_index(j).configuration(MC_tri_index(j).configuration==0) = []; %delete zero's
    
    % Find index tetrahedrons of sphere 
    D_MC_tet = pdist2(centers_tetrahedron,MC_center(j).configuration,'euclidean'); %Distance between all tetrahedron centers and center of MC
    D_MC_tet(:,2) = 1:size(D_MC_tet,1); % Add tetrahedron index
    D_MC_tet(D_MC_tet(:,1)>r,:) = []; % Only keep indeces of tetrahedrons with a distance shorter than 10mm, so it is in the sphere
    for ii=1:length(D_MC_tet)
        if m2(j).configuration(1).output.tetrahedron_regions(D_MC_tet(ii,2))==2 % Only de GM within the sphere
           MC_tet_index(j).configuration(ii,1) = D_MC_tet(ii,2);
        end
    end
    MC_tet_index(j).configuration(MC_tet_index(j).configuration==0) = []; %delete zero's
 
    % Calculate mean and max E of the tetrahedrons in MC and in lesion
    % Model 1:
    MC_meanE_m1(j).configuration = mean(m1.configuration(j).element_data{2, 1}.tetdata(MC_tet_index(j).configuration));
    MC_maxE_m1(j).configuration = max(m1.configuration(j).element_data{2, 1}.tetdata(MC_tet_index(j).configuration));
    
    % Model 2:
    for i=1:length(cond_lesion)
        % Motor Cortex
        MC_meanE_m2(j).configuration(i) = mean(m2(j).configuration(i).output.element_data{2, 1}.tetdata(MC_tet_index(j).configuration));
        m2(j).configuration(i).MC_maxE = max(m2(j).configuration(i).output.element_data{2, 1}.tetdata(MC_tet_index(j).configuration));
        % Lesion
        lesion_tet_index(j).configuration=find(m2(j).configuration(i).output.element_data{2, 1}.tetdata(m2(j).configuration(i).output.tetrahedron_regions==10));
       
        m2(j).configuration(i).lesion_meanE = mean(m2(j).configuration(i).output.element_data{2, 1}.tetdata(K_tet));
        m2(j).configuration(i).lesion_maxE = max(m2(j).configuration(i).output.element_data{2, 1}.tetdata(K_tet));
       
        % Percentage difference of mean E (M2-M1./M1*100)
        m2(j).configuration(i).MC_dif=(MC_meanE_m2(j).configuration(i) - MC_meanE_m1(j).configuration) ./ MC_meanE_m1(j).configuration * 100; 
    end
 end
%% FIGURES
%% Fig. 3
% M1 Ipsilesional Stimulation
figure(1);
j=1;
mesh_show_surface(m1.configuration(j),'colormap', jet,'showElec', true, 'scaleLimits', [0 0.5]);%'scaleLimits', [0 max_E(1).configuration(end)]);
title('Model 1 Ipsilesional Stimulation');

% M1 Contralesional Stimulation
figure(2);
j=2;
mesh_show_surface(m1.configuration(j),'colormap', jet,'showElec', true, 'scaleLimits', [0 0.5]);
title('Model 1 Contralesional Stimulation');

% M2 Ipsilesional Stimulation with lesion conductivity of 0.50 S/m
figure(3);
j=1;
i=1;
mesh_show_surface(m2(j).configuration(i).output, 'colormap', jet,'showElec', true, 'scaleLimits', [0 0.5]);
title(['Model 2 Ipsilesion Stimulation with lesion conductivity of ', num2str(LesionConductivity(i)), ' [S/m]']); %Electric field strength, through GM 

% M2 Contralesional Stimulation with lesion conductivity of 0.50 S/m
figure(4);
j=2;
i=1;
mesh_show_surface(m2(j).configuration(i).output, 'colormap', jet,'showElec', true, 'scaleLimits', [0 0.5]);
title(['Model 2 Contralesion Stimulation with lesion conductivity of ', num2str(LesionConductivity(i)), ' [S/m]']);

% M2 Ipsilesional Stimulation with lesion conductivity of 2.00 S/m
figure(5);
j=1;
i=length(LesionConductivity);
mesh_show_surface(m2(j).configuration(i).output, 'colormap', jet,'showElec', true, 'scaleLimits', [0 0.5]);
title(['Model 2 Ipsilesion Stimulation with lesion conductivity of ', num2str(LesionConductivity(i)), ' [S/m]']);

% M2 Contralesional Stimulation with lesion conductivity of 2.00 S/m
figure(6);
j=2;
i=length(LesionConductivity);
mesh_show_surface(m2(j).configuration(i).output, 'colormap', jet, 'showElec', true, 'scaleLimits', [0 0.5]);
title(['Model 2 Contralesion Stimulation with lesion conductivity of ', num2str(LesionConductivity(i)), ' [S/m]']);
 %% Fig. 4 Visualization of the GM and the target areas
tet_centers_lesion = centers_tetrahedron(K_tet,:); % Lesion
tet_centers_MC_L = centers_tetrahedron(MC_tet_index(1).configuration,:);% Ipsilesional motor cortex  
tet_centers_MC_R = centers_tetrahedron(MC_tet_index(2).configuration,:);% contralesional motor cortex
centers_tetrahedron_m2(2).configuration = mesh_get_tetrahedron_centers(m2(2).configuration(1).output);
centers_tetrahedron_GM(2).configuration = centers_tetrahedron_m2(2).configuration(m2(2).configuration(1).output.tetrahedron_regions ==2,:);
tet_centers_GM =  centers_tetrahedron_GM(2).configuration; % GM

figure(7);
a = pointCloud(tet_centers_GM);
cmatrix = ones(size(a.Location)).*[0.9 0.9 0.9];
a = pointCloud(tet_centers_GM,'Color',cmatrix);
pcshow(a)
title('Visualization of the GM and the target areas');
hold on
scatter3( tet_centers_lesion (:,1), tet_centers_lesion (:,2), tet_centers_lesion (:,3));% Lesion
hold on
scatter3(tet_centers_MC_L(:,1),tet_centers_MC_L(:,2),tet_centers_MC_L(:,3),'r');% Ipsilesional motor cortex
hold on
scatter3(tet_centers_MC_R(:,1),tet_centers_MC_R(:,2),tet_centers_MC_R(:,3),'r');% contralesional motor cortex
hold on
trisurf(m2(j).configuration(i).output.triangles(m2(j).configuration(i).output.triangle_regions==1002,:),m2(j).configuration(i).output.nodes(:,1),m2(j).configuration(i).output.nodes(:,2),m2(j).configuration(i).output.nodes(:,3),ones(size(m2(j).configuration(i).output.nodes(:,3))),'Facecolor',[192/255 192/255 192/255], 'FaceAlpha', 0.25,'EdgeAlpha', 0.01) ;
axis off

%% Fig. 5 Electric field strength in motor cortex
% Ipsilesional Stimulation
figure(8);
j=1;
plot ( [min(LesionConductivity),max(LesionConductivity)],[MC_meanE_m1(j).configuration,MC_meanE_m1(j).configuration], 'color', [0 0.4470 0.7410], 'linewidth', 3); % straight line for the mean E in mc in Model 1
title('Electric field strength in motor cortex');xlabel('Lesion conductivity [S/m]');ylabel('Mean electric field strength [V/mm]');
hold on
plot(LesionConductivity, MC_meanE_m2(j).configuration,'.', 'MarkerSize', 40,'color', [0 0.4470 0.7410] );
set(gca,'Fontsize',32);
xticks(LesionConductivity)
     for i = 1:length(cond_lesion);
         xt(i) = LesionConductivity(1,i);
         yt(i)= MC_meanE_m2(j).configuration(1,i);
         str= ['  ', num2str(round(m2(j).configuration(i).MC_dif,1)), '%'];
         text(xt(i),yt(i),str, 'FontSize', 24)
     end
     
% Contralesional Stimulation
hold on
j=2;
plot ( [min(LesionConductivity),max(LesionConductivity)],[MC_meanE_m1(j).configuration,MC_meanE_m1(j).configuration], 'color', [0.8500 0.3250 0.0980],'linewidth', 3); % straight line for the mean E in mc in Model 1
title('Electric field strength in motor cortex');xlabel('Lesion conductivity [S/m]');ylabel('Mean electric field strength [V/m]');
hold on
plot(LesionConductivity, MC_meanE_m2(j).configuration,'.', 'MarkerSize', 40, 'color', [0.8500 0.3250 0.0980]);
set(gca,'Fontsize',32);
legend({['Ipsilesional: Model 1'],'Model 2',['Contralesional: Model 1'],'Model 2', 'Location', 'best'})
% xticks(LesionConductivity)
     for i = 1:length(cond_lesion)
         xt(i) = LesionConductivity(1,i);
         yt(i)= MC_meanE_m2(j).configuration(1,i);
         str= ['  ', num2str(round(m2(j).configuration(i).MC_dif,1)), '%'];
         text(xt(i),yt(i),str, 'FontSize', 14)
     end
     
%% Fig. 7
for j=1:2
    for i = 1:length(cond_lesion)
        m2(j).configuration(i).output.triangle_regions(MC_tri_index(j).configuration)=1011;
        m2(j).configuration(i).output.tetrahedron_regions(MC_tet_index(j).configuration)=11;
    end
end
% Model 1 Ipsi
figure(9);
j=1;
m1.configuration(j).triangle_regions(MC_tri_index(j).configuration)=1011;
m1.configuration(j).tetrahedron_regions(MC_tet_index(j).configuration)=11;
mesh_show_surface(m1.configuration(j),'region_idx',1011, 'colormap', jet, 'scaleLimits', [0 0.7]);
trisurf(m2(j).configuration(i).output.triangles(m2(j).configuration(i).output.triangle_regions==1002,:),m2(j).configuration(i).output.nodes(:,1),m2(j).configuration(i).output.nodes(:,2),m2(j).configuration(i).output.nodes(:,3),ones(size(m2(j).configuration(i).output.nodes(:,3))),'Facecolor',[192/255 192/255 192/255], 'FaceAlpha', 0.25,'EdgeAlpha', 0.01) ;
title('Model 1 Ipsilesional stimulation');
      
% Model 1 Contra
figure(10);
j=2;
m1.configuration(j).triangle_regions(MC_tri_index(j).configuration)=1011;
m1.configuration(j).tetrahedron_regions(MC_tet_index(j).configuration)=11;
mesh_show_surface(m1.configuration(j),'region_idx',1011, 'colormap', jet, 'scaleLimits', [0 0.7]);
trisurf(m2(j).configuration(i).output.triangles(m2(j).configuration(i).output.triangle_regions==1002,:),m2(j).configuration(i).output.nodes(:,1),m2(j).configuration(i).output.nodes(:,2),m2(j).configuration(i).output.nodes(:,3),ones(size(m2(j).configuration(i).output.nodes(:,3))),'Facecolor',[192/255 192/255 192/255], 'FaceAlpha', 0.25,'EdgeAlpha', 0.01) ;
title('Model 1 Contralesional stimulation');

% Model 2 Ipsi with lesion conductivity of 0.50 S/m
figure(11);
j=1;
i=1;
mesh_show_surface(m2(j).configuration(i).output,'region_idx',1010, 'colormap', jet, 'scaleLimits', [0 0.7]); %lesion
trisurf(m2(j).configuration(i).output.triangles(m2(j).configuration(i).output.triangle_regions==1002,:),m2(j).configuration(i).output.nodes(:,1),m2(j).configuration(i).output.nodes(:,2),m2(j).configuration(i).output.nodes(:,3),ones(size(m2(j).configuration(i).output.nodes(:,3))),'Facecolor',[192/255 192/255 192/255], 'FaceAlpha', 0.25,'EdgeAlpha', 0.01) ;% GM
hold on
mesh_show_surface(m2(j).configuration(i).output,'region_idx',1011,'colormap', jet, 'scaleLimits', [0 0.7]);% Motor Cortex
title(['Model 2 Ipsilesional stimulation with lesion conductivity ', num2str(LesionConductivity(i))]);

% Model 2 Contra with lesion conductivity of 0.50 S/m
figure(12);
j=2;
i=1;
mesh_show_surface(m2(j).configuration(i).output,'region_idx',1010, 'colormap', jet, 'scaleLimits', [0 0.7]); %lesion
trisurf(m2(j).configuration(i).output.triangles(m2(j).configuration(i).output.triangle_regions==1002,:),m2(j).configuration(i).output.nodes(:,1),m2(j).configuration(i).output.nodes(:,2),m2(j).configuration(i).output.nodes(:,3),ones(size(m2(j).configuration(i).output.nodes(:,3))),'Facecolor',[192/255 192/255 192/255], 'FaceAlpha', 0.25,'EdgeAlpha', 0.01) ;% GM
hold on
mesh_show_surface(m2(j).configuration(i).output,'region_idx',1011,'colormap', jet, 'scaleLimits', [0 0.7]);% Motor Cortex
title(['Model 2 Contralesional stimulation with lesion conductivity ', num2str(LesionConductivity(i))]);

% Model 2 Ipsi with lesion conductivity of 2.00 S/m
figure(13);
j=1;
i=LesionConductivity(end);
mesh_show_surface(m2(j).configuration(i).output,'region_idx',1010, 'colormap', jet, 'scaleLimits', [0 0.7]); %lesion
trisurf(m2(j).configuration(i).output.triangles(m2(j).configuration(i).output.triangle_regions==1002,:),m2(j).configuration(i).output.nodes(:,1),m2(j).configuration(i).output.nodes(:,2),m2(j).configuration(i).output.nodes(:,3),ones(size(m2(j).configuration(i).output.nodes(:,3))),'Facecolor',[192/255 192/255 192/255], 'FaceAlpha', 0.25,'EdgeAlpha', 0.01) ;% GM
hold on
mesh_show_surface(m2(j).configuration(i).output,'region_idx',1011,'colormap', jet, 'scaleLimits', [0 0.7]);% Motor Cortex
title(['Model 2 Ipsilesional stimulation with lesion conductivity ', num2str(LesionConductivity(i))]);

% Model 2 Contra with lesion conductivity of 2.00 S/m
figure(14);
j=2;
i=LesionConductivity(end);
mesh_show_surface(m2(j).configuration(i).output,'region_idx',1010, 'colormap', jet, 'scaleLimits', [0 0.7]); %lesion
trisurf(m2(j).configuration(i).output.triangles(m2(j).configuration(i).output.triangle_regions==1002,:),m2(j).configuration(i).output.nodes(:,1),m2(j).configuration(i).output.nodes(:,2),m2(j).configuration(i).output.nodes(:,3),ones(size(m2(j).configuration(i).output.nodes(:,3))),'Facecolor',[192/255 192/255 192/255], 'FaceAlpha', 0.25,'EdgeAlpha', 0.01) ;% GM
hold on
mesh_show_surface(m2(j).configuration(i).output,'region_idx',1011,'colormap', jet, 'scaleLimits', [0 0.7]);% Motor Cortex
title(['Model 2 Contralesional stimulation with lesion conductivity ', num2str(LesionConductivity(i))]);

%% TABLES
 
%% Table II
Table2.Total = [length(centers_tetrahedron);length(centers_triangles)];
Table2.Lesion =[length(K_tet);length(K_tri)];
Table2.MC_i =[length(MC_tet_index(1).configuration);length(MC_tri_index(1).configuration)];
Table2.MC_c = [length(MC_tet_index(2).configuration);length(MC_tri_index(2).configuration)];
Table2=struct2table(Table2)

%% Table III
% Model 1
Table3_M1.Max_E_Ipsi = m1.configuration(1).max;
Table3_M1.Max_E_Contra = m1.configuration(2).max;
Table3_Model1=struct2table(Table3_M1) 

% Model 2
 for i = 1:length(LesionConductivity)
    Table3_M2(i).Lesion_Conductivity = LesionConductivity(i);
    Table3_M2(i).Max_E_Ipsi = m2(1).configuration(i).max;
    Table3_M2(i).Max_E_Contra = m2(2).configuration(i).max;
 end 
Table3_Model2=struct2table(Table3_M2) 

%% Table IV
% Ipsilesional Stimulation
j=1;
for i = 1:length(LesionConductivity)
    Table4_ipsilesional(i).Lesion_Conductivity = LesionConductivity(i);
    Table4_ipsilesional(i).Max_E_MotorCotex = m2(j).configuration(i).MC_maxE;
    Table4_ipsilesional(i).Mean_E_MotorCortex =  MC_meanE_m2(j).configuration(i);        
    Table4_ipsilesional(i).Percentage_Dif_Mean_E_MotorCortex =   m2(j).configuration(i).MC_dif;        
    Table4_ipsilesional(i).Mean_E_Lesion = m2(j).configuration(i).lesion_meanE;
    Table4_ipsilesional(i).Max_E_Lesion =m2(j).configuration(i).lesion_maxE ;
end
     Table4_Ipsilesional=struct2table(Table4_ipsilesional)
     
% Contralesional Stimulation
j=2;
for i = 1:length(LesionConductivity)
    Table4_contralesional(i).Lesion_Conductivity = LesionConductivity(i);
    Table4_contralesional(i).Max_E_MotorCotex = m2(j).configuration(i).MC_maxE;
    Table4_contralesional(i).Mean_E_MotorCortex =  MC_meanE_m2(j).configuration(i);        
    Table4_contralesional(i).Percentage_Dif_Mean_E_MotorCortex =   m2(j).configuration(i).MC_dif;        
    Table4_contralesional(i).Mean_E_Lesion = m2(j).configuration(i).lesion_meanE;
    Table4_contralesional(i).Max_E_Lesion =m2(j).configuration(i).lesion_maxE ;
end
     Table4_Contralesional=struct2table(Table4_contralesional)
