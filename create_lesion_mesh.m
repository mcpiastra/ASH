function [lesion_hm, idx_lesion_tet] = create_lesion_mesh(path2msh_file,path2linda,subjID,path2mri)

%% load
general_hm = mesh_load_gmsh4(path2msh_file);
lesion_mask = ft_read_mri(path2linda);

%% create lesion_hm
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Find coordinates of lesion mask
[r,c,v] = ind2sub(size(lesion_mask.anatomy),find(lesion_mask.anatomy == 1)); % get the indices of lesion voxels
lesion_vox = [r c v];                                               
lesion_coord = lesion_mask.transform*[lesion_vox ones(length(lesion_vox),1)]'; % from voxel indices to head-coordinates
lesion_coordinates = lesion_coord(1:3,:)';

% wx = 1; wy = 0.9375; wz = wy; % voxel dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Find centers of elements
centers_triangles = mesh_get_triangle_centers(general_hm); % Centers triangles
centers_tetrahedron = mesh_get_tetrahedron_centers(general_hm); % Centers tetrahedrons 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Filter on distance between voxel and element
% Tetrahedrons
[idx_vox_tet,D] = knnsearch(lesion_coordinates,centers_tetrahedron); %index voxel closest to element, D = distance between element and closest voxel
A_tet = [idx_vox_tet, D];
A_tet(:,3) = 1:size(A_tet,1); %element index
A_tet(A_tet(:,2)>2,:) = []; % Only keep indeces with a distance shorter than 2mm 
wx=1;
wy=1;
wz=1;

% Triangles
[idx_vox_tri,D] = knnsearch(lesion_coordinates,centers_triangles); %index voxel closest to element, D = distance between element and closest voxel
A_tri = [idx_vox_tri, D];
A_tri(:,3) = 1:size(A_tri,1); %element index
A_tri(A_tri(:,2)>2,:) = []; % Only keep indeces with a distance shorter than 2mm

% 4. Find elements within lesion voxels
% tetrahedrons
for i = 1:length(A_tet)
    if centers_tetrahedron((A_tet(i,3)),1) > lesion_coordinates(A_tet(i,1),1) - wx/2 && centers_tetrahedron((A_tet(i,3)),1) < lesion_coordinates(A_tet(i,1),1) + wx/2 && ...
            centers_tetrahedron((A_tet(i,3)),2) > lesion_coordinates(A_tet(i,1),2) - wy/2 && centers_tetrahedron((A_tet(i,3)),2) < lesion_coordinates(A_tet(i,1),2) + wy/2 && ...
            centers_tetrahedron((A_tet(i,3)),3) > lesion_coordinates(A_tet(i,1),3) - wz/2 && centers_tetrahedron((A_tet(i,3)),3) < lesion_coordinates(A_tet(i,1),3) + wz/2 
        K_tet(i,1) = A_tet(i,3); % gives indices of lesion elements
    end
end
K_tet(K_tet==0) = []; % delete zero's

% Triangles
for i = 1:length(A_tri)
    if centers_triangles((A_tri(i,3)),1) > lesion_coordinates(A_tri(i,1),1) - wx/2 && centers_triangles((A_tri(i,3)),1) < lesion_coordinates(A_tri(i,1),1) + wx/2 && ...
            centers_triangles((A_tri(i,3)),2) > lesion_coordinates(A_tri(i,1),2) - wy/2 && centers_triangles((A_tri(i,3)),2) < lesion_coordinates(A_tri(i,1),2) + wy/2 && ...
            centers_triangles((A_tri(i,3)),3) > lesion_coordinates(A_tri(i,1),3) - wz/2 && centers_triangles((A_tri(i,3)),3) < lesion_coordinates(A_tri(i,1),3) + wz/2 
        K_tri(i,1) = A_tri(i,3); % gives indices of lesion elements
    end
end
K_tri(K_tri==0) = []; % delete zero's

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Replace corresponding labels with lesion label (11)
lesion_hm = general_hm;

% 6. Exclude skull (1004 and 4) and scalp (1005 and 5)
noss_tri = lesion_hm.triangle_regions(K_tri) ~= 1004 & lesion_hm.triangle_regions(K_tri) ~= 1005;
noss_tet = lesion_hm.tetrahedron_regions(K_tet) ~= 4 & lesion_hm.tetrahedron_regions(K_tet) ~= 5;

idx_lesion_tri = K_tri(noss_tri);
idx_lesion_tet = K_tet(noss_tet);

lesion_hm.triangle_regions(idx_lesion_tri)=1011; 
lesion_hm.tetrahedron_regions(idx_lesion_tet)=11; 


mesh_show_surface(lesion_hm)
hold on
scatter3(centers_tetrahedron(idx_lesion_tet,1),centers_tetrahedron(idx_lesion_tet,2),centers_tetrahedron(idx_lesion_tet,3))
title(num2str(subjID))
saveas(gcf,[path2mri,'figures/', subjID, '_lesion.png'])
figure,
scatter3(centers_tetrahedron(idx_lesion_tet,1),centers_tetrahedron(idx_lesion_tet,2),centers_tetrahedron(idx_lesion_tet,3))

toc; 