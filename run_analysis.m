load([path2mri,'contra_mc.mat'])
load([path2mri, 'ipsi_mc.mat'])

r =10;        % radius in mm of sphere
max_diff = zeros(length(idx_subj),1);
max_grey = cell(length(idx_subj),1);
for subj=1:length(idx_subj)
    subjID = ['sub-0', num2str(idx_subj(subj))];
    path2linda = [path2mri, subjID,'/linda/Prediction3_native.nii.gz'];
    path2msh_folder = [path2mri, subjID,'/'];
    path2msh_file = [path2msh_folder, subjID,'.msh'];
    
    maxs = zeros(length(LesionConductivity)+1,2);
    rel_diff = zeros(length(LesionConductivity),2);
    max_grey_help = zeros(length(LesionConductivity)+1,2);

    % Make sphere representing the target area MC and compute max E norm there


    % ipsi%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mc_center_ipsi      = ipsi_mc(subj,:); %useless

    % nolesion
    results = mesh_load_gmsh4([path2mri, subjID, '/C3_Fp2_simulation/general_hm/', subjID,'_TDCS_1_scalar.msh']);
    centers_tetrahedron = mesh_get_tetrahedron_centers(results); % Centers tetrahedrons 

    % Find index tetrahedrons of sphere 
    D_MC_tet = pdist2(centers_tetrahedron,mc_center_ipsi,'euclidean'); %Distance between all tetrahedron centers and center of MC
    D_MC_tet(:,2) = 1:size(D_MC_tet,1); % Add tetrahedron index
    D_MC_tet(D_MC_tet(:,1)>r,:) = []; % Only keep indeces of tetrahedrons with a distance shorter than 10mm, so it is in the sphere
    MC_tet_index = zeros(length(D_MC_tet),1);
    for ii=1:length(D_MC_tet)
        if results.tetrahedron_regions(D_MC_tet(ii,2))==2 % Only de GM within the sphere
           MC_tet_index(ii,1) = D_MC_tet(ii,2);
        end
    end
    MC_tet_index(MC_tet_index==0) = []; %delete zero's
    if subj ~=13
        maxs(1,1)  = max(results.element_data{2,1}.tetdata(MC_tet_index));
        max_grey_help(1,1)  = max(results.element_data{2,1}.tetdata(results.tetrahedron_regions==2));
    end
    
    %yeslesion
    for ll =1:length(LesionConductivity)
        results = mesh_load_gmsh4([path2mri, subjID, '/C3_Fp2_simulation/lesion_hm_',num2str(LesionConductivity(ll)*1000),'/', subjID,'_lesion_TDCS_1_scalar.msh']);
        centers_tetrahedron = mesh_get_tetrahedron_centers(results); % Centers tetrahedrons 
        D_MC_tet = pdist2(centers_tetrahedron,mc_center_ipsi,'euclidean'); %Distance between all tetrahedron centers and center of MC
        D_MC_tet(:,2) = 1:size(D_MC_tet,1); % Add tetrahedron index
        D_MC_tet(D_MC_tet(:,1)>r,:) = []; % Only keep indeces of tetrahedrons with a distance shorter than 10mm, so it is in the sphere
        MC_tet_index = zeros(length(D_MC_tet),1);
        for ii=1:length(D_MC_tet)
            if results.tetrahedron_regions(D_MC_tet(ii,2))==2 % Only de GM within the sphere
               MC_tet_index(ii,1) = D_MC_tet(ii,2);
            end
        end
        MC_tet_index(MC_tet_index==0) = []; %delete zero's
        if subj ~=13
            maxs(1+ll,1)  = max(results.element_data{2,1}.tetdata(MC_tet_index));
            max_grey_help(1+ll,1)  = max(results.element_data{2,1}.tetdata(results.tetrahedron_regions==2));
        end
    end
    
    mesh_show_surface(results,'showSurface', true)
    hold on
    scatter3(centers_tetrahedron(MC_tet_index,1),centers_tetrahedron(MC_tet_index,2),centers_tetrahedron(MC_tet_index,3))
    title(num2str(subjID))
    saveas(gcf,[path2mri,'/figures/', subjID, '_target_ipsi.png'])

    % contra%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mc_center_contra    = contra_mc(subj,:);

    % nolesion
    results = mesh_load_gmsh4([path2mri, subjID, '/C4_Fp1_simulation/general_hm/', subjID,'_TDCS_1_scalar.msh']);
    centers_tetrahedron = mesh_get_tetrahedron_centers(results); % Centers tetrahedrons 

    % Find index tetrahedrons of sphere 
    D_MC_tet = pdist2(centers_tetrahedron,mc_center_contra,'euclidean'); %Distance between all tetrahedron centers and center of MC
    D_MC_tet(:,2) = 1:size(D_MC_tet,1); % Add tetrahedron index
    D_MC_tet(D_MC_tet(:,1)>r,:) = []; % Only keep indeces of tetrahedrons with a distance shorter than 10mm, so it is in the sphere
    MC_tet_index = zeros(length(D_MC_tet),1);
    for ii=1:length(D_MC_tet)
        if results.tetrahedron_regions(D_MC_tet(ii,2))==2 % Only de GM within the sphere
           MC_tet_index(ii,1) = D_MC_tet(ii,2);
        end
    end
    MC_tet_index(MC_tet_index==0) = []; %delete zero's
    maxs(1,2)  = max(results.element_data{2,1}.tetdata(MC_tet_index));
    max_grey_help(1,2)  = max(results.element_data{2,1}.tetdata(results.tetrahedron_regions==2));

    %yeslesion
    for ll =1:length(LesionConductivity)
        results = mesh_load_gmsh4([path2mri, subjID, '/C4_Fp1_simulation/lesion_hm_',num2str(LesionConductivity(ll)*1000),'/', subjID,'_lesion_TDCS_1_scalar.msh']);
        centers_tetrahedron = mesh_get_tetrahedron_centers(results); % Centers tetrahedrons 
        D_MC_tet = pdist2(centers_tetrahedron,mc_center_contra,'euclidean'); %Distance between all tetrahedron centers and center of MC
        D_MC_tet(:,2) = 1:size(D_MC_tet,1); % Add tetrahedron index
        D_MC_tet(D_MC_tet(:,1)>r,:) = []; % Only keep indeces of tetrahedrons with a distance shorter than 10mm, so it is in the sphere
        MC_tet_index = zeros(length(D_MC_tet),1);
        for ii=1:length(D_MC_tet)
            if results.tetrahedron_regions(D_MC_tet(ii,2))==2 % Only de GM within the sphere
               MC_tet_index(ii,1) = D_MC_tet(ii,2);
            end
        end
        MC_tet_index(MC_tet_index==0) = []; %delete zero's
        maxs(1+ll,2)  = max(results.element_data{2,1}.tetdata(MC_tet_index));
        max_grey_help(1+ll,2)  = max(results.element_data{2,1}.tetdata(results.tetrahedron_regions==2));

    end

    
    max_grey{subj} = max_grey_help;
    
    mesh_show_surface(results,'showSurface', true)
    hold on
    scatter3(centers_tetrahedron(MC_tet_index,1),centers_tetrahedron(MC_tet_index,2),centers_tetrahedron(MC_tet_index,3))
    title(num2str(subjID))
    saveas(gcf,[path2mri,'/figures/', subjID, '_target_contra.png'])
    
    % rel diff and plot maxs
    rel_diff(:,1) = 100*(maxs(2:end,1) - maxs(1,1))/maxs(1,1);
    rel_diff(:,2) = 100*(maxs(2:end,2) - maxs(1,2))/maxs(1,2);


    figure, 
    plot ( [min(LesionConductivity),max(LesionConductivity)],[maxs(1,1), maxs(1,1)], 'color', [0 0.4470 0.7410], 'linewidth', 3); % straight line for the mean E in mc in Model 1
    hold on,
    plot(LesionConductivity, maxs(2:end,1), '.', 'MarkerSize', 40, 'color', [0 0.4470 0.7410], 'linewidth', 3), 
    hold on, 
    plot ( [min(LesionConductivity),max(LesionConductivity)],[maxs(1,2), maxs(1,2)], 'color', [0.8500 0.3250 0.0980], 'linewidth', 3); % straight line for the mean E in mc in Model 1
    hold on,
    plot(LesionConductivity, maxs(2:end,2),'.', 'MarkerSize', 40,'color', [0.8500 0.3250 0.0980],'linewidth', 3), 

    set(gca,'Fontsize',32)

    % add the max and min rel diff
    [a,b] = max(rel_diff(:,1));
    str= ['  ', num2str(round(a)), '%'];
    text(LesionConductivity(b),maxs(b+1,1),str, 'FontSize', 32)

    [a,b] = min(rel_diff(:,1));
    str= ['  ', num2str(round(a)), '%'];
    text(LesionConductivity(b),maxs(b+1,1),str, 'FontSize', 32)

    [a,b] = max(rel_diff(:,2));
    str= ['  ', num2str(round(a)), '%'];
    text(LesionConductivity(b),maxs(b+1,2),str, 'FontSize', 32)

    [a,b] = min(rel_diff(:,2));
    str= ['  ', num2str(round(a)), '%'];
    text(LesionConductivity(b),maxs(b+1,2),str, 'FontSize', 32)
    % text(LesionConductivity(b),0.51,str, 'FontSize', 32)


    % axis([min(LesionConductivity) max(LesionConductivity) 0.25 0.75])
    xlabel('Lesion conductivity value (S/m)');
    ylabel('Maximum electric field strength (V/m)');
    set(gcf, 'Position', get(0, 'Screensize'));
    title(num2str(subjID))
    saveas(gcf,[path2mri,'/figures/', subjID, '.png'])
    
    % save max diff in ipsi
    max_diff(subj) = max(abs(rel_diff(:,1)));
    
    
end

save max_diff max_diff
save max_grey max_grey