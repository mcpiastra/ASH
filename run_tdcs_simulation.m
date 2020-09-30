function run_tdcs_simulation(path2msh_folder, path2msh_file,subj, idx_subj)

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
    sim_folder = [path2msh_folder, num2str(configuration(j).e1), '_', num2str(configuration(j).e2),'_simulation'];
    mkdir(sim_folder); % Make folder for output
    
    %% Calculations 
    % General Head Model (no lesion)
    S.fnamehead = path2msh_file;            % Mesh to simulate stimulation
    S.pathfem = [sim_folder,'/general_hm/']; % Set path for the simulation output
         
    run_simnibs(S)                                        % Run the simulation

   
    % Lesion Head Model
    S.fnamehead = [path2msh_folder,'sub-0', num2str(idx_subj(subj)),'_lesion.msh'];                                         % Mesh to simulate stimulation
    cond_lesion = 126:100:1654;
    LesionConductivity = 0.126:0.1:1.654;
    for i = 1:length(cond_lesion)
        S.pathfem = [sim_folder,'/lesion_hm_',num2str(cond_lesion(i))];    % Folder for the simulation output                           
        S.poslist{1,1}.cond(11).value = LesionConductivity(i);             % Lesion conductivity
        S.subpath = [path2msh_folder,'/m2m_sub-0', num2str(idx_subj(subj))]; 
    
        run_simnibs(S)                                % Run the simulation

    end
 end



