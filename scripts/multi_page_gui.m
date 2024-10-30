function multi_page_gui()
    close all; clear; clc;
    disp('Starting GUI application.')
    disp(' ');

    % Set up filepaths
    mainPath = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(mainPath, '..', 'src')));

    % Global variables
    head_model = [];
    pcMri_model = [];

    head_scan = [];
    prealigned_scan = [];
    filtered_scan = [];
    aligned_scan = [];


    % Create figure
    fig = uifigure('Name', 'Multi-Page SPOT3D GUI', 'Position', [100, 100, 800, 600]);

    % Create tab group for pages
    tabGroup = uitabgroup(fig, 'Position', [20, 20, 760, 560]);



    % ---- First Page: Import MRI and Create Head Model ----
    tab1 = uitab(tabGroup, 'Title', 'Create Head Model');
    
    % Import MRI file
    uilabel(tab1, 'Position', [20, 500, 100, 22], 'Text', 'Load MRI (.nii):');
    niiPathField = uieditfield(tab1, 'text', 'Position', [120, 500, 400, 22]);
    niiButton = uibutton(tab1, 'push', 'Text', 'Browse', 'Position', [530, 500, 80, 22], ...
                         'ButtonPushedFcn', @(~, ~) loadFile(niiPathField, '*.nii'));

    % Create and display head model button
    createModelButton = uibutton(tab1, 'push', 'Text', 'Create 3D Head Model', ...
                                 'Position', [20, 450, 150, 30], ...
                                 'ButtonPushedFcn', @(~, ~) createHeadModel(niiPathField.Value, tab1));

    % Save head model
    saveHeadModelButton = uibutton(tab1, 'push', 'Text', 'Save Head Model', ...
                                   'Position', [200, 450, 110, 30], ...
                                   'Enable', 'off', ...
                                   'ButtonPushedFcn', @(~, ~) saveHeadModel());

    % Axes for displaying head model
    niiAxes = uiaxes(tab1, 'Position', [120, 20, 500, 400]);
    niiAxes.Title.String = '3D Head Model';



    % ---- Second Page: Prepare for Alignment ----
    tab2 = uitab(tabGroup, 'Title', 'Prepare Alignment');

    % Import head scan
    uilabel(tab2, 'Position', [20, 500, 100, 22], 'Text', 'Load Scan (.stl):');
    stlPathField = uieditfield(tab2, 'text', 'Position', [120, 500, 400, 22]);
    stlButton = uibutton(tab2, 'push', 'Text', 'Browse', 'Position', [530, 500, 80, 22], ...
                         'ButtonPushedFcn', @(~, ~) loadScan(stlPathField, '*.stl'));

    % Prealign and Filter Scan buttons
    prealignButton = uibutton(tab2, 'push', 'Text', 'Prealign Scan', 'Position', [20, 450, 100, 30], ...
                              'ButtonPushedFcn', @(~, ~) prealignScan(stlPathField.Value, tab2));
    filterButton = uibutton(tab2, 'push', 'Text', 'Filter Scan', 'Position', [150, 450, 100, 30], ...
                            'Enable', 'off', 'ButtonPushedFcn', @(~, ~) filterScan(stlPathField.Value, tab2));

    % Axes for displaying pre-aligned model
    prealignAxes = uiaxes(tab2, 'Position', [120, 20, 500, 400]);
    prealignAxes.Title.String = 'Pre-aligned MRI and Scan';
    
    
    
    % ---- Third Page: Final Alignment ----
    tab3 = uitab(tabGroup, 'Title', 'Final Alignment');

    % Align button and save aligned scan
    alignButton = uibutton(tab3, 'push', 'Text', 'Align', 'Position', [20, 500, 100, 30], ...
                           'ButtonPushedFcn', @(~, ~) finalAlign());
    saveAlignedScanButton = uibutton(tab3, 'push', 'Text', 'Save Aligned Scan', 'Position', [150, 500, 150, 30], ...
                                     'Enable', 'off', ...
                                     'ButtonPushedFcn', @(~, ~) saveAlignedScan());

    % Axes for displaying final alignment
    finalAlignAxes = uiaxes(tab3, 'Position', [120, 20, 500, 400]);
    finalAlignAxes.Title.String = 'Final Aligned Model';
    
    

    % ---- Fourth Page: Electrode Registration ----
    tab4 = uitab(tabGroup, 'Title', 'Electrode Registration');

    % Import electrodes
    uilabel(tab4, 'Position', [20, 500, 100, 22], 'Text', 'Load Electrodes (.elc):');
    
    


    % ---- Nested functions for each operation ----
    function loadFile(field, filter)
        [file, path] = uigetfile(filter);
        if isequal(file, 0)
            return;
        end
        field.Value = fullfile(path, file);
    end

    % ---- Tab1: Create Head Model ----

    function createHeadModel(mriPath, parentTab)
        % Load MRI and create 3D head model using your function
        if ~isempty(mriPath)
            [mri_t1, T_matrix] = load_files(mriPath);
            
            disp('Create 3D Head Model'); disp(' ');
            disp('Creating the 3D head model from MRI... Please wait');
            head_model = create_3d_model(mri_t1, T_matrix);  % Creates a point cloud of the MRI scan
            disp('Finished creating the final head model.');
            disp(' ');
            
            saveHeadModelButton.Enable = 'on'; % Enable filter button after prealignment
            visualizeHeadModel(head_model, niiAxes);
        end
    end

    function saveHeadModel()
        save_stlfile_user(head_model, 'Save Head Model As');
        disp(' ');
    end
    

    % ---- Tab2: Prepare Alignment ----
    function loadScan(field, filter)
        loadFile(field, filter);
        head_scan = stlread(field.Value);
        visualizeScan(pointCloud(head_scan.Points), prealignAxes, 'Head scan');
    end

    function prealignScan(scanPath, parentTab)
        if ~isempty(scanPath)
            prealigned_scan = prealign_scan(head_model);
            filterButton.Enable = 'on'; % Enable filter button after prealignment
            visualizePrealigned(prealigned_scan, pcMri_model, prealignAxes);
        end
    end

    function prealigned_scan = prealign_scan(head_model)
        disp('Starting ICP registration...'); 
    
        pcMri_model = pointCloud(head_model.Points);    % Loads MRI model and scan point cloud
        pcHead_scan = pointCloud(head_scan.Points);
        
        pcHead_scan_translated = prepare_scan_for_icp(pcMri_model, pcHead_scan);    % Does preregistration
        pcSphere = create_sphere(pcMri_model);
        plot_pair(pcSphere, pcHead_scan, pcMri_model);
        alignedToSphere = register_to_sphere(pcHead_scan_translated, pcSphere);     % Pre-registers onto a sphere
        
        prealigned_scan = user_rotate(alignedToSphere, pcMri_model);     % Makes user pre-align the scan
        
        disp('Finished ICP registration.')
    end

    function filterScan(scanPath, parentTab)
        if ~isempty(scanPath)
            % Call your filtering function here
            filtered_scan_temp  = user_filter_gui(prealigned_scan, parentTab);
            if ~isempty(filtered_scan_temp)
                filtered_scan = filtered_scan_temp;
                fig = ancestor(parentTab, 'figure');  % Get the main figure from the tab
                fig.Name = 'Filtered scan';     % Set the new title
                visualizeScan(filtered_scan, prealignAxes, 'Filtered scan');
            end
        end
    end


    % ---- Tab3: Final Alignment ----
    function finalAlign()
        % Call your final alignment function here
        aligned_scan = register_to_mri(filtered_scan, pcMri_model);      % Does final registraion onto MRI ptCloud
        saveAlignedScanButton.Enable = 'on';
        visualizeFinalAlignment(aligned_scan, pcMri_model, finalAlignAxes);
    end

    function saveAlignedScan()
        save_stlfile_user(aligned_scan, 'Aligned Scan');
    end

    function colors = computeDistances(scan, model)
        % Compute distances
        distances = pdist2(model.Location, scan.Location, 'euclidean', 'Smallest', 1);
        
        % Normalize distances to [0, 1]
        normalizedDistances = (distances - min(distances)) / (max(distances) - min(distances));
        
        % Map to colormap (e.g., 'jet')
        colormapJet = jet(256);
        colorIndices = round(normalizedDistances * 255) + 1;
        disp('Indices');
        disp(size(colorIndices))
        disp('Scan');
        disp(size(scan.Location))
        % Ensure colorIndices matches the number of points
        % numPoints = size(aligned_scan.Location, 1);
        % colorIndices = min(colorIndices, numel(colormapJet));  % Limit indices to colormap size
        
        % Generate RGB colors array
        colors = colormapJet(colorIndices, :);  % Match color rows to points
        size(colors)
    end

    % ---- Visualization functions ----
    function visualizeHeadModel(head_model, ax)
        cla(ax);

        trimesh(head_model, 'Parent', ax);
        title(ax, '3D Head Model');
        axis(ax, 'equal');
        colormap(ax, 'gray'); 
        view(ax, 3);
    end

    function visualizeScan(scan, ax, plot_title)
        cla(ax);

        pcshow(scan.Location, [0.6350 0.0780 0.1840], 'Parent', ax);
        rotate3d(ax, 'off')
        rotate3d(ax, 'on')
        title(ax, plot_title)
        axis(ax, 'equal');
        colormap(ax, 'gray');
        view(ax, 3);
        grid(ax, 'on');
        set(ax, 'Color', 'white', ...
                'XColor', [0.15 0.15 0.15], ...
                'YColor', [0.15 0.15 0.15], ...
                'ZColor', [0.15 0.15 0.15])

        fig = ancestor(ax, 'figure');
        fig.Color = 'white';
    end

    function visualizePrealigned(scan, model, ax)
        cla(ax);

        hold(ax, 'on');
        pcshow(scan.Location, [0.6350 0.0780 0.1840], 'Parent', ax, 'MarkerSize', 40);
        pcshow(model.Location, 'Parent', ax, 'MarkerSize', 40);
        title(ax, 'Pre-aligned MRI and Scan');
        axis(ax, 'equal');
        colormap(ax, 'gray');
        view(ax, 3);
        set(ax, 'Color', 'white', ...
                'XColor', [0.15 0.15 0.15], ...
                'YColor', [0.15 0.15 0.15], ...
                'ZColor', [0.15 0.15 0.15])
        
        fig = ancestor(ax, 'figure');
        fig.Color = 'white';
        
        hold(ax, 'off');
    end
    
    function visualizeFinalAlignment(aligned_scan, model, ax)
        cla(ax);
        
        hold(ax, 'on');
        colors = computeDistances(aligned_scan, model);
        pcshow(pointCloud(aligned_scan.Location, 'Color', uint8(colors * 255)), 'Parent', ax, 'MarkerSize', 20);
        pcshow(model, 'Parent', ax, 'MarkerSize', 10);
        title(ax, 'Final Aligned Model');
        axis(ax, 'equal');
        colormap(ax, 'gray');
        view(ax, 3);
        set(ax, 'Color', 'white', ...
                'XColor', [0.15 0.15 0.15], ...
                'YColor', [0.15 0.15 0.15], ...
                'ZColor', [0.15 0.15 0.15])

        fig = ancestor(ax, 'figure');
        fig.Color = 'white';

        hold(ax, 'off');
    end

end
