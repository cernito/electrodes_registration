    function main_gui()
        
        % Clear matlab environment
        close all; clear; clc;
        all_fig = findall(groot, 'type', 'figure');
        close(all_fig)
    
        disp('Starting GUI application.')
        disp(' ');
        
        % Set up filepaths
        mainPath = fileparts(mfilename('fullpath'));
        dataPath = fullfile(mainPath,'..','data');
        resultsPath = fullfile(mainPath, '..', 'results');
        addpath(genpath(fullfile(mainPath, '..', 'src')));

        % Create results directory if it doesn't exist
        if ~exist(resultsPath, 'dir')
            mkdir(resultsPath);
        end
    
        % Global variables
        head_model = [];    % TR of head model created from MRI
        pcMri_model = [];   % PointCloud of that TR
    
        head_scan = [];              % TR of head scan
        head_scan_prealigned = [];   % TR of prealigned scan
        head_scan_fidu_aligned = []; % TR of scan aligned with fiducial lendmarks
        head_scan_aligned = [];      % TR of final aligned scan
        
        pcHead_scan = [];       % ptCloud of head scan
        prealigned_scan = [];   % ptCloud of prealigned and filtered scan
    
        fidu_aligned_scan = [];     % ptCloud of scan aligned with fiducial lendmarks
        aligned_scan = [];          % ptCloud of final aligned scan
        tform = {};                 % tform array holding all tforms made on head scan
        scan_fiducials = [];         % Coordinates of fiducial points on scan
    
        electrodes = [];                % electrodes struct read from elc_read
        pcElectrodes = [];              % ptCloud of electrode position
        pcElectrodes_aligned = [];      % ptCloud of aligned electrodes to scan (single position for electrode)
        pcElectrodes_refined = [];      % ptCloud of refined electrode positions
        pcElectrodes_templates = [];    % ptCloud of aligned tempalte electrodes (whole electrode ptClouds)
    
    
        % Define main figure dimensions
        figWidth = 1200;
        figHeight = 700;
        figXPos = 200;
        figYPos = 150;
        
        % Define padding and margins
        horizPadding = 20;
        vertPadding = 20;
        buttonHeight = 30;
        editFieldHeight = 22;
    
    
        
        % Create main GUI - UI figure
        GUIfig = uifigure('Name', 'Multi-Page Electrode Registration GUI', ...
                       'Position', [figXPos, figYPos, figWidth, figHeight], ...
                       'WindowButtonDownFcn', @(src,event)uiresume(src), ...
                       'WindowStyle','modal');
        
        % Create tab group for pages
        tabGroup = uitabgroup(GUIfig, 'Position', [horizPadding, vertPadding, ...
                                                figWidth - 2*horizPadding, ...
                                                figHeight - 2*vertPadding]);
        
    
    
        % ---- First Page: Import MRI and Create Head Model ----
        tab1 = uitab(tabGroup, 'Title', 'Create Head Model');
        
        % Import MRI file
        uilabel(tab1, 'Position', [horizPadding, figHeight-100, 100, editFieldHeight], ...
                'Text', 'Load MRI (.nii):');
        niiPathField = uieditfield(tab1, 'text', ...
                        'Position', [horizPadding+120, figHeight-100, figWidth-300, editFieldHeight]);
        niiButton = uibutton(tab1, 'push', 'Text', 'Browse', ...
                             'Position', [figWidth-150, figHeight-100, 80, editFieldHeight], ...
                             'ButtonPushedFcn', @(~, ~) loadFile(niiPathField, '*.nii'));
        
        % Create and display head model button
        createModelButton = uibutton(tab1, 'push', 'Text', 'Create 3D Head Model', ...
                                     'Position', [horizPadding, figHeight-150, 200, buttonHeight], ...
                                     'ButtonPushedFcn', @(~, ~) createHeadModel(niiPathField.Value));
        
        % Save head model
        saveHeadModelButton = uibutton(tab1, 'push', 'Text', 'Save Head Model', ...
                                       'Position', [horizPadding+230, figHeight-150, 170, buttonHeight], ...
                                       'Enable', 'off', ...
                                       'ButtonPushedFcn', @(~, ~) saveHeadModel());
        
        % Load head model
        loadHeadModelButton = uibutton(tab1, 'push', 'Text', 'Load Head Model', ...
                                       'Position', [horizPadding+430, figHeight-150, 170, buttonHeight], ...
                                       'Enable', 'on', ...
                                       'ButtonPushedFcn', @(~, ~) loadHeadModel());
        
        % Axes for displaying head model
        niiAxes = uiaxes(tab1, 'Position', [horizPadding+80, vertPadding, ...
                                            figWidth-2*horizPadding-160, figHeight-250]);
        niiAxes.Title.String = '3D Head Model';
        
    
    
        % ---- Second Page: Prepare for Alignment ----
        tab2 = uitab(tabGroup, 'Title', 'Prepare Alignment');
        
        % Import head scan
        uilabel(tab2, 'Position', [horizPadding, figHeight-100, 100, editFieldHeight], ...
                'Text', 'Load Scan (.stl):');
        stlPathField = uieditfield(tab2, 'text', ...
                        'Position', [horizPadding+120, figHeight-100, figWidth-300, editFieldHeight]);
        stlButton = uibutton(tab2, 'push', 'Text', 'Browse', ...
                             'Position', [figWidth-150, figHeight-100, 80, editFieldHeight], ...
                             'ButtonPushedFcn', @(~, ~) loadScan(stlPathField, '*.stl'));
        
        % Prealign and Filter Scan buttons
        prealignButton = uibutton(tab2, 'push', 'Text', 'Prealign Scan', ...
                                'Position', [horizPadding, figHeight-150, 150, buttonHeight], ...
                                'ButtonPushedFcn', @(~, ~) prealignScan(stlPathField.Value));
        filterButton = uibutton(tab2, 'push', 'Text', 'Filter Scan', ...
                                'Position', [horizPadding+180, figHeight-150, 150, buttonHeight], ...
                                'Enable', 'off', 'ButtonPushedFcn', @(~, ~) filterScan(stlPathField.Value, tab2));
        denoiseButton = uibutton(tab2, 'push', 'Text', 'Denoise Scan', ...
                                'Position', [horizPadding+360, figHeight-150, 150, buttonHeight], ...
                                'ButtonPushedFcn', @(~, ~) denoiseScan(stlPathField.Value, tab2));
        
        % Axes for displaying pre-aligned model
        prealignAxes = uiaxes(tab2, 'Position', [horizPadding+80, vertPadding, ...
                                                 figWidth-2*horizPadding-160, figHeight-250]);
        prealignAxes.Title.String = 'Pre-aligned MRI and Scan';




        % ---- Third Page: Final Alignment ----
        tab3 = uitab(tabGroup, 'Title', 'Final Alignment');
        
        % Align button and save aligned scan
        fiducialsButton = uibutton(tab3, 'push', 'Text', 'Select Fiducials', ...
                                'Position', [horizPadding, figHeight-100, 150, buttonHeight], ...
                                'ButtonPushedFcn', @(~, ~) selectFiducials());
        alignButton = uibutton(tab3, 'push', 'Text', 'Align', ...
                               'Position', [horizPadding+180, figHeight-100, 150, buttonHeight], ...
                               'Enable', 'off', 'ButtonPushedFcn', @(~, ~) finalAlign());
        warpingButton = uibutton(tab3, 'push', 'Text', 'Warp', ...
                                'Position', [horizPadding+360, figHeight-100, 150, buttonHeight], ...
                                 'Enable', 'off', ...
                                 'ButtonPushedFcn', @(~, ~) warpToMri());
        saveAlignedScanButton = uibutton(tab3, 'push', 'Text', 'Save Aligned Scan', ...
                                         'Position', [horizPadding+540, figHeight-100, 150, buttonHeight], ...
                                         'Enable', 'off', ...
                                         'ButtonPushedFcn', @(~, ~) saveAlignedScan());
        
        % Axes for displaying final alignment
        finalAlignAxes = uiaxes(tab3, 'Position', [horizPadding+80, vertPadding, ...
                                                   figWidth-2*horizPadding-160, figHeight-250]);
        finalAlignAxes.Title.String = 'Final Aligned Model';
        
    
        
    
        % ---- Fourth Page: Electrode Registration ----
        tab4 = uitab(tabGroup, 'Title', 'Electrode Registration');
        
        % Import electrodes (.elc, (taky .fcsv, .csv ?) )
        uilabel(tab4, 'Position', [horizPadding, figHeight-100, 130, editFieldHeight], ...
                'Text', 'Load Electrodes (.elc):');
        electrodesPathField = uieditfield(tab4, 'text', ...
                                'Position', [horizPadding+150, figHeight-100, figWidth-330, editFieldHeight]);
        electrodesButton = uibutton(tab4, 'push', 'Text', 'Browse', ...
                             'Position', [figWidth-150, figHeight-100, 80, editFieldHeight], ...
                             'ButtonPushedFcn', @(~, ~) loadElectrodes(electrodesPathField, '*.elc'));
        
        % Align and detect electrodes buttons
        alignElcButton = uibutton(tab4, 'push', 'Text', 'Align Electrodes', ...
                                'Position', [horizPadding, figHeight-150, 150, buttonHeight], ...
                                'Enable', 'off', ...
                                'ButtonPushedFcn', @(~, ~) alignElectrodes());
        alignTemplateButton = uibutton(tab4, 'push', 'Text', 'Align Template with Electrodes', ...
                                'Position', [horizPadding+180, figHeight-150, 250, buttonHeight], ...
                                'Enable', 'off', ...
                                'ButtonPushedFcn', @(~, ~) alignTemplateElectrodes());
        refineAlignmentButton = uibutton(tab4, 'push', 'Text', 'Refine Alignment', ...
                                'Position', [horizPadding+460, figHeight-150, 150, buttonHeight], ...
                                'Enable','off', ...
                                'ButtonPushedFcn', @(~,~) refineAlignment());
        templateMatchingButton = uibutton(tab4, 'push', 'Text', 'Template Matching', ...
                            'Position', [horizPadding+640, figHeight-150, 170, buttonHeight], ...
                            'Enable','off', ...
                            'ButtonPushedFcn', @(~,~) templateMatching());
    
        
        % Axes for displaying final alignment
        electrodesAxes = uiaxes(tab4, 'Position', [horizPadding+80, vertPadding, ...
                                                   figWidth-2*horizPadding-160, figHeight-250]);
        electrodesAxes.Title.String = 'Electrodes';




        % ---- Fifth Page: Electrode Revision ----
        tab5 = uitab(tabGroup, 'Title', 'Electrode Revision');
        
        % Axes for displaying the trisurf and electrodes
        revisionAxes = uiaxes(tab5, 'Position', [horizPadding+80, vertPadding, ...
                                               figWidth-2*horizPadding-160, figHeight-250]);
        revisionAxes.Title.String = 'Electrode Revision';
        axis(revisionAxes, 'equal');
        view(revisionAxes, 3);
        
        % Electrode revision buttons
        % Add a "Select Electrode" button
        selectElectrodeButton = uibutton(tab5, 'push', 'Text', 'Select Electrode', ...
                                 'Position', [horizPadding, figHeight-300, 120, buttonHeight], ...
                                 'ButtonPushedFcn', @(~, ~) startElectrodeSelection());
        refreshButton = uibutton(tab5, 'push', 'Text', 'Refresh Plot', ...
                                 'Position', [horizPadding, figHeight-250, 120, buttonHeight], ...
                                 'ButtonPushedFcn', @(~, ~) visualizeElectrodesForRevision());
        displayLabelsButton = uibutton(tab5, 'push', 'Text', 'Display Labels', ...
                         'Position', [horizPadding, figHeight-200, 120, buttonHeight], ...
                         'ButtonPushedFcn', @(~, ~) displayLabelsRevision());
        % Add an "Undo Move" button (optional, if you want to include it)
        undoMoveButton = uibutton(tab5, 'push', 'Text', 'Undo Move', ...
                         'Position', [horizPadding, figHeight-150, 120, buttonHeight], ...
                         'ButtonPushedFcn', @(~, ~) undoLastMove());

        % Initialize selected electrode index
        selectedElectrodeIdx = [];

        % Initialize selection mode flag
        isSelecting = false;

        % Initialize a stack to store previous positions (for Undo functionality)
        previousPositions = [];
    
    
    
    
    
        
        % ==========================================================================================
        % ==========================================================================================
        
        
        % ---- Nested functions for each operation ----
        function loadFile(field, filter)
            disp('Loading file...')
    
            [file, path] = uigetfile(filter);
            if isequal(file, 0)
                field.Value = [];
                return;
            end
            field.Value = fullfile(path, file);
            disp('File loaded.')
            disp(' ')
        end
    
        % ---- Tab1: Create Head Model ----
    
        function createHeadModel(mriPath)
            % Load MRI and create 3D head model using your function
            if ~isempty(mriPath)
                
                % Clear all plots
                %             cla(niiAxes);
                %             cla(prealignAxes);
                %             cla(finalAlignAxes);
                %             cla(electrodesAxes);
    
                [mri_t1, T_matrix] = load_files(mriPath);
                
                disp('Create 3D Head Model'); disp(' ');
                disp('Creating the .stl 3D head model from MRI... Please wait');
                
                % Creation of 3d head model
                head_model = create_3d_model(mri_t1, T_matrix);  % Creates a point cloud of the MRI scan
                
                assignin('base','head_model',head_model);
                disp('Finished creating the final head model.');
                disp(' ');
                
    
                disp('Visualizing head model... Please wait')
                % Visualization of created head model
                saveHeadModelButton.Enable = 'on'; % Enable filter button after prealignment
                visualizeHeadModel(head_model, niiAxes);
                disp('Completed visualization.')
                disp(' ')
    
                % Loads MRI model
                pcMri_model = pointCloud(head_model.Points);
            end
        end
    
        function saveHeadModel()
            disp('Saving head model as .stl file...')
            
            
            save_stlfile_user(head_model, 'Save Head Model As');
            
            disp('completed.')
            disp(' ');
        end
    
        function loadHeadModel()
            disp('Loading .stl head model...')
    
            head_model_tmp = load_stl();
            if isempty(head_model_tmp)
                disp('loading canceled.')
                disp(' ')
                return;
            end
            
            % Clear all plots
            %         cla(niiAxes);
            %         cla(prealignAxes);
            %         cla(finalAlignAxes);
            %         cla(electrodesAxes);
    
            head_model = head_model_tmp;
            visualizeHeadModel(head_model, niiAxes);
            pcMri_model = pointCloud(head_model.Points);
            
            disp('completed.')
            disp(' ')
        end
        
    
        % ---- Tab2: Prepare Alignment ----
        function loadScan(field, filter)
            loadFile(field, filter);
            if isempty(field.Value)
                return
            end
            disp('Loading new .stl scan file...')
    
            % Reset figure axes
            cla(prealignAxes);
            cla(finalAlignAxes);
            cla(electrodesAxes);
    
            % Reset scan variables
            prealigned_scan = [];       % Point Clouds
            fidu_aligned_scan = [];
            aligned_scan = [];
            
            tform = {}; 
            head_scan_prealigned = [];  % Triangulations
            head_scan_fidu_aligned = [];
            head_scan_aligned = [];
    
            head_scan = stlread(field.Value);
            pcHead_scan = pointCloud(head_scan.Points);
            
            disp('scan loaded.'); disp(' ');
    
    
            % Visualization
            disp('Visualizing loaded head scan... Please wait')
    
            %visualizeScan(pcHead_scan, prealignAxes, 'Head scan');
            visualizeTriangulation(head_scan, prealignAxes, 'Head scan');
            view(prealignAxes,3);
    
            alignButton.Enable = 'on'; % Enable alignment
            filterButton.Enable = 'on'; % Enable filter button
            alignTemplateButton.Enable = 'on';
            warpingButton.Enable = 'off';
    
            disp('Completed visualization.')
            disp(' ')
        end
    
        function prealignScan(scanPath)
            if ~isempty(scanPath)
                if isempty(head_model)
                    uialert(GUIfig, 'Head model not found. Create the 3D Head Model before performing alignment.', 'Alignment Error');
                    return
                end
    
                if ~isempty(prealigned_scan)
                    use_scan = prealigned_scan;
                    use_scanTR = head_scan_prealigned;
                else
                    use_scan = pcHead_scan;
                    use_scanTR = head_scan;
                end
    
                prealigned_scan = prealign_scan(use_scan);
                head_scan_prealigned = triangulation(use_scanTR.ConnectivityList, prealigned_scan.Location);
    
                disp('Completed.')
                disp(' ')
                
                disp('Visualizing the prealigned... Please wait');
                
                if ~isempty(head_model)
                    visualizeHeadModel(head_model,prealignAxes);
                    hold(prealignAxes,'on');
                else
                    % cla(prealignAxes);
                end
    
                visualizePrealigned(prealigned_scan, pcMri_model, prealignAxes);

                disp('Complteted visualization.'); disp(' ')
            end
        end
    
        function prealigned_scan = prealign_scan(use_scan)
            disp('Starting the prealignment...');
            
            tic
            % 1. PREALIGNMENT: Translate scan onto MRI centroid
            disp('1st PREALIGNMENT: Translating scan onto MRI centroid... Please wait')
            [pcHead_scan_translated, ~] = prepare_scan_for_icp(pcMri_model, use_scan);    % Does preregistration
            disp('Completed.')
        
            % 2. PREALIGNMENT: Align with a MRI sized sphere
            disp('2nd PREALIGNMENT: Aligning with a MRI sized sphere... Please wait')
            pcSphere = create_sphere(pcMri_model);
            plot_pair(pcSphere, use_scan, pcMri_model);
            [alignedToSphere, ~] = register_to_sphere(pcHead_scan_translated, pcSphere);     % Pre-registers onto a sphere
            disp('Completed.')
            
            % 3. PREALIGNMENT: Align principal axis
            disp('3nd PREALIGNMENT: Aligning main principal axis along z-axis... Please wait');
            alignedUpright = prealign_upright(alignedToSphere);
            disp('Completed.')
            toc

            % 4. PREALIGNMENT: Take user made prealignment
            disp('4rd PREALIGNMENT: Aligning the head scan manually')
            [prealigned_scan, ~] = user_rotate(alignedUpright, pcMri_model);     % Makes user pre-align the scan
            disp('Compelted.')
    
            disp('Finished ICP registration.'); disp(' ');
        end
    
        function filterScan(scanPath, parentTab)
            if ~isempty(scanPath)
                disp('Starting filtering of head scan.')
                disp('Choose filtering method and apply filter...')
                
                if ~isempty(prealigned_scan)
                    curr_scan = prealigned_scan;
                    curr_scanTR = head_scan_prealigned;
                else 
                    curr_scan = pcHead_scan;
                    curr_scanTR = head_scan;
                end
                
                % Calling filtering function
                [filtered_scan_temp, keepMask] = user_filter_gui(curr_scan, parentTab);
                
                if ~isempty(filtered_scan_temp)
                    prealigned_scan = filtered_scan_temp;
                    
                    % Filter the triangulation head scan
                    head_scan_prealigned = filterTriangulation(curr_scanTR, keepMask);
                    
                    if ~isempty(head_model)
                        visualizeHeadModel(head_model,prealignAxes);
                        hold(prealignAxes,'on');
                    else
                        % cla(prealignAxes);
                    end
                    
                    visualizePrealigned(prealigned_scan, pcMri_model, prealignAxes);
                    %visualizeTriangulation(head_scan_prealigned, prealignAxes, 'Filtered scan');
    
                    disp('Filter applied.');
                end
                disp('Completed filtering.'); disp(' ')
            end
        end
    
        function filteredTR = filterTriangulation(TR, keepMask)
            disp('Filtering head scan triangulation...')
            
            % Ensure keepMask is logical and within bounds
            keepMask = logical(keepMask);
    
            if numel(keepMask) > size(TR.Points,1)
                keepMask = keepMask(1:size(TR.Points,1));
            elseif numel(keepMask) < size(TR.Points,1)
                % Pad keepMask with false for any points beyond the mask
                tmpMask = false(size(TR.Points,1),1);
                tmpMask(1:numel(keepMask)) = keepMask;
                keepMask = tmpMask;
            end
            
            % Indices of points to keep
            keepIndices = find(keepMask);
            
            % Filter the vertex list
            filteredVertices = TR.Points(keepIndices, :);
        
            % Filter faces: only keep faces whose all vertices are in keepIndices
            validFacesMask = all(ismember(TR.ConnectivityList, keepIndices), 2);
            filteredFaces = TR.ConnectivityList(validFacesMask, :);
            
            % Remap old vertex indices to new indices
            [~, newIndexMap] = ismember(filteredFaces, keepIndices);
            filteredFaces = newIndexMap;
            
            % Create the new triangulation with updated faces and vertices
            filteredTR = triangulation(filteredFaces, filteredVertices);
            
            disp('Completed triangulation filtering.')
            disp(' ')
        end
    
        function denoiseScan(scanPath, parentTab)
            disp('denoising head scan...')
    
            if ~isempty(scanPath)
                if ~isempty(prealigned_scan)
                    curr_scan = prealigned_scan;
                    curr_scanTR = head_scan_prealigned;
                else
                    curr_scan = pcHead_scan;
                    curr_scanTR = head_scan;
                end
    
                keepMask = false(size(curr_scan.Location,1),1);
                [prealigned_scan,inlierIndices,~] = pcdenoise(curr_scan);
                
                keepMask(inlierIndices) = true;
                head_scan_prealigned = filterTriangulation(curr_scanTR, keepMask);
    
                fig = ancestor(parentTab, 'figure');  % Get the main figure from the tab
    
                if ~isempty(head_model)
                    visualizeHeadModel(head_model, prealignAxes);
                    hold(prealignAxes,'on');
                else
                    % cla(prealignAxes);
                end
    
                %visualizeScan(prealigned_scan, prealignAxes, 'Denoised scan');
                visualizePrealigned(prealigned_scan, pcMri_model, prealignAxes);
                %visualizeTriangulation(head_scan_prealigned, prealignAxes, 'Denoised scan')
                
                disp('Scan denoised'); 
                disp(' ');
            end
        end
    
    
        % ---- Tab3: Final Alignment ----
        function selectFiducials()
            
            if isempty(head_model)
                uialert(GUIfig, 'Import head model for fiducial selection.', 'Select Fiducial Error');
                return;
            end
            if isempty(head_scan)
                uialert(GUIfig, 'Import head scan for fiducial selection.', 'Select Fiducial Error');
                return;
            end
            
            disp('Starting fiducial selection.')
            
            if ~isempty(prealigned_scan)
                pcScan = prealigned_scan;
                scanTR = head_scan_prealigned;
            else
                pcScan = pcHead_scan;
                scanTR = head_scan;
            end
           
            % Create the main figure
            mainFig = uifigure('Name', 'Fiducial Selection', ...
                               'Position', [100, 100, 1600, 800]);
        
            % Create two axes within the figure
            ax1 = uiaxes(mainFig, 'Units', 'normalized', 'Position', [0.05 0.1 0.4 0.8]);
            ax1.Title.String = 'Head Scan';
        
            ax2 = uiaxes(mainFig, 'Units', 'normalized', 'Position', [0.55 0.1 0.4 0.8]);
            ax2.Title.String = 'Head Model';
            
            visualizeInFiducialPlot(head_model,ax2);
        
            % Call get_fiducials_stl on the first axes
            disp('Select fiducials on the head scan (Figure 1):');
            F2 = get_fiducials_stl(scanTR, ax1);
            disp('Fiducials chosen for head scan:');
            disp(F2);
    
            visualizeInFiducialPlot(scanTR,ax1);
        
            % Now call get_fiducials_stl on the second axes
            disp('Select fiducials on the head model (Figure 2):');
            F1 = get_fiducials_stl(head_model, ax2);
            disp('Fiducials chosen for head model:');
            disp(F1);
    
            
            % ========================
            % FIDUCIAL ALIGNMENT
            % ========================
    
            % ALIGNMENT: Transform based on fiducial lendmarks
            tform_fid = transformWithFiducials(F2,F1);
    
            % Apply the transformation to align ptCloud1 with ptCloud2
            fidu_aligned_scan = pctransform(pointCloud(pcScan.Location), tform_fid); % (R * pcScan.Location' + T)';
            
            % Align Triangulation head scan
            head_scan_fidu_aligned = triangulation(scanTR.ConnectivityList,fidu_aligned_scan.Location);
            
            % Save aligned scan to results folder
            savePath = fullfile(resultsPath, 'head_scan_fidu_aligned.mat');
            save(savePath, "head_scan_fidu_aligned");

            % Remove final aligned scan. 
            % Because we could want so use just fidu aligned scan as
            % aligned enough scan.
            head_scan_aligned = [];
            aligned_scan = [];
            
            % Transform Fiducial points
            fidu_ptcloud = pctransform(pointCloud(F2), tform_fid);
            scan_fiducials = fidu_ptcloud.Location;
    
            visualizeFiducials(F2,F1,tform_fid);
    
            % In case of bad performance of ICP, fiducial alignment can
            % be selected as sufficient alignment.
            saveAlignedScanButton.Enable = 'on';
            alignTemplateButton.Enable = 'on';
            warpingButton.Enable = 'on';
            
            
            % Compute RMSE
            points_static = pcMri_model.Location;
            points_moved = fidu_aligned_scan.Location;
        
            kdtree = KDTreeSearcher(points_static);
            nearestIdx = knnsearch(kdtree, points_moved);
            points_closest = points_static(nearestIdx,:);
        
            distance = sqrt(sum((points_closest - points_moved).^2, 2));
        
            % Compute the RMSE.
            rmse = sqrt(sum(distance.^2,'all') / numel(distance));
        
            fprintf('Fiducial alignment rmse:\n    %.4f\n', rmse);
            disp(' ');



            disp('Completed fiducial selection.')
            disp(' ')
            disp('Visualizing alignment results... Please wait')
    
            % Plot to Final Alignment Tab
            visualizeHeadModel(head_model,finalAlignAxes);
            hold(finalAlignAxes,'on');
            visualizeFinalAlignment(fidu_aligned_scan, pcMri_model, finalAlignAxes);
            hold(finalAlignAxes,'off');
            
            close(mainFig)
    
            % Plot to Prepare Alignment Tab
            % visualizeHeadModel(head_model, prealignAxes);
            % hold(prealignAxes,'on');
            % visualizeTriangulation(head_scan_fidu_aligned, prealignAxes, 'Filtered scan');
    
    
            disp('Finished visualization.')
            disp('Now perform Alignment.')
            disp(' ')
    
        end
    
        function tform_fid = transformWithFiducials(F2,F1)
            % F2: Moving fiducial points (scan)
            % F1: Fixed fiducial points (MRI head model)
            [tform_fid,~,~] = pcregistericp(pointCloud(F2),pointCloud(F1),"Metric","pointToPoint");
        end
    
        function simplifiedTR = simplifyTR(TR, reductionFactor)
            % reductionFactor: fraction of faces to keep (e.g., 0.1 for 10%)
            FV.faces = TR.ConnectivityList;
            FV.vertices = TR.Points;
        
            % Use reducepatch to simplify the mesh
            [faces, vertices] = reducepatch(FV, reductionFactor);
        
            simplifiedTR = triangulation(faces, vertices);
        end
    
        function finalAlign()
            % Final alignment onto MRI
            disp('Final alignment onto MRI... Please wait');
    
            if ~isempty(fidu_aligned_scan)
                use_scan = fidu_aligned_scan;
                use_scanTR = head_scan_fidu_aligned;
            elseif ~isempty(prealigned_scan)
                use_scan = prealigned_scan;
                use_scanTR = head_scan_prealigned;
            else
                use_scan = pcHead_scan;
                use_scanTR = head_scan;
            end
            
            % Does final registraion onto MRI ptCloud
            [aligned_scan, tform_aligned] = register_to_mri(use_scan, pcMri_model);
            prealigned_scan_aligned = pctransform(pointCloud(use_scanTR.Points), tform_aligned);
            head_scan_aligned = triangulation(use_scanTR.ConnectivityList, prealigned_scan_aligned.Location);
            
            disp('Completed alignment.'); disp(' ')
    
            % Assign to workspace
            assignin('base', 'head_scan_aligned', head_scan_aligned);
            
            % Save to results folder
            savePath = fullfile(resultsPath, 'head_scan_aligned.mat');
            save(savePath, "head_scan_aligned");           
            
            saveAlignedScanButton.Enable = 'on';
            alignTemplateButton.Enable = 'on';
            warpingButton.Enable = 'on';
            
            disp('Displaying results... Please wait')
            visualizeHeadModel(head_model,finalAlignAxes);
            hold(finalAlignAxes,'on');
            visualizeFinalAlignment(aligned_scan, pcMri_model, finalAlignAxes);
            hold(finalAlignAxes,'off')
            disp('Completed displaying results.'); disp(' ')
        end
    
        function combined_tform = getCombinedTransform(tform)
            % Create Final Transform Matrix - postmultiply
            combined_tform = eye(4);
            for i=length(tform):-1:1
                combined_tform = combined_tform * tform{i};
            end
        end
    
        function new_TR = transformTriangulation(TR, T)
            %TR_mesh = extendedObjectMesh(TR.Points, TR.ConnectivityList);
            %TR_mesh = applyTransform(TR_mesh, tform);
            
            points = TR.Points;
            homPoints = [points, ones(size(points,1),1)]';
            transformedPoints = T * homPoints;
            transformedPoints = transformedPoints(1:3,:)';
            new_TR = triangulation(TR.ConnectivityList, transformedPoints);
    
            %         R = tform(1:3,1:3);
            %         t = tform(1:3,4)';
            %         new_points = transformPointsForward(rigid3d(R,t),TR.Points);
        end
        
        function warpToMri()
            disp('Warping onto MRI... Please wait');

            if ~isempty(aligned_scan)
                use_scan = aligned_scan;
                use_scanTR = head_scan_aligned;
            elseif ~isempty(fidu_aligned_scan)
                use_scan = fidu_aligned_scan;
                use_scanTR = head_scan_fidu_aligned;
            elseif ~isempty(prealigned_scan)
                use_scan = prealigned_scan;
                use_scanTR = head_scan_prealigned;
            else
                use_scan = pcHead_scan;
                use_scanTR = head_scan;
            end
            
            [aligned_scan, ~] = warp_to_mri(use_scan, pcMri_model);
            head_scan_aligned = triangulation(use_scanTR.ConnectivityList, aligned_scan.Location);
            
            disp('Finished warping.'); disp(' ')
            

            % Visualization
            assignin('base','head_scan_aligned',head_scan_aligned);
            savePath = fullfile(resultsPath, 'head_scan_aligned.mat');
            save(savePath, "head_scan_aligned");

            saveAlignedScanButton.Enable = 'on';
            alignTemplateButton.Enable = 'on';
            
            disp('Displaying results... Please wait')
            visualizeHeadModel(head_model, finalAlignAxes);
            hold(finalAlignAxes,'on');
            visualizeFinalAlignment(aligned_scan, pcMri_model, finalAlignAxes);
            hold(finalAlignAxes,'off')
            disp('Completed displaying results.'); disp(' ')
        end

    
        function saveAlignedScan()
            if isempty(head_scan_aligned) && isempty(fidu_aligned_scan)
                uialert(GUIfig, 'Align the head scan first.', 'Save Aligned Error');
                return;
            end
            if ~isempty(head_scan_aligned)
                save_stlfile_user(head_scan_aligned, 'Save Aligned Scan as');
                verts = head_scan_aligned.Points;
                faces = head_scan_aligned.ConnectivityList;
                plyPath = fullfile(resultsPath, 'aligned_scan.ply');
                plywrite(plyPath,faces, verts);
            else
                save_stlfile_user(head_scan_fidu_aligned, 'Save Aligned Scan as');
                verts = head_scan_fidu_aligned.Points;
                faces = head_scan_fidu_aligned.ConnectivityList;
                plyPath = fullfile(resultsPath, 'aligned_scan.ply');
                plywrite(plyPath, faces, verts);
            end
        end
    
        function [colors, minDistance, maxDistance] = computeDistances(scan, model)
            
            indices = knnsearch(model.Location, scan.Location);
            distances = sqrt(sum((model.Location(indices, :) - scan.Location).^2, 2));
    
            %FV.faces = head_model.ConnectivityList;
            %FV.vertices = head_model.Points;
            % [distances, ~] = point2trimesh(FV, 'QueryPoints', scan.Location);
            maxDistance = max(distances);
            minDistance = min(distances);
            clampedDistances = min(distances, maxDistance);
            

            % Normalize distances to [0, 1]
            normalizedDistances = clampedDistances / maxDistance;
            
            % Map to colormap (e.g., 'jet')
            colormapJet = jet(256);
            colorIndices = round(normalizedDistances * 255) + 1;
            
            % Generate RGB colors array
            if all(isnan(colorIndices))
                colorIndices = ones(size(colorIndices)) * 128;
                colors = colormapJet(colorIndices, :);
            else
                colors = colormapJet(colorIndices, :);  % Match color rows to points
            end
        end
        
        % ---- tab4: Electrode Registration ----
        function loadElectrodes(field, filter)
            loadFile(field, filter);
            if isempty(field.Value)
                return
            end
            
            disp('Loading electrode .elc file...')
            [electrodes, pcElectrodes] = load_electrodes(field.Value);
            disp('File loaded.'); disp(' ')
    
            disp('Visualizing electrodes...')
            visualizeElectrodes(pcElectrodes, head_scan_aligned, electrodesAxes)
            disp('Finished visualization.'); disp(' ')
    
            alignElcButton.Enable = 'on';
        end
    
        function alignElectrodes()
            pcElectrodes_aligned = register_electrodes(electrodes, pcElectrodes, aligned_scan);
            visualizeElectrodes(pcElectrodes_aligned, head_scan_aligned, electrodesAxes)
        end
    
        function alignTemplateElectrodes()
            if isempty(head_scan_aligned)
                if ~isempty(head_scan_fidu_aligned)
                    head_scan_aligned = head_scan_fidu_aligned;
                    disp('Using scan aligned with fiducial lendmarks.')
                else
                    uialert(GUIfig, 'Align the head scan first.', 'Align Template Error');
                    return;
                end
            end
            
            % Reset electrode positions
            pcElectrodes_refined = [];
            
            
            disp('Alignment of template HD-EEG cap with head scan...')
            [pcElectrodes_aligned_tmp] = align_template_electrodes(head_scan_aligned,[],electrodesAxes); % [] = scan_fiducials
            if isempty(pcElectrodes_aligned_tmp)
                disp('User canceled the process.'); disp(' ')
                return;
            end
            pcElectrodes_aligned = pcElectrodes_aligned_tmp;
    
            assignin('base','pcElectrodes_aligned',pcElectrodes_aligned)
            disp('Finished alignment.'); disp(' ')
    
            disp('Visualizing results... Please wait')
            visualizeElectrodes(pcElectrodes_aligned, head_scan_aligned, electrodesAxes);
            
            % Save aligned electrode positions
            saveElectrodes();
    
            disp('Completed visualization.'); disp(' ')
            refineAlignmentButton.Enable = 'on';
            templateMatchingButton.Enable = 'on';
        end
    
        function refineAlignment()
            % electrodes_found = {}
            % electrodes_not_found = {}
            %
            % for each electrode
            %   select neighboring area in head_scan point cloud
            %   template match in that area to find a electrode
            %   if reasonable rmse set as 'found' otherwise 'not_found'
            %       store the found and not_found electrodes...
            %       ...for further visualization and possible manual refinement
            %   select the centroid of matched electrode as the electrode position
            %   store the electrode position
            % end
            %
            % return electrodes_found, electrodes_not_found
            if isempty(head_scan_aligned)
                if ~isempty(head_scan_fidu_aligned)
                    head_scan_aligned = head_scan_fidu_aligned;
                    disp('Using scan aligned with fiducial lendmarks.')
                else
                    uialert(GUIfig, 'Align the head scan first.', 'Refine Alignment Error');
                    return;
                end
            end

            assignin('base','scan',head_scan_aligned);
            assignin('base','pcElectrodes_aligned',pcElectrodes_aligned)
    
            set(GUIfig, 'HandleVisibility', 'off');
            pcElectrodes_refined = refine_alignment_with_gc(head_scan_aligned, []);  % [] = head_scan_prealigned [] = pcElectrodes_aligned
            set(GUIfig, 'HandleVisibility', 'on');
    
            visualizeElectrodes(pcElectrodes_refined, head_scan_aligned, electrodesAxes)
        end
        
        function templateMatching()
            disp('Starting template matching process...');
            
            if isempty(head_scan_aligned)
                uialert(GUIfig, 'Head scan must be aligned before template matching.', 'Template Matching Error');
                return;
            end
            
            [pcElectrodes_templates] = template_matching(head_scan_aligned, []); % [] = pcElectrodes_refined
            
            % Save results if needed
            assignin('base', 'pcElectrodes_templates', pcElectrodes_templates);
            disp('Template matching completed.');
            
            % Visualize results
            visualizeElectrodes(pcElectrodes_templates, head_scan_aligned, electrodesAxes);
            
            disp('Results visualized.');
        end
        
        
        function saveElectrodes()
            disp('Saving aligned template electrodes as .stl file...')
    
            electrodesFile = fullfile(resultsPath, 'electrodes.mat');
            save(electrodesFile, "pcElectrodes_aligned");
            
            disp('completed.')
            disp(' ');
        end

        % ---- tab5: Electrode Revision ----
        function visualizeElectrodesForRevision()
            cla(revisionAxes);
            
            if isempty(head_scan_aligned)
                uialert(GUIfig, 'Head scan is not aligned. Please align the head scan before revising electrodes.', 'Alignment Error');
                return;
            end
            if isempty(head_model)
                uialert(GUIfig, 'Head model does not exist. Please create or load the head model before revising electrodes.', 'Alignment Error');
                return;
            end
            
            % Plot the head model trisurf
            trisurf(head_model.ConnectivityList, head_model.Points(:,1), ...
                    head_model.Points(:,2), head_model.Points(:,3), ...
                    'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none', 'Parent', revisionAxes,...
                    'DisplayName', 'Scan');
            hold(revisionAxes, 'on');
            
            % Plot electrodes
            if ~isempty(pcElectrodes_refined)
                electrodesPos = pcElectrodes_refined.Location;
            elseif ~isempty(pcElectrodes_aligned)
                electrodesPos = pcElectrodes_aligned.Location;
            else
                electrodesPos = pcElectrodes.Location;
            end
            
            % Store electrode scatter object for interaction
            electrodeScatter = scatter3(electrodesPos(:,1), electrodesPos(:,2), electrodesPos(:,3), ...
                                        40, 'r', 'filled', 'MarkerEdgeColor', 'k', ...
                                'Parent', revisionAxes, 'DisplayName', 'Electrodes', ...
                                'HitTest', 'off', 'PickableParts', 'none');

            electrodesPos
            displayLabels(revisionAxes, pointCloud(electrodesPos));


            % Improve Visualization
            camlight(revisionAxes,'headlight');
            lighting(revisionAxes, 'gouraud');
            grid(revisionAxes, 'on');
            hold(revisionAxes, 'off');
            axis(revisionAxes, 'equal')

            
            % Add labels if necessary
            xlabel(revisionAxes, 'X');
            ylabel(revisionAxes, 'Y');
            zlabel(revisionAxes, 'Z');
            
            legend(revisionAxes, 'show');


            % Same code but in independent figure

            % figure
            % ax = gca;
            % % Plot the head model trisurf
            % trisurf(head_model.ConnectivityList, head_model.Points(:,1), ...
            %         head_model.Points(:,2), head_model.Points(:,3), ...
            %         'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none', 'Parent', ax,...
            %         'DisplayName', 'Scan');
            % hold(ax, 'on');
            % 
            % % Plot electrodes
            % if ~isempty(pcElectrodes_refined)
            %     electrodesPos = pcElectrodes_refined.Location;
            % elseif ~isempty(pcElectrodes_aligned)
            %     electrodesPos = pcElectrodes_aligned.Location;
            % else
            %     electrodesPos = pcElectrodes.Location;
            % end
            % 
            % % Store electrode scatter object for interaction
            % electrodeScatter = scatter3(electrodesPos(:,1), electrodesPos(:,2), electrodesPos(:,3), ...
            %                             40, 'r', 'filled', 'MarkerEdgeColor', 'k', ...
            %                     'Parent', ax, 'DisplayName', 'Electrodes', ...
            %                     'HitTest', 'off', 'PickableParts', 'none');
            % 
            % electrodesPos
            % displayLabels(ax, pointCloud(electrodesPos));
            % 
            % 
            % % Improve Visualization
            % camlight(ax,'headlight');
            % lighting(ax, 'gouraud');
            % grid(ax, 'on');
            % hold(ax, 'off');
            % axis(ax, 'equal')
            % 
            % 
            % % Add labels if necessary
            % xlabel(ax, 'X');
            % ylabel(ax, 'Y');
            % zlabel(ax, 'Z');
        end

        function visualizeElectrodesForRevisionScan()
            cla(revisionAxes);
            
            if isempty(head_scan_aligned)
                uialert(GUIfig, 'Head scan is not aligned. Please align the head scan before revising electrodes.', 'Alignment Error');
                return;
            end
            if isempty(head_model)
                uialert(GUIfig, 'Head scan is not aligned. Please align the head scan before revising electrodes.', 'Alignment Error');
                return;
            end
            
            % Plot the head model trisurf
            trisurf(head_scan_aligned.ConnectivityList, head_scan_aligned.Points(:,1), ...
                    head_scan_aligned.Points(:,2), head_scan_aligned.Points(:,3), ...
                    'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none', 'Parent', revisionAxes,...
                    'DisplayName', 'Scan');
            hold(revisionAxes, 'on');
            
            % Plot electrodes
            if ~isempty(pcElectrodes_refined)
                electrodesPos = pcElectrodes_refined.Location;
            elseif ~isempty(pcElectrodes_aligned)
                electrodesPos = pcElectrodes_aligned.Location;
            else
                electrodesPos = pcElectrodes.Location;
            end
            
            % Store electrode scatter object for interaction
            electrodeScatter = scatter3(electrodesPos(:,1), electrodesPos(:,2), electrodesPos(:,3), ...
                                        40, 'r', 'filled', 'MarkerEdgeColor', 'k', ...
                                'Parent', revisionAxes, 'DisplayName', 'Electrodes', ...
                                'HitTest', 'off', 'PickableParts', 'none');


            % Improve Visualization
            camlight(revisionAxes,'headlight');
            lighting(revisionAxes, 'gouraud');
            grid(revisionAxes, 'on');
            hold(revisionAxes, 'off');
            axis(revisionAxes, 'equal')

            
            % Add labels if necessary
            xlabel(revisionAxes, 'X');
            ylabel(revisionAxes, 'Y');
            zlabel(revisionAxes, 'Z');
            
            legend(revisionAxes, 'show');
        end

        function startElectrodeSelection()
            % Function to start electrode selection
            if isempty(head_scan_aligned)
                uialert(GUIfig, 'Head scan is not aligned. Please align the head scan before revising electrodes.', 'Alignment Error');
                return;
            end
            
            disp('Electrode selection mode activated. Click near an electrode to select it.');
            isSelecting = true;
            
            % Change the cursor to indicate selection mode
            
            set(GUIfig, 'Pointer', 'crosshair');
            
            % Set the axes ButtonDownFcn to handle the selection click
            set(revisionAxes, 'ButtonDownFcn', @selectElectrodeCallback);
        end

        function selectElectrodeCallback(src, event)
            % Callback function for selecting an electrode
            if ~isSelecting
                return;
            end
            
            % Get the click point in data coordinates
            clickPoint = revisionAxes.CurrentPoint(1,1:3);
            
            % Get electrode positions
            if ~isempty(pcElectrodes_refined)
                electrodesPos = pcElectrodes_refined.Location;
            elseif ~isempty(pcElectrodes_aligned)
                electrodesPos = pcElectrodes_aligned.Location;
            else
                electrodesPos = pcElectrodes.Location;
            end
            
            % Calculate distances to all electrodes
            distances = sqrt(sum((electrodesPos - clickPoint).^2, 2));
            
            % Find the closest electrode
            [minDist, idx] = min(distances);
            
            % Define a threshold (adjust based on data scale)
            threshold = 10; % Example value; adjust as needed
            
            if minDist <= threshold
                selectedElectrodeIdx = idx;
                disp(['Selected Electrode ID: ', num2str(idx)]);
                
                % Highlight the selected electrode
                highlightSelectedElectrode(idx);
                
                % Prompt user to select a new position
                uialert(GUIfig, 'Now click on the new position for the selected electrode.', 'Move Electrode');
                
                % Change the ButtonDownFcn to handle the new position click
                set(revisionAxes, 'ButtonDownFcn', @moveElectrodeCallback);
            else
                disp('No electrode nearby.');
                uialert(GUIfig, 'No electrode found near the clicked location. Please try again.', 'Selection Error');
            end
            
            % Reset selection mode
            isSelecting = false;
            set(GUIfig, 'Pointer', 'arrow');
        end

        function moveElectrodeCallback(src, event)
            % Callback function for moving the selected electrode
            if isempty(selectedElectrodeIdx)
                return;
            end
            
            % Get the new click point in data coordinates
            newClickPoint = revisionAxes.CurrentPoint(1,1:3);
            
            % Update the electrode position in the data
            if ~isempty(pcElectrodes_refined)
                previousPositions = [previousPositions; pcElectrodes_refined.Location(selectedElectrodeIdx, :)];
                pcElectrodes_refined.Location(selectedElectrodeIdx, :) = newClickPoint;
            elseif ~isempty(pcElectrodes_aligned)
                previousPositions = [previousPositions; pcElectrodes_aligned.Location(selectedElectrodeIdx, :)];
                pcElectrodes_aligned.Location(selectedElectrodeIdx, :) = newClickPoint;
            else
                previousPositions = [previousPositions; pcElectrodes.Location(selectedElectrodeIdx, :)];
                pcElectrodes.Location(selectedElectrodeIdx, :) = newClickPoint;
            end
            
            disp(['Electrode ', num2str(selectedElectrodeIdx), ' moved to new position: (', ...
                  num2str(newClickPoint(1), '%.2f'), ', ', ...
                  num2str(newClickPoint(2), '%.2f'), ', ', ...
                  num2str(newClickPoint(3), '%.2f'), ')']);
            
            % Refresh the plot to show the updated position
            visualizeElectrodesForRevision();
            
            % Reset the ButtonDownFcn
            set(revisionAxes, 'ButtonDownFcn', '');
        end
    
        function highlightSelectedElectrode(idx)
            % Function to highlight the selected electrode

            % Re-plot electrodes with the selected one highlighted
            cla(revisionAxes);
            
            % Plot the head model trisurf
            trisurf(head_scan_aligned.ConnectivityList, head_scan_aligned.Points(:,1), ...
                    head_scan_aligned.Points(:,2), head_scan_aligned.Points(:,3), ...
                    'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none', 'Parent', revisionAxes, 'DisplayName', 'Scan');
            hold(revisionAxes, 'on');
            
            % Plot electrodes
            if ~isempty(pcElectrodes_refined)
                scatter3(pcElectrodes_refined.Location(:,1), pcElectrodes_refined.Location(:,2), ...
                         pcElectrodes_refined.Location(:,3), 40, 'r', 'filled', 'MarkerEdgeColor', 'k', ...
                         'Parent', revisionAxes, 'DisplayName', 'Electrodes', ...
                         'HitTest', 'off', 'PickableParts', 'none');
            elseif ~isempty(pcElectrodes_aligned)
                scatter3(pcElectrodes_aligned.Location(:,1), pcElectrodes_aligned.Location(:,2), ...
                         pcElectrodes_aligned.Location(:,3), 40, 'r', 'filled', 'MarkerEdgeColor', 'k', ...
                         'Parent', revisionAxes, 'DisplayName', 'Electrodes', ...
                         'HitTest', 'off', 'PickableParts', 'none');
            else
                scatter3(pcElectrodes.Location(:,1), pcElectrodes.Location(:,2), ...
                         pcElectrodes.Location(:,3), 40, 'r', 'filled', 'MarkerEdgeColor', 'k', ...
                         'Parent', revisionAxes, 'DisplayName', 'Electrodes', ...
                         'HitTest', 'off', 'PickableParts', 'none');
            end
            
            % Highlight the selected electrode in yellow
            if ~isempty(pcElectrodes_refined)
                scatter3(pcElectrodes_refined.Location(idx,1), pcElectrodes_refined.Location(idx,2), ...
                         pcElectrodes_refined.Location(idx,3), 100, 'y', 'filled', 'MarkerEdgeColor', 'k', ...
                         'Parent', revisionAxes, 'DisplayName', 'Selected Electrode', ...
                         'HitTest', 'off', 'PickableParts', 'none');
            elseif ~isempty(pcElectrodes_aligned)
                scatter3(pcElectrodes_aligned.Location(idx,1), pcElectrodes_aligned.Location(idx,2), ...
                         pcElectrodes_aligned.Location(idx,3), 100, 'y', 'filled', 'MarkerEdgeColor', 'k', ...
                         'Parent', revisionAxes, 'DisplayName', 'Selected Electrode', ...
                         'HitTest', 'off', 'PickableParts', 'none');
            else
                scatter3(pcElectrodes.Location(idx,1), pcElectrodes.Location(idx,2), ...
                         pcElectrodes.Location(idx,3), 100, 'y', 'filled', 'MarkerEdgeColor', 'k', ...
                         'Parent', revisionAxes, 'DisplayName', 'Selected Electrode', ...
                         'HitTest', 'off', 'PickableParts', 'none');
            end
            
            % Improve Visualization
            camlight('headlight', 'Parent', revisionAxes);
            lighting(revisionAxes, 'gouraud');
            grid(revisionAxes, 'on');
            hold(revisionAxes, 'off');
            
            % Add labels if necessary
            xlabel(revisionAxes, 'X');
            ylabel(revisionAxes, 'Y');
            zlabel(revisionAxes, 'Z');
            
            legend(revisionAxes, 'show');
        end

        function undoLastMove()
            if isempty(previousPositions)
                uialert(GUIfig, 'No moves to undo.', 'Undo Error');
                return;
            end
            
            % Retrieve the last position
            lastPosition = previousPositions(end, :);
            previousPositions(end, :) = []; % Remove the last position from the stack
            
            % Move the electrode back to the last position
            if ~isempty(pcElectrodes_refined)
                pcElectrodes_refined.Location(selectedElectrodeIdx, :) = lastPosition;
            elseif ~isempty(pcElectrodes_aligned)
                pcElectrodes_aligned.Location(selectedElectrodeIdx, :) = lastPosition;
            else
                pcElectrodes.Location(selectedElectrodeIdx, :) = lastPosition;
            end
            
            disp(['Electrode ', num2str(selectedElectrodeIdx), ' moved back to position: (', ...
                  num2str(lastPosition(1), '%.2f'), ', ', ...
                  num2str(lastPosition(2), '%.2f'), ', ', ...
                  num2str(lastPosition(3), '%.2f'), ')']);
            
            % Refresh the plot
            visualizeElectrodesForRevision();
        end

        
        % --- Visualization ---
        function visualizeHeadModel(head_model, ax)
            cla(ax);
    
            trimesh(head_model.ConnectivityList, ...
                head_model.Points(:,1),head_model.Points(:,2),head_model.Points(:,3),'FaceAlpha','0.3','Parent', ax);
            title(ax, '3D Head Model');
            axis(ax, 'equal');
            colormap(ax, 'gray'); 
            view(ax, 3);
        end
        
        function visualizeScan(scan, ax, plot_title)
            %cla(ax);
            points = scan.Location;
            scatter3(points(:,1),points(:,2),points(:,3),'MarkerFaceColor',[0.6350 0.0780 0.1840], 'Parent', ax)
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
    
        function visualizeTriangulation(tri, ax, plot_title)
            %cla(ax);
            
            % Plot triangulation
            hold(ax,'on');
            trisurf(tri.ConnectivityList, tri.Points(:,1), tri.Points(:,2), tri.Points(:,3), ...
                   'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none','Parent',ax);
            title(ax, plot_title)
            
            % Improve Visualization
            grid(ax, 'on');
            axis(ax, 'equal');
            camlight(ax,'headlight');
            lighting(ax,'gouraud');
            set(ax, 'Color', 'white', ...
                'XColor', [0.15 0.15 0.15], ...
                'YColor', [0.15 0.15 0.15], ...
                'ZColor', [0.15 0.15 0.15])
        
            fig = ancestor(ax, 'figure');
            fig.Color = 'white';
            hold(ax,'off');
            view(ax, 0,0)
        end
    
        function visualizePrealigned(scan, model, ax)
            % cla(ax);
            
            scan_points = scan.Location;
            s1 = scatter3(scan_points(:,1),scan_points(:,2),scan_points(:,3), ...
                'MarkerFaceColor', [0.6350 0.0780 0.1840], 'MarkerEdgeColor','none','Parent', ax);
            s1.SizeData = 10;
            hold(ax, 'on');
            
            model = pcdownsample(model, "nonuniformGridSample", 10);
            model_points = model.Location;
            s2 = scatter3(model_points(:,1),model_points(:,2),model_points(:,3), ...
                'MarkerFaceColor','k','MarkerEdgeColor','none','Parent', ax);
            s2.SizeData = 10;
            
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
        
        function visualizeFiducials(X, Y, tform)  
            % Help plot -- only fiducials
            figure(55)
            scatter3(Y(:,1), Y(:,2), Y(:,3), 50, 'r', 'filled', 'MarkerEdgeColor', 'k')
            hold on
            scatter3(X(:,1),X(:,2),X(:,3), 50, 'b', 'filled', 'MarkerEdgeColor', 'k')
            hold off
    
            figure(66)
            scatter3(Y(:,1), Y(:,2), Y(:,3), 50, 'r', 'filled', 'MarkerEdgeColor', 'k')
            hold on
            Xm = pctransform(pointCloud(X), tform);
            Xm = Xm.Location;
            scatter3(Xm(:,1),Xm(:,2),Xm(:,3), 50, 'b', 'filled', 'MarkerEdgeColor', 'k')
            hold off
        end
    
        function visualizeInFiducialPlot(head_model, ax)
            % Plot the STL model using trisurf
            hs = trisurf(head_model.ConnectivityList, head_model.Points(:,1), head_model.Points(:,2), head_model.Points(:,3), ...
                    'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none','Parent',ax);
            cl = camlight(ax, 'right'); 
            lighting(ax, 'flat'); % gouraud
            xlabel(ax, 'X');
            ylabel(ax, 'Y');
            zlabel(ax, 'Z');
            axis(ax, 'tight');
            axis(ax, 'equal');
        
            % Enable interactive rotation
            view(ax,0,0);
            hold(ax, 'on');
        end
    
    
        function visualizeFinalAlignment(aligned_scan, model, ax)
            %cla(ax);

            hold(ax, 'on');
            [colors, minDist, maxDist] = computeDistances(aligned_scan, model);
            %h = pcshow(pointCloud(aligned_scan.Location, 'Color', uint8(colors * 255)), 'Parent', ax, 'MarkerSize', 30);
            scatter3(aligned_scan.Location(:, 1), aligned_scan.Location(:, 2), aligned_scan.Location(:, 3), ...
                        30, colors, 'filled', 'Parent', ax);
            title(ax, 'Final Aligned Model');
            axis(ax, 'equal');

            % Create a secondary (invisible) axis for the colorbar
            if exist('cb_ax', 'var') && isvalid(cb_ax)
                delete(cb_ax); % Remove the secondary axes for the colorbar
            end
            fig = ancestor(ax, 'figure');
            cb = findobj(fig, 'Type', 'Colorbar'); % Find existing colorbars
            if ~isempty(cb)
                delete(cb); % Remove the colorbar
            end
            cb_ax = axes('Position', ax.Position, 'Color', 'none', 'XTick', [], 'YTick', [], 'ZTick', [], ...
                 'Visible', 'off', 'Parent', ax.Parent); 
            cla(cb_ax);
            colormap(cb_ax, 'jet');
            clim(cb_ax, [minDist, maxDist]);
            cb = colorbar(cb_ax, 'Position', [0.9, 0.2, 0.02, 0.6]);
            cb.Label.String = 'Distance';

            view(ax, 3);
            grid(ax, 'on');
    
            set(ax, 'Color', 'white', ...
                    'XColor', [0.15 0.15 0.15], ...
                    'YColor', [0.15 0.15 0.15], ...
                    'ZColor', [0.15 0.15 0.15])
    
            fig = ancestor(ax, 'figure');
            fig.Color = 'white';
    
            hold(ax, 'off');
        end
    
        function visualizeElectrodes(pcElec, head_scan, ax)
            cla(ax);
            title(ax,'Positions of electrodes on head scan.')
            
            trisurf(head_scan.ConnectivityList, head_scan.Points(:,1), head_scan.Points(:,2), head_scan.Points(:,3), ...
                  'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none','Parent',ax);
            hold(ax,'on')
            % Improve Visualization
            grid(ax, 'on');
            axis(ax, 'equal');
            lighting(ax,'gouraud');

            % --- Add multiple light sources ---
            % A light to from camera view:
            camlight(ax, 'headlight');    
            % camlight(ax,180,30);
            
            % A light from the left:
            % camlight(ax, -30, 30);        % 'camlight(az, el)' sets azimuth/elevation relative to camera
            
            % A light from the right:
            % camlight(ax, 30, 30);
            
            % Add more lights at different positions:
            %light(ax, 'Position',[ 1  1  0],'Style','infinite');
            %light(ax, 'Position',[-5  0  0],'Style','infinite');
            %light(ax, 'Position',[ 0  5  0],'Style','infinite');
            %light(ax, 'Position',[ 0 -5  0],'Style','infinite');
            %light(ax, 'Position',[ 0  0  5],'Style','infinite');

            electrode_positions = pcElec.Location;
            scatter3(electrode_positions(:,1),electrode_positions(:,2),electrode_positions(:,3), ...
                40,'r','filled','MarkerEdgeColor','none','Parent',ax)
    
            %displayLabels(ax, pcElec);

            set(ax, 'Color', 'white', ...
                'XColor', [0.15 0.15 0.15], ...
                'YColor', [0.15 0.15 0.15], ...
                'ZColor', [0.15 0.15 0.15])
        
            fig = ancestor(ax, 'figure');
            fig.Color = 'white';
        
            hold(ax,'off');
        end
    
    
        function displayLabels(ax, pcElectrodePositions)
    
            % Display labels
            if isempty(electrodes)
                templateElcPath = fullfile(dataPath, 'template_elc.elc');
                electrodes = elc_read(templateElcPath);
            end
            
            labels = electrodes.labels;
            labelOffset = 1;
            labelPositions = pcElectrodePositions.Location;
            electrodesPositions = pcElectrodePositions.Location;
          
            labelPositions = labelPositions(1:128,:); % last could be nasion
    
            for i = 1:128 % size(electrodesPositions, 1) but last is nasion
                pos = electrodesPositions(i,:);
                offset = labelOffset * sign(pos);
                    
                dr = 3;
            
                if pos(1) > 0
                    offset(1) = offset(1) + dr;
                elseif pos(1) < 0
                    offset(1) = offset(1) - dr;
                end
            
                if pos(2) > 0
                    offset(2) = offset(2) + dr;
                elseif pos(2) < 0
                    offset(2) = offset(2) - dr;
                end
            
                labelPositions(i,:) = pos + offset;
            end
            
            text(labelPositions(:,1), labelPositions(:,2), labelPositions(:,3), ...
                labels, "HorizontalAlignment", "center", "VerticalAlignment", "bottom", ...
                "Color", "k", "FontWeight","bold",'Parent',ax);
    
        end
    
    end
