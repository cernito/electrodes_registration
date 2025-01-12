function [pcAligned_labels] = align_template_electrodes(head_scan, F_head_scan, ax)
    
    % Align template_scan (that includes the electrodes and their positions) 
    % with the head_scan.
    % Possible to locate probable positions of electrodes...
    % If not, could be used as a good initial positions for where
    % to look for the electrodes on the head scan.

    % Input:
    %   head_scan ... already aligned head scan with MRI
    %   F_head_scan ... 3 fiducial lendmark points NAS/LAS/RAS
    % Output:
    %   pcAligned_labels ... exact positions of the aligned electrodes
    %   pcAligned_electrodes ... aligned big electrode templates
    
    % ***************************************************************
    % Load files
    % ***************************************************************
    if isempty(head_scan)
        full_file_path = get_user_file_path('*.stl', 'Select the scan of patients head');
        head_scan = stlread(full_file_path);
        head_points = head_scan.Points;
        pcHead_scan = pointCloud(head_points);
    else
        % If STL:
        head_points = head_scan.Points;
        pcHead_scan = pointCloud(head_points);

        % If PointCloud:
        % head_points = head_scan.Location;
        % pcHead_scan = head_scan;
    end

    % Get the main path (where the script is located)
    mainPath = fileparts(mfilename('fullpath'));
    
    % Define paths for data and results folders
    dataPath = fullfile(mainPath,'..','..','..', 'data');
    resultsPath = fullfile(mainPath,'..','..', '..', 'results');
    
    % Create results directory if it doesn't exist
    if ~exist(resultsPath, 'dir')
        mkdir(resultsPath);
    end
    

    % Load Template Head Scan
    %template = stlread("C:\ČVUT\Bakalarka\P311_template.stl");
    templatePath = fullfile(dataPath, 'P311_template.stl');
    template = stlread(templatePath);
    %template = stlread(get_user_file_path('.stl','Load template scan. Add the filepath to the template instead of this code.'));
    template_points = template.Points;
    pcTemplate = pointCloud(template_points);

    % Load Template Cap
    %template_cap = stlread("C:\ČVUT\Bakalarka\P311_template_cap.stl");
    templateCapPath = fullfile(dataPath, 'P311_template_cap.stl');
    template_cap = stlread(templateCapPath);
    %template_cap = stlread(get_user_file_path('.stl','Load template cap. Add the filepath to the template instead of this code.'));
    template_cap_points = template_cap.Points;
    pcTemplate_cap = pointCloud(template_cap_points);
    
    % Load Electrodes
    %file_path = get_user_file_path('*.elc', 'elc');
    %electrodes = elc_read("C:\ČVUT\Bakalarka\template_elc.elc");
    templateElcPath = fullfile(dataPath, 'template_elc.elc');
    electrodes = elc_read(templateElcPath);
    %electrodes = elc_read(get_user_file_path('.elc','Load tempalte electrodes. Add the filepath to the template instead of this code.'));
    labels = electrodes.labels;
    pcElectrodes = pointCloud(electrodes.pos);
    pcwrite(pcElectrodes, 'electrodes.pcd', 'Encoding', 'ascii');

    
    % ***************************************************************
    % Initial Alignments using Fiducials
    % ***************************************************************
    
    if isempty(F_head_scan)
        % Create the main figure
        mainFig = uifigure('Name', 'Fiducial Selection', ...
                           'Position', [100, 100, 1600, 800]);
    
        % Create two axes within the figure
        ax1 = uiaxes(mainFig, 'Units', 'normalized', 'Position', [0.05 0.1 0.4 0.8]);
        ax1.Title.String = 'Head Scan';
    
        ax2 = uiaxes(mainFig, 'Units', 'normalized', 'Position', [0.55 0.1 0.4 0.8]);
        ax2.Title.String = 'Template Scan';
        
        visualizeInFiducialPlot(template,ax2);
        
        % Call get_fiducials_stl on the first axes
        disp('Select fiducials on the head scan (Figure 1):');
        F_head_scan = get_fiducials_stl(head_scan, ax1);
        disp('Fiducials chosen for head scan:');
        disp(F_head_scan);
        
        visualizeInFiducialPlot(head_scan,ax1);
        
        % Now call get_fiducials_stl on the second axes
        disp('Select fiducials on the head model (Figure 2):');
        F_template = get_fiducials_stl(template, ax2);
        disp('Fiducials chosen for head model:');
        disp(F_template);
    else
        F_template = get_fiducials_stl(template,[]);
        %F_template = [-1.9005, -80.5127, -16.5057;
        %              72.4563,  10.6135, -51.7010;
        %              -71.3962,  8.4479, -51.7655];
    end
    
    pcF_template = pointCloud(F_template);
    pcF_head_scan = pointCloud(F_head_scan);

    % Perform initial ICP alignment using fiducials
    [tform,~] = pcregistericp(pcF_template, pcF_head_scan);

    % =======
    % Apply the initial transformation
    % =======
    pcFidu_Aligned_electrodes = pctransform(pcElectrodes, tform);
    pcFidu_Aligned_template = pctransform(pcTemplate, tform);
    pcFidu_Aligned_template_cap = pctransform(pcTemplate_cap, tform);
    pcF_Fidu_Aligned_template = pctransform(pcF_template, tform);

    
    % figure(1); pcshowpair(pcFidu_Aligned_template, pcHead_scan)

    % ***************************************************************
    % Plot Initial Alignment
    % ***************************************************************

    figure(55); cla;
    title('Alignment of template to head scan.')
    pcshow(pcHead_scan.Location,[0.6350 0.0780 0.1840])
    hold on
    ax1 = pcshow(pcFidu_Aligned_template.Location,[.7 .7 .7]);

    pcshow(pcFidu_Aligned_electrodes.Location, 'r', 'MarkerSize', 200,'Parent',ax1)
    displayLabels(ax1, pcFidu_Aligned_electrodes, electrodes);
    
    grid(ax1, 'on');
    set(ax1, 'Color', 'white', ...
        'XColor', [0.15 0.15 0.15], ...
        'YColor', [0.15 0.15 0.15], ...
        'ZColor', [0.15 0.15 0.15])

    fig = ancestor(ax1, 'figure');
    fig.Color = 'white';
    hold off


    % ***************************************************************
    % Precise alignment using CPD with Refined Parameters
    % ***************************************************************
    
    % Downsample
    pcFidu_Aligned_template_down = pcdownsample(pcFidu_Aligned_template,'gridAverage',10);
    pcHead_scan_down = pcdownsample(pcHead_scan,'gridAverage',10);

    % Estimating Surface Normals - for 'pointToPlane' metric
    %     pcAligned_template_down.Normal = pcnormals(pcAligned_template_down, 30);
    %     pcHead_scan_down.Normal = pcnormals(pcHead_scan_down, 30);
    %
    % Perform ICP refinement
    %     maxIterations = 30;
    %     tolerance = [1e-7, 1e-7];
    %     inlierRatio = 0.95;
    % 
    %     [tform,~,rmse] = pcregistericp(pcAligned_template_down, pcHead_scan_down, ...
    %                                     "Metric","pointToPlane",...
    %                                     "MaxIterations",maxIterations,...
    %                                     "Tolerance",tolerance,...
    %                                     "InlierRatio",inlierRatio,...
    %                                     "Verbose",true);
    % 
    %     pcAligned_electrodes = pctransform(pcAligned_electrodes, tform);
    %     pcAligned_template = pctransform(pcAligned_template, tform);
    % 
    %     disp(['RMSE: ',num2str(rmse)])

    % Compute bounding box in the z-axis
    headZBounds = [min(pcHead_scan.Location(:, 3)), max(pcHead_scan.Location(:, 3))];
    templateZBounds = [min(pcTemplate.Location(:, 3)), max(pcTemplate.Location(:, 3))];
    
    % Calculate the z-axis size
    headZSize = abs(headZBounds(2) - headZBounds(1));
    templateZSize = abs(templateZBounds(2) - templateZBounds(1));
    
    % Compute the z-axis size ratio
    zSizeRatio = headZSize / templateZSize;
    
    % Set CPD parameters based on z-axis size difference
    if zSizeRatio > 1.5
        % Head scan is much taller than the template
        opt.beta = 25;    % Smoother deformation
        opt.lambda = 10;  % More flexible
    elseif zSizeRatio > 1.2
        % Moderate height difference
        opt.beta = 20;
        opt.lambda = 15;
    else
        % Head scan height is close to or smaller than the template
        opt.beta = 15;    % Finer adjustments
        opt.lambda = 19;  % More rigid
    end


    % Perform CPD refinement
    X = [pcHead_scan_down.Location; pcF_head_scan.Location];
    Y = [pcFidu_Aligned_template_down.Location; pcF_Fidu_Aligned_template.Location];

    for b=18
        for l=20.5
            % Init full set of options
            opt.method='nonrigid'; % use nonrigid registration
            
            %opt.beta=b;            % the width of Gaussian kernel (smoothness)
            %opt.lambda=l;          % regularization weight
            
            opt.viz=1;              % show every iteration
            opt.outliers=0.9;       % noise weight
            opt.fgt=0;              % do not use FGT (default)
            opt.normalize=1;        % normalize to unit variance and zero mean before registering (default)
            opt.corresp=0;          % compute correspondence vector at the end of registration (not being estimated by default)
            
            opt.max_it=120;         % max number of iterations
            opt.tol=1e-5;          % tolerance
        
            [Transform, C]=cpd_register(X,Y, opt);
        
            figure,cpd_plot_iter(X, Y); title(['Before | b=',num2str(b),' l=',num2str(l)]);
            figure,cpd_plot_iter(X, Transform.Y);  title(['After registering Y to X | b=',num2str(b),' l=',num2str(l)]);
        end 
    end


    % ======
    % Apply CPD transformation on objects
    % ======
    disp('Applying transformation to whole template and electrodes...')
    transformed_labels = cpd_transform(pcFidu_Aligned_electrodes.Location,Transform);
    transformed_F_template = cpd_transform(pcF_Fidu_Aligned_template.Location,Transform);

    % figure,cpd_plot_iter(X,transformed_electrodes); title('Registered electrodes')
    pcAligned_labels = pointCloud(transformed_labels);
    pcF_Aligned_template = pointCloud(transformed_F_template);
    disp('Applied transformation.')


    % ***************************************************************
    % Plot Refined Alignment
    % ***************************************************************

    % Plot results
    figure(66);cla;
    title('Positions of electrodes on head scan.')
    trisurf(head_scan.ConnectivityList, head_scan.Points(:,1), head_scan.Points(:,2), head_scan.Points(:,3), ...
          'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none');
    hold on
    ax2 = gca;
    s = scatter3(pcAligned_labels.Location(:,1),pcAligned_labels.Location(:,2),pcAligned_labels.Location(:,3), 'r.','Parent',ax2);
    s.SizeData = 600;
    displayLabels(ax2,pcAligned_labels, electrodes);

    % Improve Visualization
    grid(ax2, 'on');
    axis(ax2, 'equal');
    camlight(ax2,'headlight');
    lighting(ax2,'gouraud');
    set(ax2, 'Color', 'white', ...
        'XColor', [0.15 0.15 0.15], ...
        'YColor', [0.15 0.15 0.15], ...
        'ZColor', [0.15 0.15 0.15])

    fig = ancestor(ax2, 'figure');
    fig.Color = 'white';

    hold off;

    % ========
    % Refine alignment
    %
    % - Filter the HD-EEG cap from the scan using aligned template electrode positions.
    % - Creates an alpha shape around electrodes and filters points inside
    %   the shape.
    % ========

    % Segment the HD-EEG cap from scan
    nasion_scan = pcF_Aligned_template.Location(1,:);
    aligned_labels_with_nasion = pointCloud([pcAligned_labels.Location; nasion_scan]);
    segmented_cap = segment_head_cap(head_scan, aligned_labels_with_nasion);

    % Aligned template cap
    template_cap = pcFidu_Aligned_template_cap.Location; % pcFidu_Aligned_template_cap

    figure(99); clf;
    pcshow(segmented_cap.Points,'r');
    hold on
    pcshow(template_cap,'g');

    % Downsample
    dwnScan = pcdownsample(pointCloud(segmented_cap.Points),"gridAverage",8)
    dwnTemplate = pcdownsample(pointCloud(template_cap),"gridAverage",8)

    % Extract Points
    X = [dwnScan.Location; pcF_head_scan.Location];
    Y = [dwnTemplate.Location; pcF_Aligned_template.Location]; % template cap

    % Init full set of options %%%%%%%%%%
    opt.method='nonrigid';  % use nonrigid registration

    %opt.beta=15;            % the width of Gaussian kernel (smoothness)
    %opt.lambda=19.5;           % regularization weight

    opt.viz=1;              % show every iteration
    opt.outliers=0.3;       % noise weight
    opt.fgt=0;              % do not use FGT (default)
    opt.normalize=1;        % normalize to unit variance and zero mean before registering (default)
    opt.corresp=0;          % compute correspondence vector at the end of registration (not being estimated by default)

    opt.max_it=100;         % max number of iterations
    opt.tol=1e-5;           % tolerance

    [Transform, C]=cpd_register(X,Y, opt);

    figure,cpd_plot_iter(X, Y); title('Before');
    figure,cpd_plot_iter(X, Transform.Y);  title('After registering Y to X');

    % =======
    % Apply CPD transformation to objects
    % =======
    transformed_labels = cpd_transform(pcFidu_Aligned_electrodes.Location,Transform);
    transformed_F_template = cpd_transform(pcF_Aligned_template.Location, Transform);

    % figure,cpd_plot_iter(X,transformed_electrodes); title('Registered electrodes')

    pcAligned_labels = pointCloud([transformed_labels; transformed_F_template(1,:)]); % electrodes + nasion
    %pcF_Aligned_template = pointCloud(transformed_F_template);

    if isvalid(mainFig) && ishghandle(mainFig)
        close(mainFig);
    end


    % =========================
    % NESTED FUNCTIONS
    % =========================
    function visualizeInFiducialPlot(head_model, ax)
        % Plot the STL model using trisurf
        hs = trisurf(head_model.ConnectivityList, head_model.Points(:,1), head_model.Points(:,2), head_model.Points(:,3), ...
                'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none','Parent',ax);
        cl = camlight(ax, 'right'); 
        lighting(ax, 'flat'); % gouraud
        xlabel(ax, 'X');
        ylabel(ax, 'Y');
        zlabel(ax, 'Z');
        rotate3d(ax,'off');
        axis(ax, 'tight');
        axis(ax, 'equal');
    
        % Enable interactive rotation
        view(ax,0,0);
        hold(ax, 'on');
    end
    

    function transformedMesh = applyTransformationToMesh(mesh, tform, Transform, useCPD)
        % Extract points and connectivity
        points = mesh.Points;
        connectivityList = mesh.ConnectivityList;
    
        % Apply transformation
        if useCPD
            % Use CPD non-rigid transformation
            transformed_points = cpd_transform(points, Transform);
        else
            % Use rigid transformation
            transformed_points = pctransform(pointCloud(points), tform).Location;
        end
    
        % Reconstruct triangulation with transformed points
        transformedMesh = triangulation(connectivityList, transformed_points);
    end

    
    function displayLabels(ax, pcElectrodes, electrodes)
        % Display electrode labels
        labels = electrodes.labels;
        labelOffset = 1;
        labelPositions = pcElectrodes.Location(1:128,:); % last is nasion
        
        electrodesPositions = pcElectrodes.Location;
        for i = 1:128 %(size(electrodesPositions, 1) - 1) % last is nasion
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
            labels, "HorizontalAlignment", "center", "VerticalAlignment", "bottom", "Color", "red", ...
            'Parent',ax);
    end
end