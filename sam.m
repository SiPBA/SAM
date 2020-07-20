%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Statistical Agnostic Mappping (SAM) %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Please cite:
%   J.M. Gorriz et al. for the Alzheimers Disease Neuroimaging Initiative 
%   (ADNI) and for the Parkinsons Progression Markers Initiative (PPMI). 
%   STATISTICAL AGNOSTIC MAPPING: A FRAMEWORK IN NEUROIMAGING BASED ON 
%   CONCENTRATION INEQUALITIES.
%   Doi: https://doi.org/10.1101/2019.12.27.889436.
%
function out = sam(appData)

    % For non-GUI executions, the input argument should be:
    %   appdata.class0.files:       Files with images of class 0 (cell array)
    %   appdata.class0.files:       Files with images of class 1 (cell array)
    %   appdata.method.fs:          Feature selection method: Posible values: 
    %                               'anova', 'ttest','entropy','wilcoxon','roc'
    %   appdata.method.fsMaxReg:    Number of features seleted
    %   appdata.method.fe:          Feature extraction method: Posible
    %                               values: 'pls'
    %   appdata.method.feComp:      Number of features extracted
    %
    
    %% Initial checks
    
    % Check LOADNII toolbox
    if exist('load_nii.m','file') ~= 2
        error(['NIFTI toolbox required. Please download and install it from: ' ...
               'es.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image']);
    end

    % If no input arguments, launch GUI
    if nargin<1, sam_gui(); return; end

    if ~isfield(appData, 'gui'), appData.gui = 0; end
    if ~isfield(appData, 'verbose'), appData.verbose = 1; end
    if ~isfield(appData, 'stat'), appData.stat = []; end
    if ~isfield(appData.stat, 'boundMethod'), appData.stat.boundMethod = ''; end
    if ~isfield(appData.stat, 'alpha'), appData.stat.alpha = 0.05; end
    
    % Show waitbar if GUI was loaded
    if appData.gui, bar = waitbar(0, 'Starting SAM...'); else bar = 0; end
        
    %% Save self path
    appData.samPath = fileparts(which('sam'));

    %% Load images and atlas
    if ~isfield(appData, 'images')
        appData = loadData(appData, bar);
        if isempty(appData.images.stack) || isempty(appData.atlas.nii)
            if bar ~= 0, close(bar), end
            return
        end
    end
    %% Estimate accuracy
    
    predictedLabels = [];
    acc = [];
    if appData.verbose, fprintf('Analyzing region:                     \n'); end
    for reg = 1:appData.atlas.numReg
        if appData.verbose
            fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
            fprintf('%-20s\n', appData.atlas.nameReg{reg}(1:min(end,20)))
        end
        if bar~=0, waitbar(reg/appData.atlas.numReg, bar, 'Analyzing regions'); end

        voxelReg = appData.atlas.nii.img == reg;    % ROI selection mask
        trnData = appData.images.stack(:,voxelReg); % ROI selection
        trnLabels = appData.images.labels;
        
        % Feature selection
        featIdx = sam_featureSelection(trnData, trnLabels, ...
                                appData.method.fs, appData.method.fsMaxReg);
        trnData = trnData(:, featIdx);

        if isempty(featIdx)
            acc(reg) = 0;
            continue
        end
        
        % Feature extraction
        feats = sam_featureExtraction(trnData, trnLabels, ...
                                appData.method.fe, appData.method.feComp);
        
        try
            % Binary fitting on feature space
            t = templateSVM('KernelFunction','linear','Standardize',1);
            mdl = fitcecoc(feats, trnLabels, 'Learners', t, ...
                        'ClassNames', unique(trnLabels));

             % Empirical error, in-sample estimate
            [oofLabels, ~] = predict(mdl, feats);
            predictedLabels(:, reg) = oofLabels;
            acc(reg) = sum(oofLabels==trnLabels)/numel(trnLabels);
        catch
            acc(reg) = 0;
        end

    end

    if appData.verbose, 
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bAll\n');
        disp(['Averaged accuracy for resubs: ' num2str(mean(acc))]);
    end


    %% Statistical Inference

    if bar~=0, waitbar(1, bar, 'Performing statistical inference'); end

	%Comute the theoretical bounds
    [boundVC, boundG, boundGZ] = sam_bounds(appData.images.num, ...
                appData.method.feComp, appData.stat.alpha);

    % Compute z statistic
    bound = boundG;
    if strcmpi(appData.stat.boundMethod, 'Vapnik'), bound = boundVC; end
    [~,p] = sam_zStatistic(acc, bound, appData.stat.alpha);

    % Compute indices of significative regions
    sigReg = p < appData.stat.alpha;
    
    %% Obtaining SAM map

    if bar~=0, waitbar(1, bar, 'Generating SAM map'); end

    % Initialize SAM map y p-Value map
    samMap = zeros(size(appData.atlas.nii.img));

    % Estimate all ROIs. Non-releavant ROIs are set to NaN
    for reg = 1:appData.atlas.numReg
        if sigReg(reg), activ = acc(reg); else, activ = nan; end
        samMap(appData.atlas.nii.img == reg) = activ;
    end

    % Close wait bar
    if bar~=0, close(bar); end

    % Show results
    res.acc = acc;
    res.p = p;
    res.boundVC = boundVC;
    res.boundG = boundG;
    res.boundGZ = boundGZ;
    res.sigReg = sigReg;
    res.mapBg = squeeze(mean(appData.images.stack));  % Background
    res.map = samMap;
    %res.pMap = pMap;
    if appData.gui, sam_showResults(res, appData.atlas); end
    if nargout > 0, out = res; end
       
    
end


function appData = loadData(appData, bar)
    num0 = numel(appData.class0.files);
    num1 = numel(appData.class1.files);
    num = num0 + num1;
    files = [appData.class0.files appData.class1.files];

    % Initialize variables
    stack = [];
    sz = [];
    vx = [];
    atlasNii = [];

    % Load images
    fprintf('Loading images        \n')
    for j=1:num
        fprintf('\b\b\b\b\b\b\b\b%3d/%3d\n', j, num)
        if bar~=0, waitbar(j/num, bar, 'Loading images'); end
        nii = load_nii(files{j});
        if isempty(stack)
            sz = size(nii.img);
            vx = nii.hdr.dime.pixdim(2:4);
            stack = zeros(num, sz(1), sz(2), sz(3));
        elseif ~all(sz==size(nii.img)) || ~all(vx==nii.hdr.dime.pixdim(2:4))
            msg = ['ERROR! Images have different sizes. All of them ' ...
                    'should registered using the same template.'];
            msgbox(msg)
            disp(msg)
            stack = [];
            break
        end
        stack(j,:,:,:) = nii.img;
    end

    % Load atlas
    if ~isempty(stack)
        disp('Loading atlas');
        if bar~=0, waitbar(1, bar, 'Loading atlas'); end
        atlasFile = fullfile(appData.samPath, 'atlas', ...
            ['atlas' num2str(sz,'_%i') num2str(round(vx*10),'_%i') '.nii']);
        if exist(atlasFile, 'file')
            atlasNii = load_nii(atlasFile);
        else
            msg = ['ERROR! Invalid image size. Images must be of one ' ...
                   'of the following sizes:\n'];
            ls = dir(fullfile(appData.samPath, 'atlas', 'atlas*.nii'));
            for i=1:numel(ls)
                [~, n, ~] = fileparts(ls(i).name);
                p = strsplit(n, '_');
                if length(p)<7, continue; end
                vxa = str2num([p{5} ' ' p{6} ' ' p{7}])/10;
                msg = [msg '- Volume: ' p{2} ' x ' p{3} ' x ' p{4} ...
                           '  Voxel size: ' num2str(vxa(1)) ' x ' ...
                           num2str(vxa(2)) ' x ' num2str(vxa(3)) '\n'];
            end
            msgbox(sprintf(msg))
            fprintf(msg)
        end
    end
    
    nameReg = {'Precentral_L','Precentral_R','Frontal_Sup_L','Frontal_Sup_R', ...
        'Frontal_Sup_Orb_L','Frontal_Sup_Orb_R','Frontal_Mid_L','Frontal_Mid_R',...
        'Frontal_Mid_Orb_L', 'Frontal_Mid_Orb_R','Frontal_Inf_Oper_L', ...
        'Frontal_Inf_Oper_R','Frontal_Inf_Tri_L','Frontal_Inf_Tri_R', ...
        'Frontal_Inf_Orb_L','Frontal_Inf_Orb_R','Rolandic_Oper_L', ...
        'Rolandic_Oper_R','Supp_Motor_Area_L','Supp_Motor_Area_R','Olfactory_L',...
        'Olfactory_R','Frontal_Sup_Medial_L','Frontal_Sup_Medial_R', ...
        'Frontal_Mid_Orb_L','Frontal_Mid_Orb_R','Rectus_L','Rectus_R','Insula_L',...
        'Insula_R','Cingulum_Ant_L','Cingulum_Ant_R','Cingulum_Mid_L', ...
        'Cingulum_Mid_R','Cingulum_Post_L','Cingulum_Post_R','Hippocampus_L',...
        'Hippocampus_R','ParaHippocampal_L','ParaHippocampal_R','Amygdala_L',...
        'Amygdala_R','Calcarine_L','Calcarine_R','Cuneus_L','Cuneus_R',...
        'Lingual_L','Lingual_R','Occipital_Sup_L','Occipital_Sup_R', ...
        'Occipital_Mid_L','Occipital_Mid_R','Occipital_Inf_L','Occipital_Inf_R',...
        'Fusiform_L','Fusiform_R','Postcentral_L','Postcentral_R','Parietal_Sup_L',...
        'Parietal_Sup_R','Parietal_Inf_L','Parietal_Inf_R','SupraMarginal_L',...
        'SupraMarginal_R','Angular_L','Angular_R','Precuneus_L','Precuneus_R',...
        'Paracentral_Lobule_L','Paracentral_Lobule_R','Caudate_L','Caudate_R',...
        'Putamen_L', 'Putamen_R','Pallidum_L','Pallidum_R','Thalamus_L',...
        'Thalamus_R','Heschl_L','Heschl_R','Temporal_Sup_L','Temporal_Sup_R',...
        'Temporal_Pole_Sup_L','Temporal_Pole_Sup_R','Temporal_Mid_L',...
        'Temporal_Mid_R','Temporal_Pole_Mid_L','Temporal_Pole_Mid_R',...
        'Temporal_Inf_L','Temporal_Inf_R','Cerebelum_Crus1_L','Cerebelum_Crus1_R',...
        'Cerebelum_Crus2_L','Cerebelum_Crus2_R','Cerebelum_3_L','Cerebelum_3_R',...
        'Cerebelum_4_5_L','Cerebelum_4_5_R','Cerebelum_6_L','Cerebelum_6_R',...
        'Cerebelum_7b_L','Cerebelum_7b_R','Cerebelum_8_L','Cerebelum_8_R',...
        'Cerebelum_9_L','Cerebelum_9_R','Cerebelum_10_L','Cerebelum_10_R',...
        'Vermis_1_2','Vermis_3','Vermis_4_5','Vermis_6','Vermis_7','Vermis_8',...
        'Vermis_9','Vermis_10'};
    
    
    % Store loaded data
    appData.images.stack = stack;
    appData.images.labels = [zeros(num0, 1); ones(num1, 1)];
    appData.images.num0 = num0;
    appData.images.num1 = num1;
    appData.images.num = num;
    appData.atlas.nii = atlasNii;
    appData.atlas.nameReg = nameReg';
    appData.atlas.numReg = length(nameReg);
    
end
