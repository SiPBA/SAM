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
function sam(appData)

    % For non-GUI executions, the input argument should be:
    %   appdata.class0.files:       Files with images of class 0 (cell array)
    %   appdata.class0.path:        Path of files with images of class 0
    %   appdata.class0.files:       Files with images of class 1 (cell array)
    %   appdata.class0.path:        Path of files with images of class 1
    %   appdata.atlas.file:         Name of the atlas file (NII format)
    %   appdata.atlas.path:         Path of the atlas file
    %   appdata.outpus.file:        Name of the file that will be generated
    %                               (SAM map)
    %   appdata.output.path:        Path of the file that will be generated
    %   appdata.method.fs:          Feature selection method: Posible values: 
    %                               'anova', 'ttest','entropy','wilcoxon','roc'
    %   appdata.method.fsMaxReg:    Number of features seleted
    %   appdata.method.fe:          Feature extraction method: Posible
    %                               values: 'pls
    %   appdata.method.feComp:      Number of features extracted
    %
    
    if nargin<1
        sam_gui()
        return
    end

    if appData.gui, bar = waitbar(0, 'Staring SAM...'); else bar = 0; end
    
    %% Load images and atlas
    appData = loadData(appData, bar);

    %% Estimate accuracy
    
    predictedLabels = [];
    acc = [];
    for reg = 1:appData.atlas.numReg
        disp(['Parcellation of region: ' appData.atlas.nameReg{reg}])
        if bar~=0, waitbar(reg/appData.atlas.numReg, bar, 'Parcellating regions'); end

        voxelReg = appData.atlas.nii.img == reg;    % ROI selection mask
        trnData = appData.images.stack(:,voxelReg); % ROI selection
        trnLabels = appData.images.labels;

        % Feature selection
        disp('Feature selection.....');
        featIdx = sam_featureSelection(trnData, trnLabels, ...
                                appData.method.fs, appData.method.fsMaxReg);
        trnData = trnData(:, featIdx);

        % Feature extraction
        disp('Feature _Extraction.....');
        feats = sam_featureExtraction(trnData, trnLabels, ...
                                appData.method.fe, appData.method.feComp);
        
        % Binary fitting on Feature Space
        disp('Fitting linear classifer (SVM).....');
        t = templateSVM('KernelFunction','linear','Standardize',1);
        mdl = fitcecoc(feats, trnLabels, 'Learners', t, ...
                        'ClassNames', unique(trnLabels));

        % Empirical error, in-sample estimate
        [oofLabels, ~] = predict(mdl, feats);

        predictedLabels(:, reg) = oofLabels;
        acc(reg) = sum(oofLabels==trnLabels)/numel(trnLabels);
    end

    disp(['Averaged accuracy for resubs: ' num2str(mean(acc))]);
    


    %% Statistical Inference

    if bar~=0, waitbar(1, bar, 'Performing statistical inference'); end

	%Comute the theoretical bounds
    alpha = 0.05;
    [boundVC, boundG, boundGZ] = sam_bounds(appData.images.num, ...
                appData.method.feComp, alpha);

    % Compute z statistic
    boundAcc = acc' - boundG;
    [~,p] = sam_zStatistic(boundAcc, mean([boundAcc; acc']), alpha);

    % Compute indices of significative regions
    sigReg = p < alpha;
    
    % Plot bounds
    drawBounds(acc, boundVC, boundG, boundGZ, sigReg, appData.atlas.nameReg, p)
    

    %% Obtaining SAM map

    if bar~=0, waitbar(1, bar, 'Generating SAM map'); end

    % Initialize SAM map
    samMap = zeros(size(appData.atlas.nii.img));

    % Estimate all ROIs. Non-releavant ROIs are set to 0.1 to preserve atlas shape
    for reg = 1:appData.atlas.numReg
        if sigReg(reg), activ = acc(reg); else, activ = 0.1; end
        samMap(appData.atlas.nii.img == reg) = activ;
    end

    % Save SAM map
    origin = appData.atlas.nii.hdr.hist.originator(1:3); 
    pixdim = appData.atlas.nii.hdr.dime.pixdim(2:4);
    datatype = 16;
    samNii = make_nii(samMap, pixdim, origin, datatype);
    save_nii(samNii, fullfile(appData.output.path, appData.output.file));
    
    if bar~=0, close(bar); end
end

function drawBounds(acc, boundVC, boundG, boundGZ, sigReg, nameReg, p)
    figure;
    subplot(1, 2, 1);
    plot(acc);
    hold on
    plot(acc - boundVC, '-*');
    plot(acc - boundG, '-^');
    plot(acc - boundGZ, '-+');
    plot(0.5 * ones(size(acc)), '-.');
    legend({'$1-P_n(\hat{g_n})$','$1-P(\hat{g_n})$ worst case VC',...
        '$1-P(\hat{g_n})$ worst case G','$1-P(\hat{g_n})$ worst case GZ',...
        'random'},'Interpreter','latex');
    xticks(find(sigReg))
    xticklabels(nameReg(sigReg));
    xtickangle(45);
    ylabel('$1-P_n(g_n)$','Interpreter','latex')
    xlabel('116 Standardized Regions')
    grid on;
    set(groot,'defaultAxesTickLabelInterpreter','none')
    hold off

    % Ploting Test and upper bounds
    subplot(1, 2, 2);
    hAx = plotyy(1:length(acc), acc - boundG, 1:length(acc), log(p));
    title('Significance test for a proportion')
    xticks(find(sigReg))
    xticklabels(nameReg(sigReg));
    xtickangle(45);
    set(groot,'defaultAxesTickLabelInterpreter','none')
    xlabel('116 Standardized Regions')
    ylabel(hAx(1),'$1-P_n(g_n)$','Interpreter','latex') % left y-axis
    ylabel(hAx(2),'log(p-value)') % right y-axis
    yticks(hAx(2),[log(0.0005):1:log(0.05) log(0.05):1:log(1)]);
    yticks(hAx(1),0.3:0.1:0.9);
    legend({'$1-P(\hat{g_n})$ worst case G','log(p-value)'},...
            'Interpreter','latex');

end

function appData = loadData(appData, bar)
    stack = [];
    num0 = numel(appData.class0.files);
    num1 = numel(appData.class1.files);

    % Load class 0 images
    disp('Loading images for class 0');
    for j=1:num0
        disp(['Subject ' num2str(j) ' out of ' num2str(num0)]);
        if bar~=0, waitbar(j/num0, bar, 'Loading subjects in class 0'); end
        fileName = fullfile(appData.class0.path, appData.class0.files{j});
        nii = load_nii(fileName);
        if isempty(stack)
            sz = size(nii.img);
            stack = zeros(num0+num1, sz(1), sz(2), sz(3));
        end
        stack(j,:,:,:) = nii.img;
    end

    % Load class 1 images
    disp('Loading images for class 1');
    for j=1:num1
        disp(['Subject ' num2str(j) ' out of ' num2str(num1)]);
        if bar~=0, waitbar(j/num1, bar, 'Loading subjects in class 1'); end
        fileName = fullfile(appData.class1.path, appData.class1.files{j});
        nii = load_nii(fileName);
        if isempty(stack)
            sz = size(nii.img);
            stack = zeros(num0+num1, sz(1), sz(2), sz(3));
        end
        stack(num0+j,:,:,:) = nii.img;
    end

    % Load atlas
    disp('Loading atlas');
    if bar~=0, waitbar(1, bar, 'Loading atlas'); end
    atlasNii = load_nii(fullfile(appData.atlas.path, appData.atlas.file));
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
    appData.images.num = num0 + num1;
    appData.atlas.nii = atlasNii;
    appData.atlas.nameReg = nameReg';
    appData.atlas.numReg = length(nameReg);
    
end
