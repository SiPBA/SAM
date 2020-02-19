% Plot interactive SAM map
%
% Usage:
%   plotMap(bg, map, pMap)
% Input:
%   bg:       3-D matrix used as background.
%   map:      SAM map. It shold be a 3-D matrix of same size as bg.
%   pMap:     Map with p-Values. It should have same size as bg and map.
%
function sam_showResults(varargin)
    h = varargin{1};
    if ishandle(h)
        f = eval(['@' varargin{3}]);
        f(varargin{4:end});
    else
        showResults(varargin{:});
    end
end

function setSlices(plotData, num)
    pto = get(plotData.axs(num),'CurrentPoint');
    plotData.point(mod(num-2*(num~=2),3)+1) = round(pto(1,1));
    plotData.point(mod(num-2*(num==2),3)+1) = round(pto(1,2));
    plotMapSlices(plotData);
end

function showResults(res, atlas)

    % Binarize map
    binMap = res.map;
    binMap(binMap >  0.1) = 1;
    binMap(isnan(binMap)) = 0;
    
    % Set background between 0 and 0.5
    bg = res.mapBg;
    bg = bg - min(bg(:));
    bg = bg ./ max(bg(:));
    bg = bg .* 0.5;
    
    % Draw figure
    f = figure('Name', 'SAM Map', 'NumberTitle', 'off', ...
            'MenuBar','none','ToolBar','none','DockControls','off');
        
    pnBtn = uipanel(f, 'Units','normalized','Position',[0 0.92 1 0.08]);
    pnAxs = uipanel(f, 'Units','normalized','Position',[0 0 1 0.92], ...
                    'BorderType', 'none');
    a1 = axes(pnAxs, 'Units','normalized','Position',[0.05 0.53 0.42 0.42]);
    a2 = axes(pnAxs, 'Units','normalized','Position',[0.53 0.53 0.42 0.42]);
    a3 = axes(pnAxs, 'Units','normalized','Position',[0.05 0.05 0.42 0.42]);
    a4 = axes(pnAxs, 'Units','normalized','Position',[0.53 0.05 0.42 0.42]);
    axis(a4, 'off');
    
    btnSave = uicontrol( ...
            'Parent', pnBtn, ...
            'Style', 'pushbutton',...
            'Units', 'normalized', ...
            'Position', [0.01 0.05 0.2 0.9], ...
            'Callback', {'sam_showResults','saveMap', res, atlas}, ...
            'String', 'Save map');

    btnBounds = uicontrol( ...
            'Parent', pnBtn, ...
            'Style', 'pushbutton',...
            'Units', 'normalized', ...
            'Position', [0.22 0.05 0.2 0.9], ...
            'Callback', {'sam_showResults','drawBounds', res, atlas}, ...
            'String', 'Show bounds');

    btnTest = uicontrol( ...
            'Parent', pnBtn, ...
            'Style', 'pushbutton',...
            'Units', 'normalized', ...
            'Position', [0.43 0.05 0.2 0.9], ...
            'Callback', {'sam_showResults','drawTest', res, atlas}, ...
            'String', 'Show test');
        
    % Store data
    plotData.map = res.map;
    plotData.binMap = binMap;
    plotData.bg = bg;
    plotData.axs = [a1 a2 a3 a4];
    plotData.point = uint16(size(bg)/2);
    plotData.boundG = res.boundG;
    plotData.p = res.p;
    plotData.atlas = atlas;

    % Plot slices
    plotMapSlices(plotData);
end
    
function plotMapSlices(plotData)
    x = plotData.point(1); 
    y = plotData.point(2); 
    z = plotData.point(3);

    % Get slices from background
    s1Bg = squeeze(plotData.bg(x,:,:));
    s2Bg = squeeze(plotData.bg(:,y,:));
    s3Bg = squeeze(plotData.bg(:,:,z));

    % Get slices from map
    s1Map = squeeze(plotData.binMap(x,:,:));
    s2Map = squeeze(plotData.binMap(:,y,:));
    s3Map = squeeze(plotData.binMap(:,:,z));
    
    % Set first slice for color image
    s1c(:,:,1) = s1Bg + s1Map;
    s1c(:,:,2) = s1Bg + s1Map;
    s1c(:,:,3) = s1Bg;
    s1c = min(s1c, 1);

    % Set second slice for color image
    s2c(:,:,1) = s2Bg + s2Map;
    s2c(:,:,2) = s2Bg + s2Map;
    s2c(:,:,3) = s2Bg;
    s2c = min(s2c, 1);

    % Set third slice for color image
    s3c(:,:,1) = s3Bg + s3Map;
    s3c(:,:,2) = s3Bg + s3Map;
    s3c(:,:,3) = s3Bg;
    s3c = min(s3c, 1);
    
    % Draw slices
    im1 = image(plotData.axs(1), s1c);
    im2 = image(plotData.axs(2), s2c);
    im3 = image(plotData.axs(3), s3c);
    set(im1,'ButtonDownFcn',{'sam_showResults','setSlices',plotData, 1});
    set(im2,'ButtonDownFcn',{'sam_showResults','setSlices',plotData, 2});
    set(im3,'ButtonDownFcn',{'sam_showResults','setSlices',plotData, 3});
    set(plotData.axs(1), 'DataAspectRatio', [1 1 1]);
    set(plotData.axs(2), 'DataAspectRatio', [1 1 1]);
    set(plotData.axs(3), 'DataAspectRatio', [1 1 1]);
    axis(plotData.axs(1),'off')
    axis(plotData.axs(2),'off')
    axis(plotData.axs(3),'off')

    % Draw crosshair
    hold(plotData.axs(1),'on');
    plot(plotData.axs(1),[0 size(s1c,2)],[y y],'w');
    plot(plotData.axs(1),[z z],[0 size(s1c,1)],'w');
    hold(plotData.axs(1),'off');
    
    hold(plotData.axs(2),'on');
    plot(plotData.axs(2),[0 size(s2c,2)],[x x],'w');
    plot(plotData.axs(2),[z z],[0 size(s2c,1)],'w');
    hold(plotData.axs(2),'off');

    hold(plotData.axs(3),'on');
    plot(plotData.axs(3),[0 size(s3c,2)],[x x],'w');
    plot(plotData.axs(3),[y y],[0 size(s3c,1)],'w');
    hold(plotData.axs(3),'off');

    % Draw info
    acc = plotData.map(x,y,z);
    if isnan(acc) || acc==0
        accCor = '/';
        acc = '/'; 
    else 
        accCor = num2str(round((acc - plotData.boundG)*100, 2));
        acc = num2str(round(acc*100, 2)); 
    end
    
    regIdx = plotData.atlas.nii.img(x,y,z);
    if regIdx > 0
        regName = strrep(plotData.atlas.nameReg{regIdx}, '_', ' ');
        pValue = num2str(plotData.p(regIdx));
        pValueLog = num2str(log(plotData.p(regIdx)));
    else
        regName = '/';
        pValue = '/';
        pValueLog = '/';
    end
    pos = ['(' num2str(x) ', ' num2str(y) ', ' num2str(z) ')'];
    cla(plotData.axs(4));
    text(plotData.axs(4), 0, 0.90, ['Position = ' pos]);
    text(plotData.axs(4), 0, 0.75, ['Region = ' regName]);
    text(plotData.axs(4), 0, 0.60, ['Accuracy = ' acc]);
    text(plotData.axs(4), 0, 0.45, ['Accuracy (corrected) = ' accCor]);
    text(plotData.axs(4), 0, 0.30, ['p-Value = ' pValue]);
    text(plotData.axs(4), 0, 0.15, ['log(p-Value) = ' pValueLog]);
end

function drawBounds(res, atlas)
    acc = res.acc;
    boundVC = res.boundVC;
    boundG = res.boundG;
    boundGZ = res.boundGZ;
    sigReg = res.sigReg;
    nameReg = atlas.nameReg;
    p = res.p;
    
    figure;
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
end

function drawTest(res, atlas)
    acc = res.acc;
    boundG = res.boundG;
    sigReg = res.sigReg;
    nameReg = atlas.nameReg;
    p = res.p;
    
    figure;
    if exist('yyaxis', 'file')
        yyaxis left
        plot(1:length(acc), acc - boundG);
        yticks(0.3:0.1:0.9);
        ylabel('$1-P_n(g_n)$','Interpreter','latex')

        yyaxis right
        plot(1:length(acc), log(p));
        yticks([log(0.0005):1:log(0.05) log(0.05):1:log(1)]);
        ylabel('log(p-value)') % right y-axis
    else
        hAx = plotyy(1:length(acc), acc - boundG, 1:length(acc), log(p));
        yticks(hAx(1),0.3:0.1:0.9);
        ylabel(hAx(1),'$1-P_n(g_n)$','Interpreter','latex') % left y-axis
        yticks(hAx(2),[log(0.0005):1:log(0.05) log(0.05):1:log(1)]);
        ylabel(hAx(2),'log(p-value)') % right y-axis
    end
    
    title('Significance test for a proportion')
    xticks(find(sigReg))
    xticklabels(nameReg(sigReg));
    xtickangle(45);
    set(groot,'defaultAxesTickLabelInterpreter','none')
    xlabel('116 Standardized Regions')
    legend({'$1-P(\hat{g_n})$ worst case G','log(p-value)'},...
            'Interpreter','latex');
end

function saveMap(res, atlas)
    [file, path] = uiputfile(['./*.nii'],'Select map file');
    if ~isequal(file, 0)
        outFile = fullfile(path, file);
        origin = atlas.nii.hdr.hist.originator(1:3); 
        pixdim = atlas.nii.hdr.dime.pixdim(2:4);
        datatype = 16;
        samNii = make_nii(res.map, pixdim, origin, datatype);
        save_nii(samNii, outFile);
        disp(['SAM map saved: ' outFile]);
    end
end