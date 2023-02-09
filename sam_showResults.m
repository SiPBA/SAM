%
% Plot interactive SAM map
% This function is called from sam.m
%
function sam_showResults(varargin)
    h = varargin{1};
    if ishandle(h)
        f = eval(['@' varargin{3}]);
        f(varargin{4:end});
    else
        createFigure(varargin{:});
    end
end

function setSlices(resData, num)
    pto = get(resData.axs(num),'CurrentPoint');
    resData.point(mod(num-2*(num~=2),3)+1) = round(pto(1,1));
    resData.point(mod(num-2*(num==2),3)+1) = round(pto(1,2));
    plotMapSlices(resData);
end

function createFigure(resData)

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
    
    ui.btnSave = uicontrol( ...
            'Parent', pnBtn, ...
            'Style', 'pushbutton',...
            'Units', 'normalized', ...
            'Position', [0.01 0.05 0.2 0.9], ...
            'Callback', {'sam_showResults','saveMap', f}, ...
            'String', 'Save map');

    ui.btnBounds = uicontrol( ...
            'Parent', pnBtn, ...
            'Style', 'pushbutton',...
            'Units', 'normalized', ...
            'Position', [0.22 0.05 0.2 0.9], ...
            'Callback', {'sam_showResults','drawBounds', f}, ...
            'String', 'Show bound');

    ui.btnTest = uicontrol( ...
            'Parent', pnBtn, ...
            'Style', 'pushbutton',...
            'Units', 'normalized', ...
            'Position', [0.43 0.05 0.2 0.9], ...
            'Callback', {'sam_showResults','drawTest', f}, ...
            'String', 'Show test');
        
    ui.lblBound = uicontrol( ...
            'Parent', pnBtn, ...
            'Style', 'text',...
            'HorizontalAlignment','left', ...
            'Units', 'normalized', ...
            'Position', [0.64 0.05 0.15 0.9], ...
            'String', 'Bound');

    ui.pumBound = uicontrol( ...
            'Parent', pnBtn, ...
            'Style', 'popupmenu',...
            'String', resData.bound.list, ...
            'Units', 'normalized', ...
            'Position', [0.70 0.05 0.2 0.9], ...
            'Callback', {'sam_showResults','setBound', f}, ...
            'Value',  find(strcmp(resData.bound.list, resData.bound.name)));

    % Store data
    resData.ui = ui;
    resData.axs = [a1 a2 a3 a4];
    guidata(f, resData);

    % Create map and plot slices
    createMapVolume(f);
end
   
function setBound(f)
    resData = guidata(f);

    % Get name of selected bound
    bounds = get(resData.ui.pumBound, 'String');
    resData.bound.name = bounds{get(resData.ui.pumBound, 'Value')};
    resData.bound.value = sam_bound(resData.bound.name, resData.n, ...
                                    resData.dim, resData.alpha, ...
                                    resData.model, resData.acc, ...
                                    resData.bound.dropout);

    % Compute SAM map
    [resData.map, resData.p, resData.sigReg] = sam_map(resData.acc, ...
                        resData.bound.value, resData.alpha, resData.atlas);

    guidata(f, resData);

    % Create map and plot slices
    createMapVolume(f);
end

function createMapVolume(f)
    resData = guidata(f);

    % Binarize map
    binMap = resData.map;
    binMap(binMap >  0.1) = 1;
    binMap(isnan(binMap)) = 0;
    
    % Set background between 0 and 0.5
    bg = resData.mapBg;
    bg = bg - min(bg(:));
    bg = bg ./ max(bg(:));
    bg = bg .* 0.5;
    
    resData.binMap = binMap;
    resData.bg = bg;
    resData.point = uint16(size(bg)/2);
    guidata(f, resData);

    plotMapSlices(resData)
end

function plotMapSlices(resData)

    x = resData.point(1); 
    y = resData.point(2); 
    z = resData.point(3);

    % Get slices from background
    s1Bg = squeeze(resData.bg(x,:,:));
    s2Bg = squeeze(resData.bg(:,y,:));
    s3Bg = squeeze(resData.bg(:,:,z));

    % Get slices from map
    s1Map = squeeze(resData.binMap(x,:,:));
    s2Map = squeeze(resData.binMap(:,y,:));
    s3Map = squeeze(resData.binMap(:,:,z));
    
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
    im1 = image(resData.axs(1), s1c);
    im2 = image(resData.axs(2), s2c);
    im3 = image(resData.axs(3), s3c);
    set(im1,'ButtonDownFcn',{'sam_showResults','setSlices',resData, 1});
    set(im2,'ButtonDownFcn',{'sam_showResults','setSlices',resData, 2});
    set(im3,'ButtonDownFcn',{'sam_showResults','setSlices',resData, 3});
    set(resData.axs(1), 'DataAspectRatio', [1 1 1]);
    set(resData.axs(2), 'DataAspectRatio', [1 1 1]);
    set(resData.axs(3), 'DataAspectRatio', [1 1 1]);
    axis(resData.axs(1),'off')
    axis(resData.axs(2),'off')
    axis(resData.axs(3),'off')

    % Draw crosshair
    hold(resData.axs(1),'on');
    plot(resData.axs(1),[0 size(s1c,2)],[y y],'w');
    plot(resData.axs(1),[z z],[0 size(s1c,1)],'w');
    hold(resData.axs(1),'off');
    
    hold(resData.axs(2),'on');
    plot(resData.axs(2),[0 size(s2c,2)],[x x],'w');
    plot(resData.axs(2),[z z],[0 size(s2c,1)],'w');
    hold(resData.axs(2),'off');

    hold(resData.axs(3),'on');
    plot(resData.axs(3),[0 size(s3c,2)],[x x],'w');
    plot(resData.axs(3),[y y],[0 size(s3c,1)],'w');
    hold(resData.axs(3),'off');

    % Draw info
    regIdx = resData.atlas.nii.img(x,y,z);
    if regIdx > 0
        regName = strrep(resData.atlas.nameReg{regIdx}, '_', ' ');
        pValue = num2str(resData.p(regIdx));
        pValueLog = num2str(log(resData.p(regIdx)));
    else
        regName = '/';
        pValue = '/';
        pValueLog = '/';
    end

    acc = resData.map(x,y,z);
    if isnan(acc) || acc==0
        accCor = '/';
        acc = '/'; 
    else 
        boundReg = resData.bound.value;
        if numel(boundReg) > 1, boundReg = boundReg(regIdx); end
        accCor = num2str(round((acc - boundReg)*100, 2));
        acc = num2str(round(acc*100, 2)); 
    end

    pos = ['(' num2str(x) ', ' num2str(y) ', ' num2str(z) ')'];
    cla(resData.axs(4));
    text(resData.axs(4), 0, 0.90, ['Position = ' pos]);
    text(resData.axs(4), 0, 0.75, ['Region = ' regName]);
    text(resData.axs(4), 0, 0.60, ['Accuracy = ' acc]);
    text(resData.axs(4), 0, 0.45, ['Accuracy (corrected) = ' accCor]);
    text(resData.axs(4), 0, 0.30, ['p-Value = ' pValue]);
    text(resData.axs(4), 0, 0.15, ['log(p-Value) = ' pValueLog]);
end

function drawBounds(f)
    resData = guidata(f);

    acc = resData.acc;
    bound = resData.bound.value;
    sigReg = resData.sigReg;
    nameReg = resData.atlas.nameReg;
    
    figure;
    plot(acc, 'DisplayName', 'Empirical');
    hold on
    plot(acc - bound , 'DisplayName', 'Corrected');
    hold off
    legend('show','Interpreter','latex')
    xticks(find(sigReg))
    xticklabels(nameReg(sigReg));
    xtickangle(45);
    ylabel('$1-P_n(g_n)$','Interpreter','latex')
    xlabel('116 Standardized Regions')
    grid on;
    set(groot,'defaultAxesTickLabelInterpreter','none')
end

function drawTest(f)
    resData = guidata(f);

    acc = resData.acc;
    bound = resData.bound.value;
    sigReg = resData.sigReg;
    nameReg = resData.atlas.nameReg;
    p = resData.p;
    
    figure;
    if exist('yyaxis', 'file')
        yyaxis left
        plot(1:length(acc), acc - bound);
        yticks(0.3:0.1:0.9);
        ylabel('$1-P_n(g_n)$','Interpreter','latex')

        yyaxis right
        plot(1:length(acc), log(p));
        yticks([log(0.0005):1:log(0.05) log(0.05):1:log(1)]);
        ylabel('log(p-value)') % right y-axis
    else
        hAx = plotyy(1:length(acc), acc - bound, 1:length(acc), log(p));
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
    legend({['$1-P(\hat{g_n})$ worst case ' resData.bound.name],'log(p-value)'},...
            'Interpreter','latex');
end

function saveMap(f)
    resData = guidata(f);
    [file, path] = uiputfile(['./*.nii'],'Select map file');
    if ~isequal(file, 0)
        outFile = fullfile(path, file);
        origin = resData.atlas.nii.hdr.hist.originator(1:3); 
        pixdim = resData.atlas.nii.hdr.dime.pixdim(2:4);
        datatype = 16;
        samNii = make_nii(resData.map, pixdim, origin, datatype);
        save_nii(samNii, outFile);
        disp(['SAM map saved: ' outFile]);
    end
end

