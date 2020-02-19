%
% Implement a GUI for SAM toolbox.
% This function is called from sam.m
%
function sam_gui(handle, event, action, params)
    if nargin < 1
        drawGUI()
    else
        f = eval(['@' action]);
        f(handle, params{:});
    end
end

function drawGUI()
    fsMethods = {'anova','ttest','entropy','wilcoxon','roc','bhattacharyya'};
    feMethods = {'pls'};
    
    wWin = 600;
    hWin = 650;
    screenSize = get (0, 'ScreenSize');
    mainWinPos = [round((screenSize(3) - wWin)/2) ...
                  round((screenSize(4) - hWin)/2)]; 
              
    uiObj.mainWin = figure('Name', 'SAM', ...
            'NumberTitle', 'off', ...
            'Menubar', 'none', ...
            'Resize', 'off', ...
            'Position',[mainWinPos wWin hWin], ...
            'DockControls', 'off');
    
    uiObj.panelTitle  = uipanel(uiObj.mainWin, ...
            'Units', 'pixels', ...
            'Position',[5 hWin-5-50 wWin-10 50]);
        
    uiObj.title = uicontrol(uiObj.panelTitle, ...
            'Style','text',...
            'String', 'STATISTICAL AGNOSTIC MAPPING', ...
            'FontSize', 20, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment','center', ...
            'Units','normalized',...
            'Position', [0.05 0 0.9 0.85]);

    %%%%%%%%%%%%%%%%%%%%%
    % Panel for class 0 %
    %%%%%%%%%%%%%%%%%%%%%

    uiObj.panelClass0  = uipanel(uiObj.mainWin, ...
            'Units', 'pixels', ...
            'Title', 'Class 0', ...
            'FontSize', 11, ...
            'FontWeight', 'bold', ...
            'Position',[5 hWin-60-200 wWin-10 200]);
        
    uiObj.class0NameLbl = uicontrol(uiObj.panelClass0, ...
            'Style','text', ...
            'String', 'Name:', ...
            'HorizontalAlignment','left', ...
            'FontSize', 10, ...
            'FontWeight', 'normal', ...
            'Units', 'pixels', ...
            'Position', [10, 150, 60, 20]);

    uiObj.class0Name = uicontrol(uiObj.panelClass0, ...
            'Style','edit', ...
            'String', 'class0', ...
            'HorizontalAlignment','left', ...
            'FontSize', 10, ...
            'FontWeight', 'normal', ...
            'Units', 'pixels', ...
            'Position', [60, 150, 390, 25]);

	uiObj.class0Btn = uicontrol(uiObj.panelClass0, ...
            'Style','pushbutton', ...
            'String', 'Select images', ...
            'HorizontalAlignment','left', ...
            'FontSize', 10, ...
            'FontWeight', 'normal', ...
            'Callback',{'sam_gui','selectFiles',{'0'}}, ...
            'Units', 'pixels', ...
            'Position', [460, 150, 120, 28]);
        
	uiObj.class0Files = uicontrol(uiObj.panelClass0, ...
            'Style','list', ...
            'String', '', ...
            'HorizontalAlignment','left', ...
            'FontSize', 10, ...
            'FontWeight', 'normal', ...
            'Units', 'pixels', ...
            'Position', [10, 35, wWin-30, 105]);
        
    uiObj.class0NumFiles = uicontrol(uiObj.panelClass0, ...
            'Style','text', ...
            'String', '0 images selected', ...
            'HorizontalAlignment','left', ...
            'FontSize', 10, ...
            'FontWeight', 'normal', ...
            'Units', 'pixels', ...
            'Position', [10, 5, wWin-30, 25]);

    %%%%%%%%%%%%%%%%%%%%%
    % Panel for class 1 %
    %%%%%%%%%%%%%%%%%%%%%
        
    uiObj.panelClass1  = uipanel(uiObj.mainWin, ...
            'Units', 'pixels', ...
            'Title', 'Class 1', ...
            'FontSize', 11, ...
            'FontWeight', 'bold', ...
            'Position',[5 hWin-270-200 wWin-10 200]);
        
    uiObj.class1NameLbl = uicontrol(uiObj.panelClass1, ...
            'Style','text', ...
            'String', 'Name:', ...
            'HorizontalAlignment','left', ...
            'FontSize', 10, ...
            'FontWeight', 'normal', ...
            'Units', 'pixels', ...
            'Position', [10, 150, 60, 20]);

    uiObj.class1Name = uicontrol(uiObj.panelClass1, ...
            'Style','edit', ...
            'String', 'class1', ...
            'HorizontalAlignment','left', ...
            'FontSize', 10, ...
            'FontWeight', 'normal', ...
            'Units', 'pixels', ...
            'Position', [60, 150, 390, 25]);

	uiObj.class1Btn = uicontrol(uiObj.panelClass1, ...
            'Style','pushbutton', ...
            'String', 'Select images', ...
            'HorizontalAlignment','left', ...
            'FontSize', 10, ...
            'FontWeight', 'normal', ...
            'Callback',{'sam_gui','selectFiles',{'1'}}, ...
            'Units', 'pixels', ...
            'Position', [460, 150, 120, 28]);
        
	uiObj.class1Files = uicontrol(uiObj.panelClass1, ...
            'Style','list', ...
            'String', '', ...
            'HorizontalAlignment','left', ...
            'FontSize', 10, ...
            'FontWeight', 'normal', ...
            'Units', 'pixels', ...
            'Position', [10, 35, wWin-30, 105]);
        
    uiObj.class1NumFiles = uicontrol(uiObj.panelClass1, ...
            'Style','text', ...
            'String', '0 images selected', ...
            'HorizontalAlignment','left', ...
            'FontSize', 10, ...
            'FontWeight', 'normal', ...
            'Units', 'pixels', ...
            'Position', [10, 5, wWin-30, 25]);
        


    %%%%%%%%%%%%%%%%%%%%%%%%
    % Panel for parameters %
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    hPanel = 105;
    uiObj.panelParams  = uipanel(uiObj.mainWin, ...
            'Units', 'pixels', ...
            'Title', 'Parameters', ...
            'FontSize', 12, ...
            'FontWeight', 'bold', ...
            'Position',[5 hWin-480-hPanel wWin-10 hPanel]);
        
%     uiObj.atlasNameLbl = uicontrol(uiObj.panelParams, ...
%             'Style','text', ...
%             'String', 'Atlas:', ...
%             'HorizontalAlignment','left', ...
%             'FontSize', 10, ...
%             'FontWeight', 'normal', ...
%             'Units', 'pixels', ...
%             'Position', [10, hPanel-55, 50, 20]);
% 
%     uiObj.atlasFile = uicontrol(uiObj.panelParams, ...
%             'Style','edit', ...
%             'String', '(no file selected)', ...
%             'HorizontalAlignment','left', ...
%             'FontSize', 10, ...
%             'FontWeight', 'normal', ...
%             'Enable','off', ...
%             'Units', 'pixels', ...
%             'Position', [60, hPanel-55, 390, 25]);
% 
% 	uiObj.selectAltasBtn = uicontrol(uiObj.panelParams, ...
%             'Style','pushbutton', ...
%             'String', 'Select atlas', ...
%             'HorizontalAlignment','left', ...
%             'FontSize', 10, ...
%             'FontWeight', 'normal', ...
%             'Callback',{'sam_gui','selectAtlas',{}}, ...
%             'Units', 'pixels', ...
%             'Position', [460, hPanel-55, 120, 28]);

    uiObj.fsNameLbl = uicontrol(uiObj.panelParams, ...
            'Style','text', ...
            'String', 'Feature selection method:', ...
            'HorizontalAlignment','left', ...
            'FontSize', 10, ...
            'FontWeight', 'normal', ...
            'Units', 'pixels', ...
            'Position', [10, hPanel-55, 200, 20]);

    uiObj.fsMethod = uicontrol(uiObj.panelParams, ...
            'Style','popupmenu', ...
            'String', fsMethods, ...
            'Value', 2, ...
            'HorizontalAlignment','left', ...
            'FontSize', 10, ...
            'FontWeight', 'normal', ...
            'Units', 'pixels', ...
            'Position', [180, hPanel-55, 140, 25]);

    uiObj.fsMaxRegLbl = uicontrol(uiObj.panelParams, ...
            'Style','text', ...
            'String', '# Regions (max.):', ...
            'HorizontalAlignment','right', ...
            'FontSize', 10, ...
            'FontWeight', 'normal', ...
            'Units', 'pixels', ...
            'Position', [330, hPanel-55, 110, 20]);

    uiObj.fsMaxReg = uicontrol(uiObj.panelParams, ...
            'Style','edit', ...
            'String', '50', ...
            'HorizontalAlignment','right', ...
            'FontSize', 10, ...
            'FontWeight', 'normal', ...
            'Units', 'pixels', ...
            'Position', [450, hPanel-55, 130, 25]);

    uiObj.feNameLbl = uicontrol(uiObj.panelParams, ...
            'Style','text', ...
            'String', 'Feature extraction method:', ...
            'HorizontalAlignment','left', ...
            'FontSize', 10, ...
            'FontWeight', 'normal', ...
            'Units', 'pixels', ...
            'Position', [10, hPanel-90, 200, 20]);

    uiObj.feMethod = uicontrol(uiObj.panelParams, ...
            'Style','popupmenu', ...
            'String', feMethods, ...
            'HorizontalAlignment','left', ...
            'FontSize', 10, ...
            'FontWeight', 'normal', ...
            'Units', 'pixels', ...
            'Position', [180, hPanel-90, 140, 25]);

    uiObj.feCompLbl = uicontrol(uiObj.panelParams, ...
            'Style','text', ...
            'String', '# Components:', ...
            'HorizontalAlignment','right', ...
            'FontSize', 10, ...
            'FontWeight', 'normal', ...
            'Units', 'pixels', ...
            'Position', [330, hPanel-90, 110, 20]);

    uiObj.feComp = uicontrol(uiObj.panelParams, ...
            'Style','edit', ...
            'String', '1', ...
            'HorizontalAlignment','right', ...
            'FontSize', 10, ...
            'FontWeight', 'normal', ...
            'Units', 'pixels', ...
            'Position', [450, hPanel-90, 130, 25]);
       
%     uiObj.outputNameLbl = uicontrol(uiObj.panelParams, ...
%             'Style','text', ...
%             'String', 'Save SAM map as:', ...
%             'HorizontalAlignment','left', ...
%             'FontSize', 10, ...
%             'FontWeight', 'normal', ...
%             'Units', 'pixels', ...
%             'Position', [10, hPanel-155, 120, 20]);
% 
%     uiObj.outputFile = uicontrol(uiObj.panelParams, ...
%             'Style','edit', ...
%             'String', '(no file selected)', ...
%             'HorizontalAlignment','left', ...
%             'FontSize', 10, ...
%             'FontWeight', 'normal', ...
%             'Enable','off', ...
%             'Units', 'pixels', ...
%             'Position', [130, hPanel-155, 320, 25]);
% 
% 	uiObj.selectOutputBtn = uicontrol(uiObj.panelParams, ...
%             'Style','pushbutton', ...
%             'String', 'Select file', ...
%             'HorizontalAlignment','left', ...
%             'FontSize', 10, ...
%             'FontWeight', 'normal', ...
%             'Callback',{'sam_gui','selectOutput',{}}, ...
%             'Units', 'pixels', ...
%             'Position', [460, hPanel-155, 120, 28]);
% 
% 	uiObj.chkShowBounds = uicontrol(uiObj.panelParams, ...
%             'Style','checkbox', ...
%             'String', 'Show bounds', ...
%             'HorizontalAlignment','left', ...
%             'Value',1, ...
%             'FontSize', 10, ...
%             'FontWeight', 'normal', ...
%             'Units', 'pixels', ...
%             'Position', [10, hPanel-180, 120, 20]);
% 
%     uiObj.chkShowMap = uicontrol(uiObj.panelParams, ...
%             'Style','checkbox', ...
%             'String', 'Show map', ...
%             'HorizontalAlignment','left', ...
%             'Value',1, ...
%             'FontSize', 10, ...
%             'FontWeight', 'normal', ...
%             'Units', 'pixels', ...
%             'Position', [130, hPanel-180, 120, 20]);

    %%%%%%%%%%%%%%%%%%%%%
    % Start button      %
    %%%%%%%%%%%%%%%%%%%%%
        
	uiObj.startBtn = uicontrol(uiObj.mainWin, ...
            'Style','pushbutton', ...
            'String', 'Start analysis', ...
            'HorizontalAlignment','left', ...
            'FontSize', 12, ...
            'FontWeight', 'bold', ...
            'Callback',{'sam_gui','start',{}}, ...
            'Units', 'pixels', ...
            'Position', [10, 10, wWin-20, 40]);
        
        
    appData = [];
    appData.gui = 1;
    appData.fsMethods = fsMethods;
    appData.feMethods = feMethods;
    appData.uiObj = uiObj;
	guidata(uiObj.mainWin, appData);
        
end

function selectFiles(h, classNum)
    appData = guidata(h);
    if isfield(appData, 'lastPath'), p = appData.lastPath; else, p = pwd; end
    [files, path] = uigetfile([p filesep '*.nii'], ...
                    'Select Class 0 images','MultiSelect', 'on');
                
    if ~isequal(files, 0)
        hList = eval(['appData.uiObj.class' classNum 'Files']);
        lis = {};
        for i=1:numel(files), lis{i} = fullfile(path, files{i}); end
        set(hList, 'String', lis, 'Value',1);
        hText = eval(['appData.uiObj.class' classNum 'NumFiles']);
        set(hText, 'String', [num2str(numel(files)) ' files selected']);
        
        appData.(['class' classNum]).files = lis;
        appData.lastPath = path;
        guidata(h, appData);
    end
end

function start(h)
    appData = guidata(h);
    
    % Check class 0 is defined
    if ~isfield(appData, 'class0')
        msgbox('No files selected for class 0')
        return
    end

    % Check class 1 is defined
    if ~isfield(appData, 'class1')
        msgbox('No files selected for class 1')
        return
    end
    
    % Check FS params are correct
    appData.method.fs = appData.fsMethods{get(appData.uiObj.fsMethod, 'Value')};
    fsMaxReg = str2double(get(appData.uiObj.fsMaxReg, 'String'));
    if isnan(fsMaxReg)
        msgbox('Invalid parameter: # Regions (max.)')
        return
    end
    appData.method.fsMaxReg = double(uint16(fsMaxReg));
    
    % Check FE params are correct
    appData.method.fe = appData.feMethods{get(appData.uiObj.feMethod, 'Value')};
    feComp = str2double(get(appData.uiObj.feComp, 'String'));
    if isnan(feComp)
        msgbox('Invalid parameter: # Components')
        return
    end
    appData.method.feComp = double(uint16(feComp));
        
    sam(appData);

end

