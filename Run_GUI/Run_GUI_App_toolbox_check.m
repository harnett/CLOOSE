classdef Run_GUI_App < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        figure1           matlab.ui.Figure
        uipanel6_6        matlab.ui.container.Panel
        uibuttongroup8_4  matlab.ui.container.ButtonGroup
        radiobutton12_4   matlab.ui.control.RadioButton
        radiobutton14_4   matlab.ui.control.RadioButton
        uipanel6_5        matlab.ui.container.Panel
        uibuttongroup8_3  matlab.ui.container.ButtonGroup
        radiobutton12_3   matlab.ui.control.RadioButton
        radiobutton14_3   matlab.ui.control.RadioButton
        text8             matlab.ui.control.Label
        text13            matlab.ui.control.Label
        edit13            matlab.ui.control.EditField
        edit7             matlab.ui.control.EditField
        uipanel6_4        matlab.ui.container.Panel
        text21_3          matlab.ui.control.Label
        edit20_3          matlab.ui.control.EditField
        edit20_2          matlab.ui.control.EditField
        text21_2          matlab.ui.control.Label
        text21            matlab.ui.control.Label
        text20            matlab.ui.control.Label
        edit21            matlab.ui.control.EditField
        edit20            matlab.ui.control.EditField
        ROIActivityPanel  matlab.ui.container.Panel
        axes3             matlab.ui.control.UIAxes
        uipanel6_3        matlab.ui.container.Panel
        uibuttongroup8_2  matlab.ui.container.ButtonGroup
        radiobutton15_2   matlab.ui.control.RadioButton
        radiobutton12_2   matlab.ui.control.RadioButton
        radiobutton14_2   matlab.ui.control.RadioButton
        uipanel6_2        matlab.ui.container.Panel
        text10_5          matlab.ui.control.Label
        edit9_5           matlab.ui.control.EditField
        text10_4          matlab.ui.control.Label
        edit9_4           matlab.ui.control.EditField
        text10_3          matlab.ui.control.Label
        edit9_3           matlab.ui.control.EditField
        text10_2          matlab.ui.control.Label
        edit9_2           matlab.ui.control.EditField
        slider1           matlab.ui.control.Slider
        uipanel6          matlab.ui.container.Panel
        checkbox13        matlab.ui.control.CheckBox
        pushbutton12      matlab.ui.control.Button
        pushbutton11      matlab.ui.control.Button
        pushbutton10      matlab.ui.control.Button
        uibuttongroup8    matlab.ui.container.ButtonGroup
        radiobutton15     matlab.ui.control.RadioButton
        radiobutton12     matlab.ui.control.RadioButton
        radiobutton13     matlab.ui.control.RadioButton
        radiobutton14     matlab.ui.control.RadioButton
        text19            matlab.ui.control.Label
        text18            matlab.ui.control.Label
        edit19            matlab.ui.control.EditField
        edit18            matlab.ui.control.EditField
        text17            matlab.ui.control.Label
        edit17            matlab.ui.control.EditField
        text15            matlab.ui.control.Label
        edit14            matlab.ui.control.EditField
        pushbutton6       matlab.ui.control.Button
        uipanel4          matlab.ui.container.Panel
        edit29            matlab.ui.control.EditField
        text30            matlab.ui.control.Label
        edit26            matlab.ui.control.EditField
        text26            matlab.ui.control.Label
        text25            matlab.ui.control.Label
        edit25            matlab.ui.control.EditField
        text12            matlab.ui.control.Label
        edit11            matlab.ui.control.EditField
        checkbox6         matlab.ui.control.CheckBox
        text11            matlab.ui.control.Label
        edit10            matlab.ui.control.EditField
        uipanel2          matlab.ui.container.Panel
        text27            matlab.ui.control.Label
        edit27            matlab.ui.control.EditField
        edit30            matlab.ui.control.EditField
        text31            matlab.ui.control.Label
        checkbox12        matlab.ui.control.CheckBox
        edit28            matlab.ui.control.EditField
        text29            matlab.ui.control.Label
        checkbox11        matlab.ui.control.CheckBox
        text28            matlab.ui.control.Label
        edit23            matlab.ui.control.EditField
        pushbutton15      matlab.ui.control.Button
        pushbutton14      matlab.ui.control.Button
        text22            matlab.ui.control.Label
        edit22            matlab.ui.control.EditField
        text10            matlab.ui.control.Label
        edit9             matlab.ui.control.EditField
        text2             matlab.ui.control.Label
        pushbutton2       matlab.ui.control.Button
        axes1             matlab.ui.control.UIAxes
        axes4             matlab.ui.control.UIAxes
    end

    
    methods (Access = private)
        function pushbutton3_Callback(app, hObject, eventdata, handles)
            % hObject    handle to pushbutton2 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            % Here we call some default settings for setting up Psychtoolbox
            
            saveTrials = get(handles.checkbox6,'Value');
            if saveTrials == 0
                asw = questdlg('You are set not to save this trial. Are you sure this is right?');
                switch asw
                    case 'Yes'
                        disp ('As you command. I am only here to serve to my true master')
                    case 'No'
                        return
                    case 'cancel'
                        return
                end
            end
            directory = get(handles.edit11, 'String');
            if saveTrials == 1
                if strcmp('D:\Experiments\SessionName',directory)
                    answer = questdlg('You are using the default directory to save your data. Are you sure this is right?');
                    switch answer
                        case 'Yes'
                            disp ('As you command. I am only here to serve to my true master')
                        case 'No'
                            return
                        case 'cancel'
                            return
                    end
                end
            end
            
            Files = dir(fullfile([directory '\BCI_Trials'], '*.mat'));
            rnd = randi([1 size(Files,1)]);
            file2load = load([directory '\BCI_Trials\' Files(rnd).name]);
            
            PsychDefaultSetup(2);
            
            % Seed the random number generator. Here we use the an older way to be
            % compatible with older systems. Newer syntax would be rng('shuffle'). Look
            % at the help function of rand "help rand" for more information
            rand('seed', sum(100 * clock));
            
            % Screen Number
            screenNumber = 2;%min(Screen('Screens'));
            
            % Define black, white and grey
            white = WhiteIndex(screenNumber);
            grey = white / 2;
            black = BlackIndex(screenNumber);
            
            % Open the screen
            [window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey, [], 32, 2,...
                [], [],  kPsychNeed32BPCFloat);
            
            % Flip to clear
            Screen('Flip', window);
            
            % Query the frame duration
            ifi = Screen('GetFlipInterval', window);
            
            % Maximum priority level
            topPriorityLevel = MaxPriority(window);
            
            % Get the centre coordinate of the window
            [xCenter, yCenter] = RectCenter(windowRect);
            
            %--------------------
            % Gabor information
            %--------------------
            
            % Dimensions
            % gaborDimPix = 500;
            gaborDimPixX = windowRect(3);
            gaborDimPixY = windowRect(4);
            
            % Sigma of Gaussian
            sigma = gaborDimPixX;
            
            % Obvious Parameters
            orientation = 90;
            contrast = 1;
            aspectRatio = 5;
            
            % Spatial Frequency (Cycles Per Pixel)
            % One Cycle = Grey-Black-Grey-White-Grey i.e. One Black and One White Lobe
            if ~isempty(get(handles.edit7,'string'))
                numCycles = str2double(get(handles.edit7,'string'));
            else
                warning('Spatial Frequency automatically set to 5')
                numCycles = 5;
            end
            freq = numCycles / gaborDimPixX;
            
            
            % Build a procedural gabor texture
            gabortex = CreateProceduralGabor(window, gaborDimPixX, gaborDimPixY,...
                [], [1 1 1 0.5], 1, 100);
            
            % Positions of the Gabors
            dim = 8;
            [x, y] = meshgrid(-dim:dim, -dim:dim);
            
            % Calculate the distance in "Gabor numbers" of each gabor from the center
            % of the array
            dist = sqrt(x.^2 + y.^2);
            
            % Cut out an inner annulus
            innerDist = 3.5;
            x(dist <= innerDist) = nan;
            y(dist <= innerDist) = nan;
            
            % Cut out an outer annulus
            outerDist = 10;
            x(dist >= outerDist) = nan;
            y(dist >= outerDist) = nan;
            
            % Select only the finite values
            x = x(isfinite(x));
            y = y(isfinite(y));
            
            % Center the annulus coordinates in the centre of the screen
            xPos = x .* gaborDimPixX + xCenter;
            yPos = y .* gaborDimPixY + yCenter;
            
            % Count how many Gabors there are
            nGabors = numel(xPos);
            
            % Make the destination rectangles for all the Gabors in the array
            baseRect = [0 0 gaborDimPixX gaborDimPixY];
            allRects = nan(4, nGabors);
            for i = 1:nGabors
                allRects(:, i) = CenterRectOnPointd(baseRect, xPos(i), yPos(i));
            end
            
            % Drift speed for the 2D global motion
            degPerSec = 360 * 2;
            degPerFrame =  degPerSec * ifi;
            
            % Randomise the Gabor orientations and determine the drift speeds of each gabor.
            % This is given by multiplying the global motion speed by the cosine
            % difference between the global motion direction and the global motion.
            % Here the global motion direction is 0. So it is just the cosine of the
            % angle we use. We re-orientate the array when drawing
            gaborAngles = rand(1, nGabors) .* 180 - 90;
            degPerFrameGabors = cosd(gaborAngles) .* degPerFrame;
            
            % Randomise the phase of the Gabors and make a properties matrix. We could
            % if we want have each Gabor with different properties in all dimensions.
            % Not just orientation and drift rate as we are doing here.
            % This is the power of using procedural textures
            phaseLine = rand(1, nGabors) .* 360;
            propertiesMat = repmat([NaN, freq, sigma, contrast, aspectRatio, 0, 0, 0],...
                nGabors, 1);
            propertiesMat(:, 1) = phaseLine';
            
            % Perform initial flip to gray background and sync us to the retrace:
            vbl = Screen('Flip', window);
            
            % Numer of frames to wait before re-drawing
            waitframes = 1;
            
            % angle to show
            angle = 1:360;
            n = 0;
            % Animation loop
            % while ~KbCheck
            targetAngle = file2load.StimData.targetAngle;
            
            xx = mod((targetAngle - 90), 360);
            rotation = 0;
            cla(handles.axes1,'reset');
            h = animatedline(handles.axes1);
            axis(handles.axes1,[0,1500,0,90]);
            
            trig = get(handles.checkbox10,'Value');
            if trig == 1
                ServerSend = tcpip('10.93.6.2',50000,'NetworkRole','server', ...
                    'OutputBufferSize', 2, 'Terminator', 'CR/LF');
                fopen(ServerSend);
                fwrite(ServerSend, uint16(0), 'uint16');
                fclose(ServerSend);
                delete(ServerSend);
                ServerSend = [];
            end
            
            
            StimData.('TStampGlobal') = nan(1, length(file2load.StimData.Angle));
            StimData.('Angle')  = nan(1, length(file2load.StimData.Angle));
            StimData.('Rotation')  = nan(1, length(file2load.StimData.Angle) );
            StimData.('targetAngle') = targetAngle;
            
            timerVal = tic;
            
            for i = 1: length(file2load.StimData.Angle)
                % Set the right blend function for drawing the gabors
                Screen('BlendFunction', window, 'GL_ONE', 'GL_ZERO');
            
                % Batch draw all of the Gabors to screen
                xx = file2load.StimData.Angle(i);
                rotation = mod((targetAngle - xx), 360);
                Screen('DrawTexture', window, gabortex, [], [], xx,...
                    [], [], [], [], kPsychDontDoRotation, propertiesMat');
            
                addpoints(h, i, rotation)
                drawnow limitrate nocallbacks
            
                StimData.TStampGlobal(i) = toc(timerVal);
                StimData.Angle(i)  = xx;
                StimData.Rotation(i)  = rotation;
                Screen('Flip', window);
            end
            if xx >= targetAngle
                StimData.Reward = 1;
            elseif xx < targetAngle
                StimData.Reward = 0;
            end
            StimData.numCycles = numCycles;
            
            hold off
            if saveTrials == 1
                currentCounterValue = str2double(get(handles.edit10, 'String'));
                newString = sprintf('%d', int32(currentCounterValue +1));
                set(handles.edit10, 'String', newString );
            
                trialNumber = currentCounterValue;
                Folder = '\Playback';
                filename = ['\' num2str(trialNumber)];
                if ~exist([directory Folder], 'dir')
                    mkdir([directory Folder])
                end
                save([directory Folder filename], 'StimData')
            end
            % Clean up
            sca;
        end
        
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function Run_GUI_OpeningFcn(app, varargin)
            % Ensure that the app appears on screen when run
            movegui(app.figure1, 'onscreen');
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app); %#ok<ASGLU>
            
            % This function has no output args, see OutputFcn.
            % hObject    handle to figure
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            % varargin   command line arguments to Run_GUI (see VARARGIN)
            
            % handles.peaks=peaks(35);
            % handles.membrane=membrane;
            % [x,y] = meshgrid(-8:.5:8);
            % r = sqrt(x.^2+y.^2) + eps;
            % sinc = sin(r)./r;
            % handles.sinc = sinc;
            % % Set the current data value.
            % handles.current_data = handles.peaks;
            % surf(handles.current_data)
            
            % Choose default command line output for Run_GUI
            handles.output = hObject;
            handles.ROIdrawNum = 0;
            % Update handles structure
            guidata(hObject, handles);
        end

        % Value changed function: edit17
        function edit17_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to edit17 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: get(hObject,'String') returns contents of edit17 as text
            %        str2double(get(hObject,'String')) returns contents of edit17 as a double
            val = str2double(get(hObject,'String'));
        end

        % Button pushed function: pushbutton10
        function pushbutton10_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to pushbutton10 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            answer = questdlg('Which ROI set do you want me to load?', ...
                'Options', ...
                'Last saved set','Manual load', 'Manual load');
            directory = get(handles.edit11, 'String');
            Folder = '\ROIs';
            listing = dir([directory Folder '\*.mat']);
            table = struct2table(listing);
            sortedT = sortrows(table, 'date');
            sortedS = table2struct(sortedT);
            
            switch answer
                case 'Last saved set'
                    if sortedS(end).isdir == 0
                        load([directory Folder '\' sortedS(end).name]);
                    elseif sortedS(end).isdir == 1
                        if sortedS(end - 1).isdir == 1
                            load([directory Folder '\' sortedS(end-2).name]);
                        elseif sortedS(end - 1).isdir == 0
                            load([directory Folder '\' sortedS(end-1).name]);
                        end
            
                    end
                case 'Manual load'
                    uiload
            end
            % pixelsX = 796;
            % LinesY  = 512;
            % h = handles.axes4;
            % mask = cell(length(ROILoad.ROIpos),1);
            % for iroi = 1:length(ROILoad.ROIpos)
            %     tb = nan(pixelsX,LinesY);
            %     for ix = 1:pixelsX
            %         xq = repmat(ix, LinesY,1);
            %         yq = [1:LinesY]';
            %         tb(ix,:) = inpolygon(xq,yq, ROILoad.ROIpos(1,iroi).Position(:,1),...
            %         ROILoad.ROIpos(1,iroi).Position(:,2));
            %     end
            %     mask{iroi,1} = tb;
            %     [r, b] = find(mask{iroi, 1});
            %     scatter(h, r, b);
            %     hold(h, 'on')
            % end
            
            h = handles.axes4;
            % colormap(h, gray);
            hold (h, 'on')
            for i = 1:numel(RoiInfo.ROIpos)
                handles.ROI(i) = drawpolygon(h, 'Position',RoiInfo.ROIpos(1,i).Position,...
                    'LineWidth', 1, 'Deletable', true, 'Label', num2str(i), ...
                    'LabelVisible', 'hover', 'Color', 'y', 'FaceAlpha', 0);
            end
            handles.nR = numel(handles.ROI);
            handles.ROIdrawNum = length(listing);
            guidata(hObject, handles);
        end

        % Button pushed function: pushbutton11
        function pushbutton11_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to pushbutton11 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            h = handles.axes4;
            handles.nR = handles.nR +1;
            nR = handles.nR;
            guidata(hObject, handles);
            handles.ROI(nR) = drawpolygon(h, 'LineWidth', 1, 'Deletable', true, ...
                'Label', num2str(nR), 'LabelVisible', 'hover', 'Color', 'y', ...
                'FaceAlpha', 0);
            guidata(hObject, handles);
        end

        % Button pushed function: pushbutton12
        function pushbutton12_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to pushbutton12 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            handles.ROIdrawNum = handles.ROIdrawNum + 1;
            guidata(hObject, handles);
            
            isOk = isvalid(handles.ROI);
            RoiInfo.ROIpos = handles.ROI(1,isOk);

            RoiInfo.subpop = 1:length(RoiInfo.ROIpos);
                        
            directory = get(handles.edit11, 'String');
            if strcmp('D:\Experiments\SessionName',directory)
                answer = questdlg('You are using the default directory to save your data. Are you sure this is right?');
                switch answer
                    case 'Yes'
                        disp ('As you command. I am only here to serve to my true master')
                    case 'No'
                        return
                    case 'cancel'
                        return
                end
            end
            Folder = '\ROIs';
            nN = handles.ROIdrawNum;
            filename = ['\' num2str(nN)];
            if ~exist([directory Folder], 'dir')
                mkdir([directory Folder])
            end
            RoiInfo.dir = [directory Folder filename];
            save([directory Folder filename], 'RoiInfo')
        end

        % Button pushed function: pushbutton14
        function pushbutton14_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to pushbutton14 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA
            % hObject    handle to pushbutton2 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            % Here we call some default settings for setting up Psychtoolbox
            PsychDefaultSetup(2);
            
            % Seed the random number generator. Here we use the an older way to be
            % compatible with older systems. Newer syntax would be rng('shuffle'). Look
            % at the help function of rand "help rand" for more information
            rand('seed', sum(100 * clock));
            Screen('Preference','SkipSyncTests', 1);
            % Screen Number
            screenNumber = 2;%min(Screen('Screens'));
            
            % Define black, white and grey
            white = WhiteIndex(screenNumber);
            grey = white / 2;
            black = BlackIndex(screenNumber);
            
            % Open the screen
            [window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey, [], 32, 2,...
                [], [],  kPsychNeed32BPCFloat);
            
            % Flip to clear
            Screen('Flip', window);
            
            % Query the frame duration
            ifi = Screen('GetFlipInterval', window);
            
            % Maximum priority level
            topPriorityLevel = MaxPriority(window);
            
            % Get the centre coordinate of the window
            [xCenter, yCenter] = RectCenter(windowRect);
            
            %--------------------
            % Gabor information
            %--------------------
            
            % Dimensions
            % gaborDimPix = 500;
            gaborDimPixX = windowRect(3);
            gaborDimPixY = windowRect(4);
            
            % Sigma of Gaussian
            sigma = gaborDimPixX/8;
            
            % Obvious Parameters
            orientation = 90;
            contrast = 1;
            aspectRatio = 1;
            
            % Spatial Frequency (Cycles Per Pixel)
            % One Cycle = Grey-Black-Grey-White-Grey i.e. One Black and One White Lobe
            if ~isempty(get(handles.edit7,'string'))
                numCycles = str2double(get(handles.edit7,'string'));
            else
                warning('Spatial Frequency automatically set to 5')
                numCycles = 5;
            end
            freq = numCycles / gaborDimPixX;
            
            
            % Build a procedural gabor texture
            gabortex = CreateProceduralGabor(window, gaborDimPixX, gaborDimPixY,...
                [], [0.5 0.5 0.5 1], 1, 10);
            
            % Positions of the Gabors
            dim = 8;
            [x, y] = meshgrid(-dim:dim, -dim:dim);
            
            % Calculate the distance in "Gabor numbers" of each gabor from the center
            % of the array
            dist = sqrt(x.^2 + y.^2);
            
            % Cut out an inner annulus
            innerDist = 3.5;
            x(dist <= innerDist) = nan;
            y(dist <= innerDist) = nan;
            
            % Cut out an outer annulus
            outerDist = 10;
            x(dist >= outerDist) = nan;
            y(dist >= outerDist) = nan;
            
            % Select only the finite values
            x = x(isfinite(x));
            y = y(isfinite(y));
            
            % Center the annulus coordinates in the centre of the screen
            xPos = x .* gaborDimPixX + xCenter;
            yPos = y .* gaborDimPixY + yCenter;
            
            % Count how many Gabors there are
            nGabors = numel(xPos);
            
            % Make the destination rectangles for all the Gabors in the array
            baseRect = [0 0 gaborDimPixX gaborDimPixY];
            allRects = nan(4, nGabors);
            for i = 1:nGabors
                allRects(:, i) = CenterRectOnPointd(baseRect, xPos(i), yPos(i));
            end
            
            % Drift speed for the 2D global motion
            degPerSec = 360 * 2;
            degPerFrame =  degPerSec * ifi;
            
            % Randomise the Gabor orientations and determine the drift speeds of each gabor.
            % This is given by multiplying the global motion speed by the cosine
            % difference between the global motion direction and the global motion.
            % Here the global motion direction is 0. So it is just the cosine of the
            % angle we use. We re-orientate the array when drawing
            gaborAngles = rand(1, nGabors) .* 180 - 90;
            degPerFrameGabors = cosd(gaborAngles) .* degPerFrame;
            
            % Randomise the phase of the Gabors and make a properties matrix. We could
            % if we want have each Gabor with different properties in all dimensions.
            % Not just orientation and drift rate as we are doing here.
            % This is the power of using procedural textures
            phaseLine = rand(1, nGabors) .* 360;
            propertiesMat = repmat([NaN, freq, sigma, contrast, aspectRatio, 0, 0, 0],...
                nGabors, 1);
            propertiesMat(:, 1) = phaseLine';
            
            % Perform initial flip to gray background and sync us to the retrace:
            vbl = Screen('Flip', window);
            
            % Numer of frames to wait before re-drawing
            waitframes = 1;
            
            % angle to show
            angle = 1:360;
            n = 1;
            % Animation loop
            % while ~KbCheck
            if ~isempty(get(handles.edit13,'string'))
                targetAngle = str2double(get(handles.edit13,'string'));
            else
                error('Select a target angle')
            end
            
            xx = mod((targetAngle - 90), 360);
            angle_bins  = mod((targetAngle - 90), 360) : 10 : targetAngle;
            
            
            cla(handles.axes1,'reset');
            h = plot(nan(1,2), nan,'Parent',handles.axes1,'Color', 'k', 'LineWidth', 1);
            axis(handles.axes1,[-90,0,0,90]);
            set(h(1:end),'XData', -89:0);
            
            
            timerVal = tic;
            
            pixelsX = str2double(get(handles.edit21,'string'));
            LinesY  = str2double(get(handles.edit20,'string'));
            FramesZ = str2double(get(handles.edit22,'string'));
            tduration = FramesZ/60;
            idx = nan(1, FramesZ);
            Debug = get(handles.checkbox8,'Value');
            
            rotation = nan(1, FramesZ);
            rotation(1) = 0;
            if Debug == 0
            
                ServerSend = tcpip('10.93.6.2',50000,'NetworkRole','server', ...
                    'OutputBufferSize', 2, 'Terminator', 'CR/LF');
                fopen(ServerSend);
                fwrite(ServerSend, uint16(FramesZ), 'uint16');
                fclose(ServerSend);
            
                rawData = nan((pixelsX*LinesY),FramesZ);
                tcpipServer = tcpip('10.93.6.2',30000,'NetworkRole','client', 'Terminator', 'CR/LF');
                tcpipServer.InputBufferSize = (pixelsX*LinesY*2);
                fopen(tcpipServer);
            end
            
            if Debug == 1
                while toc(timerVal) < tduration && rotation(n) < 90
                    % Set the right blend function for drawing the gabors
                    n = n+1;
                    if n == 1
                        Screen('BlendFunction', window, 'GL_ONE', 'GL_ZERO');
                        Screen('DrawTexture', window, gabortex, [], [], 45,...
                            [], [], [], [], kPsychDontDoRotation, propertiesMat');
                    end
                    [x, ~, ~] = GetMouse(window);
                    ang = round((x/1920)*90) - mod((90-targetAngle),360);
                    idx(n) = max(find(ang >= angle_bins));
            
                    xx = angle_bins(idx(n));
                    rotation(n) = mod((targetAngle - xx), 360);
            
                    if n ~= 1 && idx(n) ~= idx(n -1)
                        Screen('DrawTexture', window, gabortex, [], [], xx,...
                            [], [], [], [], kPsychDontDoRotation, propertiesMat');
                         Screen('Flip', window, 0 ,1);
                    end
            
                    if n < 90
                        plpt = fliplr(rotation(1:n));
                        plpt(end+1 : 90) = nan;
                        plpt = fliplr(plpt);
                    else
                        plpt = rotation(n-89 : n);
                    end
                    set(h(1),'YData',plpt);
                    set(h(2),'YData',plpt.*0.5);
            %         addpoints(h, n, rotation(n))
                    drawnow limitrate nocallbacks
                    pause(0.03)
            %         if n == 10
            %             pause()
            %         end
            %         Screen('AsyncFlipBegin', window, 0, 2);
            
            %         n = n+1;
                end
            elseif Debug == 0
                while n < FramesZ && rotation <= 90
                    n = n+1;
                    Screen('BlendFunction', window, 'GL_ONE', 'GL_ZERO');
                    rawData(:,n) = fread(tcpipServer,(pixelsX*LinesY),'uint16');
            
                    x = nanmean(rawData(:,n));
                    xx = round((x/max(rawData(:,n)))*90) - mod((90-targetAngle),360);
                    rotation = mod((targetAngle - xx), 360);
                    Screen('DrawTexture', window, gabortex, [], [], xx,...
                        [], [], [], [], kPsychDontDoRotation, propertiesMat');
                    addpoints(h, n, rotation)
                    drawnow limitrate nocallbacks
            
                    Screen('Flip', window);%, vbl + (waitframes - 0.5) * ifi);
                end
                fclose(tcpipServer);
                delete(tcpipServer);
                delete(ServerSend);
                ServerSend = [];
            end
            hold off
            
            figure(10);
            scatter((1:20),rand(1,20));
            % Clean up
            sca;
        end

        % Button pushed function: pushbutton15
        function pushbutton15_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to pushbutton15 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            global tts vx vy
            
            hp = handles.axes4;
            hp.Interactions = [zoomInteraction regionZoomInteraction rulerPanInteraction];
            hp.Toolbar = [];
            colormap(hp, gray);
            
            mc = get(handles.checkbox12,'Value');
            if mc == 1
                fft_tts = nan(vy(1, 2) - vy(1, 1) + 1, vx(1, 2) - vx(1, 1) + 1, 3);
                % vx1 = [vx vx*2-1; vx vx*2-1; vx*2 vx*3-1; vx*2 vx*3-1];
                % vy1 = [vy vy*2-1; vy*2 vy*3-1; vy vy*2-1; vy*2 vy*3-1];
                vx1 = vy; %[vx vx*3-1; vx vx*3-1];
                vy1 = vx; %[vy vy*2-1; vy*2 vy*3-1];
                for i = 1:size(vx1,1)
                    fft_tts(:,:,i) = fft2(tts(vx1(i,1):vx1(i,2), vy1(i,1):vy1(i,2)));
                end
            end
            
            
            pixelsX = str2double(get(handles.edit21,'string'));
            LinesY  = str2double(get(handles.edit20,'string'));
            FramesZ = str2double(get(handles.edit27,'string'));
            oneTrial_length = str2double(get(handles.edit22,'string'));
            fold = str2double(get(handles.edit30,'string'));
            if fold ~= 0
                fold_fct = LinesY./fold;
            end
            
            answer = questdlg('Which ROI set do you want me to use to estimate Baseline?', ...
                'Options', 'Last saved set','Manual load', 'Manual load');
            directory = get(handles.edit11, 'String');
            Folder = '\ROIs';
            listing = dir([directory Folder '\*.mat']);
            lst = 0;
            switch answer
                case 'Last saved set'
                    load([directory Folder '\' listing(end).name]);
                    lst =1;
                case 'Manual load'
                    uiload
            end
            
            mask = cell(length(RoiInfo.ROIpos),1);
            meanRoiOneFrame = nan(length(RoiInfo.ROIpos), FramesZ);
            
            baseline.F0 = nan(length(RoiInfo.ROIpos), 1);
            for iroi = 1:length(RoiInfo.ROIpos)
                tb = nan(LinesY, pixelsX);
                for ix = 1:LinesY
                    yq = repmat(ix, pixelsX,1);
                    xq = [1:pixelsX]';
                    tb(ix,1:pixelsX) = inpolygon(xq, yq, RoiInfo.ROIpos(1,iroi).Position(:,1),...
                    RoiInfo.ROIpos(1,iroi).Position(:,2));
                end
                mask{iroi,1} = logical(tb);
            end
            
            plane_num = str2double(get(handles.edit18,'string'));
            plane_tot = str2double(get(handles.edit19,'string'));
            frames_sampled = plane_num:plane_tot:(FramesZ*plane_tot);
            baseline.Licks = nan(1, FramesZ);
            baseline.Opto = nan(1, FramesZ);
            baseline.timeStamp = nan(1, FramesZ);
            baseline.angle = uint8(ones(1, FramesZ));
            baseline.ImTranslateXY = nan(2, FramesZ);
            
            Screen('Preference', 'SkipSyncTests', 1);
            PsychDefaultSetup(2);
            rand('seed', sum(100 * clock));
            screenNumber = 2;%min(Screen('Screens'));
            white = WhiteIndex(screenNumber);
            grey = white / 2;
            black = BlackIndex(screenNumber);
            waitframes = 1;
            % Open the screen
            [window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey, [], 32, 2,...
                [], [],  kPsychNeed32BPCFloat);
            
            Screen('Flip', window);
            ifi = Screen('GetFlipInterval', window);
            [xCenter, yCenter] = RectCenter(windowRect);
            %--------------------
            % Gabor information
            %--------------------
            gaborDimPixX = windowRect(3);
            gaborDimPixY = windowRect(4);
            % Sigma of Gaussian
            sigma = gaborDimPixX/8 ;
            % Obvious Parameters
            orientation = 90;
            contrast = 1;
            aspectRatio = 1;
            % Spatial Frequency (Cycles Per Pixel)
            % One Cycle = Grey-Black-Grey-White-Grey i.e. One Black and One White Lobe
            if ~isempty(get(handles.edit5,'string'))
                numCycles = str2double(get(handles.edit5,'string'));
            else
                warning('Spatial Frequency automatically set to 5')
                numCycles = 5;
            end
            freq = numCycles / gaborDimPixX;
            % Build a procedural gabor texture
            gabortex = CreateProceduralGabor(window, gaborDimPixX, gaborDimPixY,...
                [], [0.5 0.5 0.5 1], 1, 10);
            % Positions of the Gabors
            dim = 8;
            [x, y] = meshgrid(-dim:dim, -dim:dim);
            % Calculate the distance in "Gabor numbers" of each gabor from the center
            % of the array
            dist = sqrt(x.^2 + y.^2);
            % Cut out an inner annulus
            innerDist = 3.5;
            x(dist <= innerDist) = nan;
            y(dist <= innerDist) = nan;
            % Cut out an outer annulus
            outerDist = 10;
            x(dist >= outerDist) = nan;
            y(dist >= outerDist) = nan;
            % Select only the finite values
            x = x(isfinite(x));
            y = y(isfinite(y));
            % Center the annulus coordinates in the centre of the screen
            xPos = x .* gaborDimPixX + xCenter;
            yPos = y .* gaborDimPixY + yCenter;
            % Count how many Gabors there are
            nGabors = numel(xPos);
            % Make the destination rectangles for all the Gabors in the array
            baseRect = [0 0 gaborDimPixX gaborDimPixY];
            allRects = nan(4, nGabors);
            for i = 1:nGabors
                allRects(:, i) = CenterRectOnPointd(baseRect, xPos(i), yPos(i));
            end
            
            % Drift speed for the 2D global motion
            if ~isempty(get(handles.edit6,'string'))
                tFreq = str2double(get(handles.edit6,'string'));
            else
                warning('Temporal Frequency automatically set to 2')
                tFreq = 2;
            end
            degPerSec = 360 * tFreq;
            degPerFrame =  degPerSec * ifi;
            % Randomise the Gabor orientations and determine the drift speeds of each gabor.
            % This is given by multiplying the global motion speed by the cosine
            % difference between the global motion direction and the global motion.
            % Here the global motion direction is 0. So it is just the cosine of the
            % angle we use. We re-orientate the array when drawing
            gaborAngles = 0; % rand(1, nGabors) .* 180 - 90;
            degPerFrameGabors = cosd(gaborAngles) .* degPerFrame;
            % Randomise the phase of the Gabors and make a properties matrix. We could
            % if we want have each Gabor with different properties in all dimensions.
            % Not just orientation and drift rate as we are doing here.
            % This is the power of using procedural textures
            phaseLine = 0.5;% rand(1, nGabors) .* 360;
            propertiesMat = repmat([NaN, freq, sigma, contrast, aspectRatio, 0, 0, 0],...
                nGabors, 1);
            propertiesMat(:, 1) = phaseLine';
            angles = 0:15:90;
            if ~isempty(get(handles.edit4,'string'))
                tdrift = str2double(get(handles.edit4,'string'));
            else
                tdrift = 1;
            end

            
            ServerSend = tcpip('10.93.6.2',50000,'NetworkRole','server', ...
                'OutputBufferSize', 12, 'Terminator', 'CR/LF');
            fopen(ServerSend);
            fwrite(ServerSend, uint32([(FramesZ*plane_tot)  plane_num  plane_tot]), 'uint32');
            fclose(ServerSend);
            
            % rawData = nan((pixelsX*LinesY),FramesZ);
            tcpipServer = tcpip('10.93.6.2',30000,'NetworkRole','client', 'Terminator', 'CR/LF');
            tcpipServer.InputBufferSize = (pixelsX*LinesY*200);
            fopen(tcpipServer);
            
            StimBaseline = get(handles.checkbox11,'Value');
            
            Lucky_angle = 1:8;
            Lucky_pick = 2;
            count_lucks = 0;
            n = 0;
            pl = 0;
            
            imin = -0.5;
            imx = 0.5;
            
            while pl < FramesZ*plane_tot
            
                if StimBaseline == 1
                    stm = 0;
                    pick_angle = angles(randperm(length(angles),1));
                    timerVal = tic;
                    trand = imin + (imx-imin).*rand;
                    while toc(timerVal) < tdrift + trand
                        pl = pl+1;
                        if pl > FramesZ*plane_tot
                            break
                        end
                        rawData = fread(tcpipServer,(pixelsX*LinesY),'uint16');
            
                        if sum(ismember(pl, frames_sampled)) == 0
                            continue
                        end
                        
                        n = n + 1;
                        stm = stm + 1;
                        if n == 1
                            timeStamp = tic;
                            rawData_tmp = uint16(reshape(rawData, pixelsX, LinesY));
                            rawData_tmp = permute(rawData_tmp, [2 1]);
                            rawData_tmp = imcomplement(rawData_tmp);
                            if fold ~= 0
                                rawData_tmp = reshape(rawData_tmp, pixelsX, fold, fold_fct);
                                rawData_tmp = uint16(mean(rawData_tmp, 3));
                            end
                            hp = imshow(rawData_tmp, [], 'InitialMagnification', 800, 'Parent', hp, 'border', 'tight');
                        end
                        rawData = reshape(uint16(rawData), pixelsX, LinesY);
                        rawData = permute(rawData, [2 1]);
                        rawData = imcomplement(rawData);
                        if fold ~= 0
                            r_tmp = uint16(zeros(fold, pixelsX, fold_fct));
                            for ifc = 1:fold_fct %#ok<BDSCI>
                                fuck_index = (fold*(ifc-1))+ 1 : fold*ifc;
                                r_tmp(:, :, ifc) = rawData(fuck_index, :, 1);
                            end
                            rawData = uint16(mean(r_tmp, 3));
                        end
            
                        if mc == 1
                            output = nan(3, 4);
                            for i = 1:size(fft_tts, 3)
                                output(i,:) = dftregistration(fft_tts(:,:,i), fft2(double(rawData(vx1(i,1):vx1(i,2), vy1(i,1):vy1(i,2)))));
                            end
                            xsh = round(nanmedian(output(:,4)));
                            ysh = round(nanmedian(output(:,3)));
                            rawData = imtranslate(rawData, [xsh ysh]);
                            baseline.ImTranslateXY(:, n)  = [xsh; ysh];
            
                        end
            
                        for iroi = 1:length(mask)
                            meanRoiOneFrame(iroi,n) = nanmean(rawData(mask{iroi,1}));
                        end
            
                        set(hp, 'CData', im2uint16(rawData));
                        try
                            drawnow limitrate nocallbacks
                        catch
                            warning('Could not plot')
                        end
            
                        if stm == 1
                            Screen('BlendFunction', window, 'GL_ONE', 'GL_ZERO');
                            % Draw texture
                            Screen('DrawTexture', window, gabortex, [], [], pick_angle,...
                                [], [], [], [], kPsychDontDoRotation, propertiesMat');
                            Screen('Flip', window, 0, 1);
                        end
                        %Screen('Flip', window);
                        %
                        baseline.timeStamp(n) = toc(timeStamp);
                        baseline.angle(n) = uint8(pick_angle);

                    end

                elseif StimBaseline == 0
                    pl = pl+1;
                    if pl > FramesZ*plane_tot
                        break
                    end
                    rawData = fread(tcpipServer,(pixelsX*LinesY),'uint16');
            
                    if sum(ismember(pl, frames_sampled)) == 0
                        continue
                    end
            
                    n = n + 1;
                    if n == 1
                        timeStamp = tic;
                        rawData_tmp = uint16(reshape(rawData, pixelsX, LinesY));
                        rawData_tmp = permute(rawData_tmp, [2 1]);
                        rawData_tmp = imcomplement(rawData_tmp);
                        if fold ~= 0
                            rawData_tmp = reshape(rawData_tmp, pixelsX, fold, fold_fct);
                            rawData_tmp = uint16(mean(rawData_tmp, 3));
                        end
                        hp = imshow(rawData_tmp, [], 'InitialMagnification', 800, 'Parent', hp, 'border', 'tight');
                    end
                    rawData = reshape(uint16(rawData), pixelsX, LinesY);
                    rawData = permute(rawData, [2 1]);
                    rawData = imcomplement(rawData);
                    if fold ~= 0
                        r_tmp = uint16(zeros(fold, pixelsX, fold_fct));
                        for ifc = 1:fold_fct %#ok<BDSCI>
                            fuck_index = (fold*(ifc-1))+ 1 : fold*ifc;
                            r_tmp(:, :, ifc) = rawData(fuck_index, :, 1);
                        end
                        rawData = uint16(mean(r_tmp, 3));
                    end
            
                    if mc == 1
                        output = nan(3, 4);
                        for i = 1:size(fft_tts, 3)
                            output(i,:) = dftregistration(fft_tts(:,:,i), fft2(double(rawData(vx1(i,1):vx1(i,2), vy1(i,1):vy1(i,2)))));
                        end
                        xsh = round(nanmedian(output(:,4)));
                        ysh = round(nanmedian(output(:,3)));
                        rawData = imtranslate(rawData, [xsh ysh]);
                        baseline.ImTranslateXY(:, n)  = [xsh; ysh];
                    end
            
                    for iroi = 1:length(mask)
                        meanRoiOneFrame(iroi,n) = nanmean(rawData(mask{iroi,1}));
                    end
            
                    set(hp, 'CData', im2uint16(rawData));
                    try
                        drawnow limitrate nocallbacks
                    catch
                        warning('Could not plot')
                    end
            
                    Screen('FillRect', window, [0 0 0]);
                    Screen('Flip', window, 0, 1);
            
                    baseline.timeStamp(n) = toc(timeStamp);
                    baseline.angle(n) = uint8(0);
                end
            end
            
            %clean up
            sca;
            
            fclose(tcpipServer);
            
            delete(tcpipServer);
            delete(ServerSend);
            ServerSend = [];
            
            meanRoiOneFrame = double(meanRoiOneFrame);
            
            figure(1)
            for i = 1:size(meanRoiOneFrame, 1)
                plot(meanRoiOneFrame(i,:)./max(meanRoiOneFrame(i,:)) + i, 'k')
                hold on
            end
            hold off
            disp('Press Any Button to Continue')
            
            try
                pause
            catch
                pause
            end
            
            subPopulations = inputdlg({'SubPop1 [x:n,y,z]','SubPop2 [x:n,y,z]'},...
                          'ROI Selection', [1 35; 1 35]);
            
            close figure 1
            
            sp1 = str2num(subPopulations{1});
            sp2 = str2num(subPopulations{2});
            subpop1_2 = cat(2, sp1, sp2);
            RoiInfo.ROIpos = RoiInfo.ROIpos(subpop1_2);
            
            RoiInfo.subpop = 1:length(sp1);
            RoiInfo.subpop2 = length(sp1)+1:length(sp1)+length(sp2);
            save(RoiInfo.dir, 'RoiInfo');
            
            baseline.subpop1 = RoiInfo.subpop;
            baseline.subpop2 = RoiInfo.subpop2;
            
            
            Roi_s1 = meanRoiOneFrame(sp1, :);
            Roi_s2 = meanRoiOneFrame(sp2, :);
            meanRoiOneFrame = cat(1, Roi_s1, Roi_s2);
            
            hp = handles.axes4;
            himage = imshow(uint16(tts), [], 'InitialMagnification', 800, 'Parent', hp, 'border', 'tight');
            set(handles.slider1, 'Min', 1, 'Max', 1, 'Value', 1, ...
                'SliderStep', [0 0]);
            hold (hp, 'on')
            for i = 1:numel(RoiInfo.ROIpos)
                if sum(ismember(i, baseline.subpop1)) > 0
                    handles.ROI(i) = drawpolygon(hp, 'Position',RoiInfo.ROIpos(1,i).Position,...
                        'LineWidth', 1, 'Deletable', true, 'Label', num2str(i), ...
                        'LabelVisible', 'hover', 'Color', 'r', 'FaceAlpha', 0);
                elseif sum(ismember(i, baseline.subpop2)) > 0
                    handles.ROI(i) = drawpolygon(hp, 'Position',RoiInfo.ROIpos(1,i).Position,...
                        'LineWidth', 1, 'Deletable', true, 'Label', num2str(i), ...
                        'LabelVisible', 'hover', 'Color', 'g', 'FaceAlpha', 0);
                end
            end
            % guidata(hObject, handles.ROI);
            
            % baseline.fullrec = rawData;
            baseline.mean = mean(meanRoiOneFrame,2);
            baseline.std  = std(meanRoiOneFrame,0,2);
            baseline.roiData = meanRoiOneFrame;
            baseline.roiMask = mask(subpop1_2);
            baseline.nroi = length(RoiInfo.ROIpos);
            baseline.stimBaseline = StimBaseline;
            
            if lst == 1
                baseline.roiFile = [directory Folder '\' listing(end).name];
            elseif lst ==0
                baseline.roiFile = 'Manual Pick';
            end
            indx = true(1, size(meanRoiOneFrame, 2));

            trial_start = randperm(size(baseline.roiData(:, indx),2), 200);
            df_oneTrial  = cell(1, length(trial_start));
            z_zz = cell(1, length(trial_start));
            F0_all = cell(1, length(trial_start));
            F_2 = repmat(baseline.roiData(:, indx), 1, 2);
            F0_tmp = nan(size(baseline.roiData(:, indx)));
            
            for i = 1:size(baseline.roiData(:, indx), 2)
                F0_tmp(:,i) = prctile(F_2(:,i:i+(oneTrial_length*4)-1), 10, 2);
            end
            F0_tmp= (repmat(F0_tmp, 1, 2));
            
            for itrial = 1:length(trial_start)
                beg_end = [trial_start(itrial)+size(baseline.roiData(:, indx),2)-749; ...
                    trial_start(itrial)+size(baseline.roiData(:, indx),2)];
            
                F0_all{itrial} = F0_tmp(:,beg_end(1):beg_end(2));
                df_oneTrial{itrial} =  (F_2(:, trial_start(itrial) : trial_start(itrial) + 749)-F0_all{itrial})./F0_all{itrial};
                if ~isempty(baseline.subpop2)
                    z_zz{itrial}= nanmean(df_oneTrial{itrial}(baseline.subpop1, :),1) - nanmean(df_oneTrial{itrial}(baseline.subpop2,:),1);
                else
                    z_zz{itrial}= nanmean(df_oneTrial{itrial},1);
                end
            end
            
            if ~isempty(baseline.subpop2)
                dfmx1 = max(cat(2, z_zz{:}));
                dfmx2 = min(cat(2, z_zz{:}));
                threshs = dfmx2: 0.005: dfmx1;
            else
                threshs = 0:0.01:max(df, [], 'all');
            end
            
            baseline.zzz = z_zz;
            Day = str2double(get(handles.edit29,'string'));
            if Day == 1
                asw = questdlg('Is this really Day 1?');
                switch asw
                    case 'Yes'
                        dayy1 = 1;
                    case 'No'
                        dayy1 = 0;
                end
            else
                asw = questdlg('Are you sure this IS NOT day 1?');
                switch asw
                    case 'Yes'
                        dayy1 = 0;
                    case 'No'
                        dayy1 = 1;
                end
            end
            
            if dayy1 == 1
                propp = nan(1,length(threshs));
                for i = 1:length(threshs)
                    total_large = nan(1,length(z_zz));
                    for ii =  1:length(z_zz)
                        mvMu = movmean(z_zz{ii}, 3);
                        %         mvMu = z_zz{ii};
                        if sum(mvMu > threshs(i)) > 0
                            total_large(ii) = 1;
                        else
                            total_large(ii) = 0;
                        end
                    end
                    propp(1,i) = sum(total_large)/length(total_large);
                    propp(2,i) = threshs(i);
                end
                baseline.propTrialsTresh = propp;
                Thresh.propTrialsTresh = baseline.propTrialsTresh;
                baseline.selected_tresh = propp(2,min(find(propp(1,:) <= 0.5)));
                [muFit, sigmaFit] = normfit(cat(2, z_zz{:}));
                baseline.muFit = muFit;
                baseline.sigmaFit = sigmaFit;
                Thresh.muFit = muFit;
                Thresh.sigmaFit = sigmaFit;
                z_sc = (baseline.selected_tresh - muFit)/sigmaFit;
                % baseline.selected_range = muFit - (z_sc*sigmaFit);
                baseline.selected_range = min(cat(2, z_zz{:}));
            
                cc =  cat(2, baseline.zzz{:});
                % Data on the left side of the distribution
                d = cc(cc <= nanmedian(cc));
                % Mirror the distribution on the right
                dd = max(d)-d;
                %concatenate the left and right side of the distribution
                ddc = cat(2,d, dd);
                %fit
                [baseline.muHalfFit, baseline.sigmaHalfFit] = normfit(ddc);
                % z score
                baseline.z_HalfFit = (baseline.selected_tresh - baseline.muHalfFit)/baseline.sigmaHalfFit;
            
            elseif dayy1 ~= 1
            
                propp = nan(1,length(threshs));
                for i = 1:length(threshs)
                    total_large = nan(1,length(z_zz));
                    for ii =  1:length(z_zz)
                        mvMu = movmean(z_zz{ii}, 3);
                        %         mvMu = z_zz{ii};
                        if sum(mvMu > threshs(i)) > 0
                            total_large(ii) = 1;
                        else
                            total_large(ii) = 0;
                        end
                    end
                    propp(1,i) = sum(total_large)/length(total_large);
                    propp(2,i) = threshs(i);
                end
                selected_tresh = propp(2,min(find(propp(1,:) <= 0.5)));
                [muFit, sigmaFit] = normfit(cat(2, z_zz{:}));
                z_sc = (selected_tresh - muFit)/sigmaFit;
                selected_range = min(cat(2, z_zz{:}));
            
                cc_2 =  cat(2, baseline.zzz{:});
                % Data on the left side of the distribution
                d_2 = cc_2(cc_2 <= nanmedian(cc_2));
                % Mirror the distribution on the right
                dd_2 = max(d_2)-d_2;
                %concatenate the left and right side of the distribution
                ddc_2 = cat(2,d_2, dd_2);
                %fit
                [muHalfFit, sigmaHalfFit] = normfit(ddc_2);
                % z score
                z_HalfFit = (selected_tresh - muHalfFit)/sigmaHalfFit;
            
                uiload
            
                if Thresh.z_HalfFit >= z_HalfFit
                    baseline.z_HalfFit = z_HalfFit;
                    baseline.selected_tresh = selected_tresh;
                    baseline.selected_range = selected_range;
                    baseline.muHalfFit = muHalfFit;
                    baseline.sigmaHalfFit = sigmaHalfFit;
                    disp ('Using today baseline')
                else
                    cc =  cat(2, baseline.zzz{:});
                    % Data on the left side of the distribution
                    d = cc(cc <= nanmedian(cc));
                    % Mirror the distribution on the right
                    dd = max(d)-d;
                    %concatenate the left and right side of the distribution
                    ddc = cat(2,d, dd);
                    %fit
                    [baseline.muHalfFit, baseline.sigmaHalfFit] = normfit(ddc);
                    baseline.z_HalfFit = Thresh.z_HalfFit;
                    baseline.selected_tresh = baseline.muHalfFit + (baseline.sigmaHalfFit*Thresh.z_HalfFit);
                    baseline.selected_range = min(cat(2, z_zz{:}));
                    disp ('Using uploaded baseline')
                end
            end
            
            baseline.F_all = F0_all;
            for iroi = 1:size(baseline.roiData,1)
                baseline.F0(iroi,1) = prctile(meanRoiOneFrame(iroi,end-869:end), 10);
            end
            
            baseline.RefImage = tts;
            
            Thresh.selected_tresh = baseline.selected_tresh;
            Thresh.selected_range = baseline.selected_range;
            Thresh.z_HalfFit = baseline.z_HalfFit;
            
            baseline.plane_num = plane_num;
            baseline.plane_tot = plane_tot;
            
            saveTrials = get(handles.checkbox6,'Value');
            if saveTrials == 0
                asw = questdlg('You are set not to save this trial. Are you sure this is right?');
                switch asw
                    case 'Yes'
                        disp ('As you command. I am only here to serve to my true master')
                    case 'No'
                        return
                    case 'cancel'
                        return
                end
            elseif saveTrials == 1
                if strcmp('D:\Experiments\SessionName',directory)
                    answer = questdlg('You are using the default directory to save your data. Are you sure this is right?');
                    switch answer
                        case 'Yes'
                            disp ('As you command. I am only here to serve to my true master')
                        case 'No'
                            return
                        case 'cancel'
                            return
                    end
                end
                currentCounterValue = str2double(get(handles.edit10, 'String'));
                % Create a new string with the number being 1 more than the current number.
                newString = sprintf('%d', int32(currentCounterValue +1));
                % Send the new string to the text control.
                set(handles.edit10, 'String', newString );
                %     directory = get(handles.edit11, 'String');
                Folder = '\Baseline';
                filename1 = ['\baselineMask_' num2str(currentCounterValue)];
                filename2 = ['\Threshold_values' num2str(currentCounterValue)];
                if ~exist([directory Folder], 'dir')
                    mkdir([directory Folder])
                end
                save([directory Folder filename1], 'baseline')
                save([directory Folder filename2], 'Thresh')
            end
        end

        % Callback function
        function pushbutton16_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to pushbutton16 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            hpp = handles.axes4;
            hpp.Interactions = [zoomInteraction regionZoomInteraction rulerPanInteraction];
            hpp.Toolbar = [];
            colormap(hpp, gray);
            
            saveTrials = get(handles.checkbox6,'Value');
            if saveTrials == 0
                asw = questdlg('You are set not to save this trial. Are you sure this is right?');
                switch asw
                    case 'Yes'
                        disp ('As you command. I am only here to serve to my true master')
                    case 'No'
                        return
                    case 'cancel'
                        return
                end
            end
            directory = get(handles.edit11, 'String');
            if saveTrials == 1
                if strcmp('D:\Experiments\SessionName',directory)
                    answer = questdlg('You are using the default directory to save your data. Are you sure this is right?');
                    switch answer
                        case 'Yes'
                            disp ('As you command. I am only here to serve to my true master')
                        case 'No'
                            return
                        case 'cancel'
                            return
                    end
                end
            end
            
            pixelsX = str2double(get(handles.edit21,'string'));
            LinesY  = str2double(get(handles.edit20,'string'));
            FramesZ = str2double(get(handles.edit22,'string'));
            
            answer = questdlg('Which ROI set do you want me to use to estimate Baseline?', ...
                'Options', 'Last saved set','Manual load', 'Manual load');
            directory = get(handles.edit11, 'String');
            Folder = '\ROIs';
            listing = dir([directory Folder '\*.mat']);
            lst = 0;
            switch answer
                case 'Last saved set'
                    load([directory Folder '\' listing(end).name]);
                    lst =1;
                case 'Manual load'
                    uiload
            end
            
            mask = cell(length(RoiInfo.ROIpos),1);
            % meanRoiOneFrame = nan(length(RoiInfo.ROIpos), FramesZ);
            % baseline.F0 = nan(length(RoiInfo.ROIpos), 1);
            for iroi = 1:length(RoiInfo.ROIpos)
                tb = nan(LinesY, pixelsX);
                for ix = 1:LinesY
                    yq = repmat(ix, pixelsX,1);
                    xq = [1:pixelsX]';
                    tb(ix,1:pixelsX) = inpolygon(xq, yq, RoiInfo.ROIpos(1,iroi).Position(:,1),...
                    RoiInfo.ROIpos(1,iroi).Position(:,2));
                end
                mask{iroi,1} = logical(tb);
            end
            roiMask = mask;
            nroi = length(RoiInfo.ROIpos);
            
            n = 0;
            
            AngleNumb = str2double(get(handles.uibuttongroup1.SelectedObject,'String'));
            if AngleNumb == 8
                angles = 0:45:315;
            elseif AngleNumb == 12
                angles = 0:30:330;
            end
            
            if get(handles.checkbox6,'Value') == 1
                blduration  = str2double(get(handles.edit3, 'String'));
            end
            if get(handles.checkbox2,'Value') == 1
                drduration  = str2double(get(handles.edit4, 'String'));
            end
            if get(handles.checkbox3,'Value') == 1
                grduration  = str2double(get(handles.edit2, 'String'));
            end
            
            
            rsmpl  = 1;
            numCycles = str2double(get(handles.edit5, 'String'));
            if isnan(numCycles)
                numCycles = 4;
            end
            % Stim.SpFreq = numCycles;
            
            
            % Stimulus presentation
            PsychDefaultSetup(2);
            screenNumber = 2;%min(Screen('Screens'));
            % Screen('Preference', 'SkipSyncTests', 1);
            
            white = WhiteIndex(screenNumber);
            grey = white / 2;
            [window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey, [], 32, 2,...
                [], [],  kPsychNeed32BPCFloat);
            Screen('Flip', window);
            
            ifi = Screen('GetFlipInterval', window);
            topPriorityLevel = MaxPriority(window);
            
            
            gaborDimPixX = 400;
            % Stim.GaborSz = gaborDimPixX;
            % Stim.type = 'Gabor';
            
            sigma = gaborDimPixX / 10;
            gabortex = CreateProceduralGabor(window, gaborDimPixX, gaborDimPixX,...
                [], [0.5 0.5 0.5 0], 1, 100);
            
            xPos = windowRect(3)/8:windowRect(3)/4:windowRect(3);
            yPos = windowRect(4)/6:windowRect(4)/3:windowRect(4);
            xPos = repmat(xPos, 1, numel(yPos));
            xPos = sort(xPos);
            yPos = repmat(yPos, 1, numel(xPos)/numel(yPos));
            % Count how many Gabors there are
            nGabors = numel(xPos);
            
            % Stim.NGabor = rsmpl;
            
            % Make the destination rectangles for all the Gabors in the array
            baseRect = [0 0 gaborDimPixX gaborDimPixX];
            allRects = nan(4, nGabors);
            for i = 1:nGabors
                allRects(:, i) = CenterRectOnPointd(baseRect, xPos(i), yPos(i));
            end
            
            contrast = 1;
            aspectRatio = 1;
            freq = numCycles / gaborDimPixX;
            
            tFreq = str2double(get(handles.edit6, 'String'));
            degPerSec = 360 * tFreq;
            degPerFrame =  degPerSec * ifi;
            
            % Randomise the Gabor orientations and determine the drift speeds of each gabor.
            % This is given by multiplying the global motion speed by the cosine
            % difference between the global motion direction and the global motion.
            % Here the global motion direction is 0. So it is just the cosine of the
            % angle we use. We re-orientate the array when drawing
            gaborAngles = rand(1, nGabors) .* 180 - 90;
            degPerFrameGabors = cosd(gaborAngles) .* degPerFrame;
            
            % Randomise the phase of the Gabors and make a properties matrix. We could
            % if we want have each Gabor with different properties in all dimensions.
            % Not just orientation and drift rate as we are doing here.
            % This is the power of using procedural textures
            phaseLine = rand(1, nGabors) .* 360;
            propertiesMat = repmat([NaN, freq, sigma, contrast, aspectRatio, 0, 0, 0],...
                nGabors, 1);
            propertiesMat(:, 1) = phaseLine';
            
            vbl = Screen('Flip', window);
            
            waitframes = 1;
            
            loop = 0;
            
            trig = get(handles.checkbox9,'Value');
            % l = 1;
            % loops = str2double(get(handles.edit2, 'String'));
            % Stim.Loops = loops;
            
            RandStream('dsfmt19937');
            
            % y = datasample(1:loops,loops, 'Replace', false);
            % p = str2double(get(handles.edit5, 'String'));
            % nv = loops*(p/100);
            % idx  = find(y <= nv);
            % angles(idx, 2) = -1;
            roi_one_frame = nan(nroi,FramesZ);
            roi_resp = nan(nroi, 500);
            x_y_gabor_coord = nan(1, 500);
            
            ServerSend = tcpip('10.93.6.2',50000,'NetworkRole','server', ...
                'OutputBufferSize', 2, 'Terminator', 'CR/LF');
            fopen(ServerSend);
            fwrite(ServerSend, uint16(FramesZ), 'uint16');
            fclose(ServerSend);
            
            tcpipServer = tcpip('10.93.6.2',30000,'NetworkRole','client', 'Terminator', 'CR/LF');
            tcpipServer.InputBufferSize = (pixelsX*LinesY*10);
            fopen(tcpipServer);
            
            Screen('BlendFunction', window, 'GL_ONE', 'GL_ZERO');
            
            while n <= FramesZ
            
                loop = loop+1;
            
                angles = datasample(angles, numel(angles), 'Replace', false);
                rnd_idx = datasample(1:size(allRects,2), rsmpl , 'Replace', false);
                rnd_gabors = allRects(:, rnd_idx);
            
                tb  = tic;
                cb = 0;
                while toc(tb) < blduration
                    n = n+1;
                    if n > FramesZ
                        break
                    end
                    rawData = uint16(fread(tcpipServer,(pixelsX*LinesY),'uint16'));
                    rawData_tmp = uint16(reshape(rawData, pixelsX, LinesY));
                    rawData_tmp = permute(rawData_tmp, [2 1]);
                    rawData_tmp = imcomplement(rawData_tmp);
                    if n == 1
                        hpp = image(im2uint16(rawData_tmp)); %, [], 'InitialMagnification', 800, 'Parent', hpp, 'border', 'tight');
                    end
                    set(hpp, 'CData',im2uint16(rawData_tmp));
                    drawnow limitrate nocallbacks
            
                    for ii = 1:nroi
                        roi_one_frame(ii,n) = double(mean(rawData_tmp(roiMask{ii, 1})));
                    end
                    cb = cb + 1;
                    if cb == 1
                        Screen('FillRect', window, [0 0 0]);
                        Screen('Flip', window, 0, 1);
                    end
                end
            
                tg  = tic;
                cg = 0;
                while toc(tg) < grduration
                    n = n+1;
                    if n > FramesZ
                        break
                    end
                    rawData = uint16(fread(tcpipServer,(pixelsX*LinesY),'uint16'));
                    rawData_tmp = uint16(reshape(rawData, pixelsX, LinesY));
                    rawData_tmp = permute(rawData_tmp, [2 1]);
                    rawData_tmp = imcomplement(rawData_tmp);
                    for ii = 1:nroi
                        roi_one_frame(ii,n) = double(mean(rawData_tmp(roiMask{ii, 1})));
                    end
                    if n == 1
                        % hpp = image(im2uint8(rawData_tmp));
                        hpp = imshow(rawData_tmp, [], 'InitialMagnification', 800, 'Parent', hpp, 'border', 'tight');
                    end
                    set(hpp, 'CData',im2uint16(rawData_tmp));
                    drawnow limitrate nocallbacks
            
            
                    cg = cg + 1;
                    if cg == 1
                        Screen('FillRect', window, [0.5 0.5 0.5]);
                        Screen('Flip', window, 0, 1);
                    end
                end
            
                count = 0;
                for i = 1: size(angles, 2)
                    timerval = tic;
                    count = count +1;
                    if count == 1
                        n_start  = n;
                    end
                    while toc(timerval) < drduration
                        n = n+1;
                        if n > FramesZ
                            break
                        end
                        rawData = uint16(fread(tcpipServer,(pixelsX*LinesY),'uint16'));
                        rawData_tmp = uint16(reshape(rawData, pixelsX, LinesY));
                        rawData_tmp = permute(rawData_tmp, [2 1]);
                        rawData_tmp = imcomplement(rawData_tmp);
                        for ii = 1:nroi
                            roi_one_frame(ii,n) = double(mean(rawData_tmp(roiMask{ii, 1})));
                        end
            
                        if n == 1
                            % hpp = image(im2uint8(rawData_tmp));
                            hpp = imshow(rawData_tmp, [], 'InitialMagnification', 800, 'Parent', hpp, 'border', 'tight');
                        end
                        set(hpp, 'CData',im2uint16(rawData_tmp));
                        drawnow limitrate nocallbacks
            
                        Screen('BlendFunction', window, 'GL_ONE', 'GL_ZERO');
                        Screen('DrawTextures', window, gabortex, [], rnd_gabors, angles(i), ...
                            [], [], [], [], kPsychDontDoRotation, propertiesMat(rnd_idx, :)');
                        Screen('Flip', window);%, vbl + (waitframes - 0.5) * ifi);
                        % Increment the phase of our Gabors
                        phaseLine = phaseLine + degPerFrameGabors;
                        propertiesMat(:, 1) = phaseLine';
                    end
                end
                if n > FramesZ
                    n_end = FramesZ;
                else
                    n_end = n;
                end
                roi_resp(:, loop) = nanmean(roi_one_frame(:, n_start:n_end), 2);
                x_y_gabor_coord(:, loop) = rnd_idx;
            end
            
            F0 = nan(size(roi_one_frame, 1), 1);
            for iroi = 1:size(roi_one_frame, 1)
                F0(iroi,1) = double(prctile(roi_one_frame(iroi,:), 10));
            end
            df = (roi_resp -F0)./F0;
            
            retinotopy.df = df; retinotopy.F0 = F0; retinotopy.gabor = allRects;
            retinotopy.roi_F = roi_resp; retinotopy.idx_gabor = x_y_gabor_coord;
            
            cla(handles.axes3);
            hp = handles.axes3;
            hp.Interactions = [zoomInteraction regionZoomInteraction rulerPanInteraction];
            hp.Toolbar = [];
            resp = nan(1, size(retinotopy.gabor, 2));
            err = nan(1, size(retinotopy.gabor, 2));
            single_neuron_resp = nan(size(df, 1), size(retinotopy.gabor, 2));
            
            for i = 1:size(retinotopy.gabor, 2)
                all_resp_one_patch = retinotopy.df(:, retinotopy.idx_gabor == i);
                resp(i) = nanmean(all_resp_one_patch, 'all');
                err(i) = std(nanmean(all_resp_one_patch, 2))/sqrt(size(all_resp_one_patch, 1));
                single_neuron_resp(:, i) = nanmean(all_resp_one_patch, 2);
                scatter(hp, repmat(i, 1, size(all_resp_one_patch, 1)), single_neuron_resp(:, i), 'k', 'Filled')
                hold (hp, 'on')
            end
            bar(hp,resp, 'r', 'FaceAlpha',0.3)
            hold on
            errorbar(hp, resp,err, 'k', 'LineWidth', 2, 'LineStyle','none')
            retinotopy.meanResp = resp;
            retinotopy.SEM = err;
            retinotopy.all_resp_perGabor = single_neuron_resp;
            
            saveTrials = get(handles.checkbox6,'Value');
            if saveTrials == 0
                asw = questdlg('You are set not to save this trial. Are you sure this is right?');
                switch asw
                    case 'Yes'
                        disp ('As you command. I am only here to serve to my true master')
                    case 'No'
                        return
                    case 'cancel'
                        return
                end
            elseif saveTrials == 1
                if strcmp('D:\Experiments\SessionName',directory)
                    answer = questdlg('You are using the default directory to save your data. Are you sure this is right?');
                    switch answer
                        case 'Yes'
                            disp ('As you command. I am only here to serve to my true master')
                        case 'No'
                            return
                        case 'cancel'
                            return
                    end
                end
                currentCounterValue = str2double(get(handles.edit10, 'String'));
                newString = sprintf('%d', int32(currentCounterValue +1));
                set(handles.edit10, 'String', newString );
                Folder = '\Retinotopy';
                filename1 = ['\retinotopy_' num2str(currentCounterValue)];
                if ~exist([directory Folder], 'dir')
                    mkdir([directory Folder])
                end
                save([directory Folder filename1], 'retinotopy')
            end
            
            sca;
        end

        % Callback function
        function pushbutton17_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to pushbutton17 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            answer = questdlg('Which Retinotopy data do you want me to load?', ...
                'Options', ...
                'Last saved set','Manual load', 'Manual load');
            directory = get(handles.edit11, 'String');
            Folder = '\Retinotopy';
            listing = dir([directory Folder '\*.mat']);
            table = struct2table(listing);
            sortedT = sortrows(table, 'date');
            sortedS = table2struct(sortedT);
            
            switch answer
                case 'Last saved set'
                    if sortedS(end).isdir == 0
                        load([directory Folder '\' sortedS(end).name]);
                    elseif sortedS(end).isdir == 1
                        if sortedS(end - 1).isdir == 1
                            load([directory Folder '\' sortedS(end-2).name]);
                        elseif sortedS(end - 1).isdir == 0
                            load([directory Folder '\' sortedS(end-1).name]);
                        end
            
                    end
                case 'Manual load'
                    uiload
            end
            
            cla(handles.axes3);
            hp = handles.axes3;
            for i = 1:size(retinotopy.all_resp_perGabor,2)
                scatter(hp, repmat(i, 1, size(retinotopy.all_resp_perGabor, 1)), ...
                    retinotopy.all_resp_perGabor(:, i), 'k', 'Filled')
                hold (hp, 'on')
            end
            bar(hp,retinotopy.meanResp, 'r', 'FaceAlpha',0.3)
            hold on
            errorbar(hp, retinotopy.meanResp, retinotopy.SEM, 'k', 'LineWidth', 2, 'LineStyle','none')
        end

        % Callback function
        function pushbutton18_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to pushbutton18 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            clear arduino
            arduino = serialport('COM3', 9600);
            flush(arduino);
            fprintf('Connection established with Arduino\n')
            
            
            hp = handles.axes4;
            hp.Interactions = [zoomInteraction regionZoomInteraction rulerPanInteraction];
            hp.Toolbar = [];
            colormap(hp, gray);
            
            pixelsX = str2double(get(handles.edit21,'string'));
            LinesY  = str2double(get(handles.edit20,'string'));
            FramesZ = str2double(get(handles.edit27,'string'));
            
            directory = get(handles.edit11, 'String');
            listing = dir([directory '\Baseline\']);
            table = struct2table(listing);
            sortedT = sortrows(table, 'date');
            sortedS = table2struct(sortedT);
            if sortedS(end).isdir == 0
                load([directory '\Baseline\' sortedS(end).name]);
                bline_loaded = [directory '\Baseline\' sortedS(end).name];
            elseif sortedS(end).isdir == 1
                if sortedS(end - 1).isdir == 1
                    load([directory '\Baseline\' sortedS(end-2).name]);
                    bline_loaded = [directory '\Baseline\' sortedS(end-2).name];
                elseif sortedS(end - 1).isdir == 0
                    load([directory '\Baseline\' sortedS(end-1).name]);
                    bline_loaded = [directory '\Baseline\' sortedS(end-1).name];
                end
            end
            if isfield(baseline, 'fullrec')
                baseline = rmfield(baseline, 'fullrec');
            end
            F0 = double(baseline.F0);
            subpop1 = baseline.subpop1;
            subpop2 = baseline.subpop2;
            nroi = baseline.nroi;
            
            mask = baseline.roiMask;
            plane_num = str2double(get(handles.edit18,'string'));
            plane_tot = str2double(get(handles.edit19,'string'));
            frames_sampled = plane_num:plane_tot:(FramesZ*plane_tot);
            
            meanRoiOneFrame = nan(nroi, FramesZ);
            
            rotating_drift.timeStamp = nan(1, FramesZ);
            rotating_drift.angle = uint8(ones(1, FramesZ));
            rotating_drift.fwd_bkwd = uint8(repmat(11, 1, FramesZ));
            rotating_drift.Licks = uint8(ones(1, FramesZ*plane_tot));
            
            Screen('Preference', 'SkipSyncTests', 1);
            PsychDefaultSetup(2);
            rand('seed', sum(100 * clock));
            
            screenNumber = 2;%min(Screen('Screens'));
            
            white = WhiteIndex(screenNumber);
            grey = white / 2;
            black = BlackIndex(screenNumber);
            waitframes = 1;
            
            % Open the screen
            [window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey, [], 32, 2,...
                [], [],  kPsychNeed32BPCFloat);
            
            Screen('Flip', window);
            ifi = Screen('GetFlipInterval', window);
            
            [xCenter, yCenter] = RectCenter(windowRect);
            
            %--------------------
            % Gabor information
            %--------------------
            
            % Dimensions
            % gaborDimPix = 500;
            gaborDimPixX = windowRect(3);
            gaborDimPixY = windowRect(4);
            
            % Sigma of Gaussian
            sigma = gaborDimPixX ;
            
            % Obvious Parameters
            orientation = 90;
            contrast = 1;
            aspectRatio = 1;
            
            % Spatial Frequency (Cycles Per Pixel)
            % One Cycle = Grey-Black-Grey-White-Grey i.e. One Black and One White Lobe
            if ~isempty(get(handles.edit5,'string'))
                numCycles = str2double(get(handles.edit5,'string'));
            else
                warning('Spatial Frequency automatically set to 5')
                numCycles = 5;
            end
            freq = numCycles / gaborDimPixX;
            
            % Build a procedural gabor texture
            gabortex = CreateProceduralGabor(window, gaborDimPixX, gaborDimPixY,...
                [], [1 1 1 0.5], 1, 100);
            
            % Positions of the Gabors
            dim = 8;
            [x, y] = meshgrid(-dim:dim, -dim:dim);
            
            % Calculate the distance in "Gabor numbers" of each gabor from the center
            % of the array
            dist = sqrt(x.^2 + y.^2);
            
            % Cut out an inner annulus
            innerDist = 3.5;
            x(dist <= innerDist) = nan;
            y(dist <= innerDist) = nan;
            
            % Cut out an outer annulus
            outerDist = 10;
            x(dist >= outerDist) = nan;
            y(dist >= outerDist) = nan;
            
            % Select only the finite values
            x = x(isfinite(x));
            y = y(isfinite(y));
            
            % Center the annulus coordinates in the centre of the screen
            xPos = x .* gaborDimPixX + xCenter;
            yPos = y .* gaborDimPixY + yCenter;
            
            % Count how many Gabors there are
            nGabors = numel(xPos);
            
            % Make the destination rectangles for all the Gabors in the array
            baseRect = [0 0 gaborDimPixX gaborDimPixY];
            allRects = nan(4, nGabors);
            for i = 1:nGabors
                allRects(:, i) = CenterRectOnPointd(baseRect, xPos(i), yPos(i));
            end
            
            % Drift speed for the 2D global motion
            if ~isempty(get(handles.edit6,'string'))
                tFreq = str2double(get(handles.edit6,'string'));
            else
                warning('Temporal Frequency automatically set to 2')
                tFreq = 2;
            end
            degPerSec = 360 * tFreq;
            degPerFrame =  degPerSec * ifi;
            
            % Randomise the Gabor orientations and determine the drift speeds of each gabor.
            % This is given by multiplying the global motion speed by the cosine
            % difference between the global motion direction and the global motion.
            % Here the global motion direction is 0. So it is just the cosine of the
            % angle we use. We re-orientate the array when drawing
            gaborAngles = 0; % rand(1, nGabors) .* 180 - 90;
            degPerFrameGabors = cosd(gaborAngles) .* degPerFrame;
            
            % Randomise the phase of the Gabors and make a properties matrix. We could
            % if we want have each Gabor with different properties in all dimensions.
            % Not just orientation and drift rate as we are doing here.
            % This is the power of using procedural textures
            phaseLine = 0.5;% rand(1, nGabors) .* 360;
            propertiesMat = repmat([NaN, freq, sigma, contrast, aspectRatio, 0, 0, 0],...
                nGabors, 1);
            propertiesMat(:, 1) = phaseLine';
            
            waitframes = 1;
            
            angles = 0:10:90;
            
            if ~isempty(get(handles.edit4,'string'))
                tdrift = str2double(get(handles.edit4,'string'));
            else
                tdrift = 1;
            end
            
            presTimeDrift = repmat(tdrift, 1, numel(angles));
            % Animation loop
            % while ~KbCheck
            
            sent_data = 'a';
            write(arduino,sent_data, 'char');
            
            ServerSend = tcpip('10.93.6.2',50000,'NetworkRole','server', ...
                'OutputBufferSize', 12, 'Terminator', 'CR/LF');
            fopen(ServerSend);
            fwrite(ServerSend, uint32([(FramesZ*plane_tot) plane_num plane_tot]), 'uint32');
            fclose(ServerSend);
            
            
            % rawData = nan((pixelsX*LinesY),FramesZ);
            tcpipServer = tcpip('10.93.6.2',30000,'NetworkRole','client', 'Terminator', 'CR/LF');
            tcpipServer.InputBufferSize = (pixelsX*LinesY*200);
            fopen(tcpipServer);
            
            Lucky_angle = 1:8;
            Lucky_pick = 2;
            count_lucks = 0;
            n = 0;
            pl = 0;
            brk = 0;
            while pl < FramesZ*plane_tot
                stm = 0;
                if brk == 1
                    break
                end
                if (count_lucks == 0 || count_lucks == 3)
                    count_lucks = 0;
                    pick_angle = angles(randperm(length(angles),1));
                    Lucky_pick  = Lucky_angle(randperm(length(Lucky_angle),1));
                    fwd_bkwd = randi([1, 2], 1);
                    if pick_angle > 70 && Lucky_pick == 1 && fwd_bkwd == 1
                        pick_angle = 50;
                    elseif pick_angle < 20 && Lucky_pick == 1 && fwd_bkwd == 2
                        pick_angle = 40;
                    end
                elseif (count_lucks > 0 && count_lucks < 3)
                    if fwd_bkwd == 1
                        pick_angle = pick_angle + 20;
                    elseif fwd_bkwd == 2
                        pick_angle = pick_angle - 20;
                    end
                end
                timerVal = tic;
                while toc(timerVal) < tdrift
                    pl = pl+1;
                    if pl > FramesZ*plane_tot
                        break
                    end
                    rawData = fread(tcpipServer,(pixelsX*LinesY),'uint16');
                    try
                        rotating_drift.Licks(pl)  = read(arduino, 1, 'uint8');
                        flush(arduino, 'input');
                    catch
                        rotating_drift.Licks(pl)  = rotating_drift.Licks(pl-1) - 1;
                        warning ('No Info From Arduino')
                    end
                    if sum(ismember(pl, frames_sampled)) == 0
                        continue
                    end
            
                    n = n + 1;
                    stm = stm + 1;
                    if n == 1
                        timeStamp = tic;
                        rawData_tmp = uint16(reshape(rawData, pixelsX, LinesY));
                        rawData_tmp = permute(rawData_tmp, [2 1]);
                        rawData_tmp = imcomplement(rawData_tmp);
                        % hp = image(im2uint8(rawData_tmp)); %
                        hp = imshow(rawData_tmp, [], 'InitialMagnification', 800, 'Parent', hp, 'border', 'tight');
                    end
                    try
                        rawData = reshape(uint16(rawData), pixelsX, LinesY);
                    catch
                        brk = 1;
                        warning('scanbox crashed. Getting out of here to save data')
                        break
                    end
                    rawData = permute(rawData, [2 1]);
                    rawData = imcomplement(rawData);
            
                    for iroi = 1:length(mask)
                        meanRoiOneFrame(iroi,n) = nanmean(rawData(mask{iroi,1}));
                    end
            
                    set(hp, 'CData',im2uint16(rawData));
                    try
                        drawnow limitrate nocallbacks
                    catch
                        warning('Could not plot')
                    end
            
                    if stm == 1
                        Screen('BlendFunction', window, 'GL_ONE', 'GL_ZERO');
                        % Draw texture
                        Screen('DrawTexture', window, gabortex, [], [], pick_angle,...
                            [], [], [], [], kPsychDontDoRotation, propertiesMat');
                        Screen('Flip', window, 0, 1);
                    end
                    %Screen('Flip', window);
                    %
                    rotating_drift.timeStamp(n) = toc(timeStamp);
                    rotating_drift.angle(n) = uint8(pick_angle);
                    if Lucky_pick == 1
                        rotating_drift.fwd_bkwd(n) = uint8(fwd_bkwd);
                    else
                        rotating_drift.fwd_bkwd(n) = uint8(0);
                    end
                end
                if Lucky_pick == 1
                    count_lucks = count_lucks + 1;
                end
            end
            %clean up
            sca;
            
            fclose(tcpipServer);
            
            delete(tcpipServer);
            delete(ServerSend);
            ServerSend = [];
            
            % Screen('Flip', window);
            % Screen('FillRect', window, [0 0 0]);
            %
            % baseline.fullrec = rawData;
            rotating_drift.mean = uint16(mean(meanRoiOneFrame,2));
            rotating_drift.std  = uint16(std(meanRoiOneFrame,0,2));
            rotating_drift.roiData = uint16(meanRoiOneFrame);
            rotating_drift.roiMask = mask;
            rotating_drift.nroi = nroi;
            rotating_drift.subpop1 = baseline.subpop1;
            rotating_drift.subpop2 = baseline.subpop2;
            rotating_drift.F0 =F0;
            
            rotating_drift.plane_num = plane_num;
            rotating_drift.plane_tot = plane_tot;
            
            saveTrials = get(handles.checkbox6,'Value');
            if saveTrials == 0
                asw = questdlg('You are set not to save this trial. Are you sure this is right?');
                switch asw
                    case 'Yes'
                        disp ('As you command. I am only here to serve to my true master')
                    case 'No'
                        return
                    case 'cancel'
                        return
                end
            elseif saveTrials == 1
                if strcmp('D:\Experiments\SessionName',directory)
                    answer = questdlg('You are using the default directory to save your data. Are you sure this is right?');
                    switch answer
                        case 'Yes'
                            disp ('As you command. I am only here to serve to my true master')
                        case 'No'
                            return
                        case 'cancel'
                            return
                    end
                end
                currentCounterValue = str2double(get(handles.edit10, 'String'));
                % Create a new string with the number being 1 more than the current number.
                newString = sprintf('%d', int32(currentCounterValue +1));
                % Send the new string to the text control.
                set(handles.edit10, 'String', newString );
                %     directory = get(handles.edit11, 'String');
                Folder = '\rotating_drift';
                filename1 = ['\rotating_drift' num2str(currentCounterValue)];
                if ~exist([directory Folder], 'dir')
                    mkdir([directory Folder])
                end
                save([directory Folder filename1], 'rotating_drift')
            end
        end

        % Callback function
        function pushbutton19_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to pushbutton19 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            clear arduino
            arduino = serialport('COM3', 9600);%, 'ErrorOccurredFcn', mycallbackFuck);
            pause(2)
            % flush(arduino);
            sent_data = 'b';
            write(arduino,sent_data, 'char');
            
            disp('Reward Delivered');
        end

        % Button pushed function: pushbutton2
        function pushbutton2_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            
            % hObject    handle to pushbutton2 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            % Here we call some default settings for setting up Psychtoolbox
            
            global tts vx vy
            
            
            
            mc = get(handles.checkbox12,'Value');
            if mc == 1
                fft_tt64s = nan(128, 199, 3);
            
                vx1 = vy; %[vx vx*3-1; vx vx*3-1];
                vy1 = vx; %[vy vy*2-1; vy*2 vy*3-1];
                for i = 1:size(vx1,1)
                    fft_tts(:,:,i) = fft2(tts(vx1(i,1):vx1(i,2), vy1(i,1):vy1(i,2)));
                end
            %     fft_tts = uint16(fft_tts);
            end
            
            
            saveTrials = get(handles.checkbox6,'Value');
            if saveTrials == 0
                asw = questdlg('You are set not to save this trial. Are you sure this is right?');
                switch asw
                    case 'Yes'
                        disp ('As you command. I am only here to serve to my true master')
                    case 'No'
                        return
                    case 'cancel'
                        return
                end
            end
            directory = get(handles.edit11, 'String');
            if saveTrials == 1
                if strcmp('D:\Experiments\SessionName',directory)
                    answer = questdlg('You are using the default directory to save your data. Are you sure this is right?');
                    switch answer
                        case 'Yes'
                            disp ('As you command. I am only here to serve to my true master')
                        case 'No'
                            return
                        case 'cancel'
                            return
                    end
                end
            end
            
            listing = dir([directory '\Baseline\']);
            table = struct2table(listing);
            sortedT = sortrows(table, 'date');
            sortedS = table2struct(sortedT);
            
            % load([directory '\Baseline\' sortedS(end-2).name]);
            % bline_loaded = [directory '\Baseline\' sortedS(end-2).name];
            
            if sortedS(end).isdir == 0
                load([directory '\Baseline\' sortedS(end).name]);
                bline_loaded = [directory '\Baseline\' sortedS(end).name];
            elseif sortedS(end).isdir == 1
                if sortedS(end - 1).isdir == 1
                    load([directory '\Baseline\' sortedS(end-2).name]);
                    bline_loaded = [directory '\Baseline\' sortedS(end-2).name];
                elseif sortedS(end - 1).isdir == 0
                    load([directory '\Baseline\' sortedS(end-1).name]);
                    bline_loaded = [directory '\Baseline\' sortedS(end-1).name];
                end
            end
            if isfield(baseline, 'fullrec')
                baseline = rmfield(baseline, 'fullrec');
            end
            baseline.mean = double(baseline.mean);
            baseline.F0 = double(baseline.F0);
            subpop1 = baseline.subpop1;
            subpop2 = baseline.subpop2;
            
            pixelsX = str2double(get(handles.edit21,'string'));
            LinesY  = str2double(get(handles.edit20,'string'));
            FramesZ = str2double(get(handles.edit22,'string'));
            ITI = str2double(get(handles.edit28,'string'));
            fold = str2double(get(handles.edit30,'string'));
            if fold ~= 0
                fold_fct = LinesY./fold;
            end
            
            
            plane_num = str2double(get(handles.edit18,'string'));
            plane_tot = str2double(get(handles.edit19,'string'));
            if ~isempty(get(handles.edit9,'string'))
                totTrials = str2double(get(handles.edit9,'string'));
            end
            tot_fr = (((FramesZ*plane_tot)+ITI)*totTrials);
            
            frames_sampled = plane_num:plane_tot:tot_fr;
            
            set(handles.edit23, 'String', baseline.selected_tresh);
            
            % cla(handles.axes4,'reset');
            hp = handles.axes4;
            hp.Interactions = [zoomInteraction regionZoomInteraction rulerPanInteraction];
            hp.Toolbar = [];
            colormap(hp, gray);
            
            % Screen('Preference','SkipSyncTests', 1);
            PsychDefaultSetup(2);
            % Seed the random number generator. Here we use the an older way to be
            % compatible with older systems. Newer syntax would be rng('shuffle'). Look
            % at the help function of rand "help rand" for more information
            rand('seed', sum(100 * clock));
            % Screen Number
            screenNumber = 2;%min(Screen('Screens'));
            % Define black, white and grey
            white = WhiteIndex(screenNumber);
            grey = white / 2;
            black = BlackIndex(screenNumber);
            % Open the screen
            [window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey, [], 32, 2,...
                [], [],  kPsychNeed32BPCFloat);
            % Flip to clear
            Screen('Flip', window);
            % Query the frame duration
            ifi = Screen('GetFlipInterval', window);
            % Maximum priority level
            topPriorityLevel = MaxPriority(window);
            % Get the centre coordinate of the window
            [xCenter, yCenter] = RectCenter(windowRect);
            %--------------------
            % Gabor information
            %--------------------
            % Dimensions
            % gaborDimPix = 500;
            gaborDimPixX = windowRect(3);
            gaborDimPixY = windowRect(4);
            % Sigma of Gaussian
            sigma = gaborDimPixX/8;
            % Obvious Parameters
            contrast = 1;
            aspectRatio = 1;
            % Spatial Frequency (Cycles Per Pixel)
            % One Cycle = Grey-Black-Grey-White-Grey i.e. One Black and One White Lobe
            if ~isempty(get(handles.edit7,'string'))
                numCycles = str2double(get(handles.edit7,'string'));
            else
                warning('Spatial Frequency automatically set to 5')
                numCycles = 5;
            end
            freq = numCycles / gaborDimPixX;
            % Build a procedural gabor texture
            gabortex = CreateProceduralGabor(window, gaborDimPixX, gaborDimPixY,...
                [], [0.5 0.5 0.5 1], 1, 10);
            % Positions of the Gabors
            dim = 8;
            [x, y] = meshgrid(-dim:dim, -dim:dim);
            % Calculate the distance in "Gabor numbers" of each gabor from the center
            % of the array
            dist = sqrt(x.^2 + y.^2);
            % Cut out an inner annulus
            innerDist = 3.5;
            x(dist <= innerDist) = nan;
            y(dist <= innerDist) = nan;
            % Cut out an outer annulus
            outerDist = 10;
            x(dist >= outerDist) = nan;
            y(dist >= outerDist) = nan;
            % Select only the finite values
            x = x(isfinite(x));
            y = y(isfinite(y));
            % Center the annulus coordinates in the centre of the screen
            xPos = x .* gaborDimPixX + xCenter;
            yPos = y .* gaborDimPixY + yCenter;
            % Count how many Gabors there are
            nGabors = numel(xPos);
            % Make the destination rectangles for all the Gabors in the array
            baseRect = [0 0 gaborDimPixX gaborDimPixY];
            allRects = nan(4, nGabors);
            for i = 1:nGabors
                allRects(:, i) = CenterRectOnPointd(baseRect, xPos(i), yPos(i));
            end
            % Drift speed for the 2D global motion
            % Randomise the Gabor orientations and determine the drift speeds of each gabor.
            % This is given by multiplying the global motion speed by the cosine
            % difference between the global motion direction and the global motion.
            % Here the global motion direction is 0. So it is just the cosine of the
            % angle we use. We re-orientate the array when drawing
            % Randomise the phase of the Gabors and make a properties matrix. We could
            % if we want have each Gabor with different properties in all dimensions.
            % Not just orientation and drift rate as we are doing here.
            % This is the power of using procedural textures
            phaseLine = rand(1, nGabors) .* 360;
            propertiesMat = repmat([NaN, freq, sigma, contrast, aspectRatio, 0, 0, 0],...
                nGabors, 1);
            propertiesMat(:, 1) = phaseLine';
            % Perform initial flip to gray background and sync us to the retrace:
            vbl = Screen('Flip', window);
            % Numer of frames to wait before re-drawing
            % angle to show
            % Animation loop
            % while ~KbCheck
            if ~isempty(get(handles.edit13,'string'))
                targetAngle = str2double(get(handles.edit13,'string'));
            else
                error('Select a target angle')
            end

            angle_bins  = mod((targetAngle - 90), 360) : 15 : targetAngle;
            act_bins  = nan(1, 7);
            act_bins(2:7)  = linspace(baseline.selected_range, baseline.selected_tresh, 6);
            act_bins(1) = act_bins(2) - (act_bins(3) - act_bins(2));
            act_bins = act_bins - baseline.selected_range;
            % xx = mod((targetAngle - 90), 360);
            
            cla(handles.axes1,'reset');
            h = line(nan, nan,'Parent',handles.axes1,'Color', 'k', 'LineWidth', 1);
            % h.Interactions = [zoomInteraction regionZoomInteraction rulerPanInteraction];
            % h.Toolbar = [];
            axis(handles.axes1,[-89,0,-5,90]);
            set(h,'XData', -89:0);
            
            cla(handles.axes3,'reset');
            h3 = line(nan(1,2), nan,'Parent',handles.axes3); %,'Color', 'k', 'LineWidth', 1);
            % h3.Interactions = [zoomInteraction regionZoomInteraction rulerPanInteraction];
            % h3.Toolbar = [];
            axis(handles.axes3,[-89,0,0,2]);
            set(h3(1:end),'XData', -89:0);
                        
            roiMask = baseline.roiMask;
            F0_moving = baseline.roiData(:,end-((FramesZ)*4)+1:end);
            nroi = baseline.nroi;
            loops = 0;
            ard = 0;
            success = 0;
            aversive = 0;
                        
            StimData.subpop1 = subpop1;
            StimData.subpop2 = subpop2;
            if ~isempty (subpop2)
                StimData.multiple_subpops = true;
                npops = 1;
                rng = baseline.selected_range;
                thresh = baseline.selected_tresh;
                N_above = baseline.selected_tresh - baseline.selected_range;
            else
                StimData.multiple_subpops = false;
                npops = 0;
                rng = 0;
                thresh = baseline.selected_tresh;
                N_above = baseline.selected_tresh;
            end
            StimData.('TStampGlobal') = nan(1, tot_fr);
            StimData.('Angle')  = nan(1, tot_fr/plane_tot);
            StimData.('targetAngle') = targetAngle;
            StimData.('z_score')  = nan(1, tot_fr/plane_tot);
            StimData.('Licks')  = ones(1, tot_fr/plane_tot);
            StimData.('Reward_frames')  = zeros(1, tot_fr/plane_tot);
            StimData.('Opto')  = zeros(1, tot_fr);
            StimData.('Aversive_frames')  = zeros(1, tot_fr/plane_tot);
            StimData.('bline_loaded')  = bline_loaded;
            StimData.('ImTranslateXY')  = nan(2, tot_fr/plane_tot);
            StimData.NewTrial = nan(1, 250);
            StimData.F0 = nan(nroi, tot_fr/plane_tot);
            StimData.Noise_exp = 0;
            StimData.Noise_df = 0;
            StimData.mapping = linspace(rng, thresh, 90);
            rotation = nan(1,tot_fr/plane_tot);
            rotation(1) = 0;
            shown = nan(1,tot_fr/plane_tot);
            shown(1) = 0;
            licks = nan(1,tot_fr/plane_tot);
            licks(1:plane_tot) = 0;
            idx = nan(1, tot_fr/plane_tot);
            x = nan(1, tot_fr/plane_tot);
            x(1:5) = 0;
            tt = 0;
            xpt = 0;
            ypt = 0;
            roi_one_frame = nan(nroi,tot_fr/plane_tot);
            n = 0;
            pl = 0;
            new_trial_fr = 0;
            rot_oneTrial = zeros(1, FramesZ);
            avv = 0;
            scc = 0;
            ITI_f = 90;
            
            
            ServerSend = tcpip('10.93.6.2',50000,'NetworkRole','server', ...
                'OutputBufferSize', 12, 'Terminator', 'CR/LF');
            fopen(ServerSend);
            fwrite(ServerSend, uint32([tot_fr plane_num plane_tot]), 'uint32');
            fclose(ServerSend);
            
            tcpipServer = tcpip('10.93.6.2',30000,'NetworkRole','client', 'Terminator', 'CR/LF');
            tcpipServer.InputBufferSize = (pixelsX*LinesY*200);
            fopen(tcpipServer);
            Screen('BlendFunction', window, 'GL_ONE', 'GL_ZERO');
            
            while pl < tot_fr
            
                if tcpipServer.BytesAvailable >= pixelsX * LinesY * 2
                    rawData = uint16(fread(tcpipServer,(pixelsX*LinesY),'uint16'));
                    rawData_tmp = uint16(reshape(rawData, pixelsX, LinesY));
                    rawData_tmp = permute(rawData_tmp, [2 1]);
                    rawData_tmp = imcomplement(rawData_tmp);

                    pl = pl + 1;
                else
                    % Not enough data; either skip this iteration or pause briefly.
                    continue;  % or use pause(0.001);
                end

                if pl == 1
                    timerVal = tic;
                end
                StimData.TStampGlobal(1, pl) = toc(timerVal);

                if sum(ismember(pl, frames_sampled)) == 0
                    continue
                end

                if fold ~= 0
                    r_tmp = uint16(zeros(fold, pixelsX, fold_fct));
                    for ifc = 1:fold_fct %#ok<BDSCI>
                        fuck_index = (fold*(ifc-1))+ 1 : fold*ifc;
                        r_tmp(:, :, ifc) = rawData_tmp(fuck_index, :, 1);
                    end
                    rawData_tmp = uint16(mean(r_tmp, 3));
                end
                n = n+1;
                new_trial_fr = new_trial_fr + 1;
            
                if new_trial_fr == 1
                    loops = loops + 1;
                    Prop_success = round((success/loops)*100);
                    newS = sprintf('%d', int32(Prop_success));
                    set(handles.edit25, 'String', newS );
                    StimData.NewTrial(loops) = n;
            
                    Prop_aversive = round((aversive/loops)*100);
                    newA = sprintf('%d', int32(Prop_aversive));
                    set(handles.edit26, 'String', newA );
            
                    newString = sprintf('%d', int32(loops));
                    set(handles.edit10, 'String', newString );
                end
            
                if mc == 1
                    output = nan(3, 4);
                    for i = 1:size(fft_tts, 3)
                        output(i,:) = dftregistration(fft_tts(:,:,i), fft2(double(rawData_tmp(vx1(i,1):vx1(i,2), vy1(i,1):vy1(i,2)))));
                    end
                    xsh = round(nanmedian(output(:,4)));
                    ysh = round(nanmedian(output(:,3)));
                    rawData_tmp = imtranslate(rawData_tmp, [xsh ysh]);
                    StimData.ImTranslateXY(:, n)  = [xsh; ysh];
            
                end
            
                if n == 1
                    hp = imshow(rawData_tmp, [], 'InitialMagnification', 800, 'Parent', hp, 'border', 'tight');
                end
            
                for i = 1:nroi
                    roi_one_frame(i,n) = double(mean(rawData_tmp(roiMask{i, 1})));
                end
                F0_moving(:,2:end) = F0_moving(:,1:end-1);
                F0_moving(:,1) = roi_one_frame(:,n);
                StimData.F0(:,n) = prctile(F0_moving, 10, 2);
            
                if npops == 0
                    x(n) = mean((roi_one_frame(:,n) - (StimData.F0))./(StimData.F0));
                    if x(n) < 0
                        x(n) = 0;
                    end
                elseif npops == 1
                    dfpop1 = mean((roi_one_frame(subpop1,n) - (StimData.F0(subpop1,n)))./(StimData.F0(subpop1,n)));
                    dfpop2 = mean((roi_one_frame(subpop2,n) - (StimData.F0(subpop2,n)))./(StimData.F0(subpop2,n)));
                    x(n) = (dfpop1 - dfpop2);
                end
            
                if n >= 3
                    z = mean(x(n-2:n)) - rng;
                else
                    z = x(n) - rng;
                end
            
                try
                    idx(n) = max(find(z >= act_bins)); %#ok<MXFND>
                catch
                    idx(n) = 1;
                end
            
                xx = angle_bins(idx(n));
            
                if new_trial_fr <= 15
                    rotation(n) = round((z/N_above)*90);
                    shown(n) = -5;
                    rot_oneTrial(new_trial_fr) = 0;
                    StimData.Angle(n)  = -5;
                    if new_trial_fr == 1
                        Screen('FillRect', window, [0.5 0.5 0.5]);
                        Screen('Flip', window, 0, 1);
                    end
                elseif new_trial_fr > 15 && new_trial_fr <= FramesZ
                    if sum(rot_oneTrial >= 90) == 0 && sum(rot_oneTrial < 0) == 0
                        ypt = ypt + 1;
                        rotation(n) = round((z/N_above)*90);
                        shown(n) = rotation(n);
                        rot_oneTrial(new_trial_fr) = rotation(n);
                        StimData.Angle(n)  = xx;
                        if idx(n) ~= idx(n -1) || ypt == 1
                            Screen('DrawTexture', window, gabortex, [], [], xx,...
                                [], [], [], [], kPsychDontDoRotation, propertiesMat');
                            Screen('Flip', window, 0, 1);
                        end
                    elseif sum(rot_oneTrial >= 90) > 0
                        tt = tt+1;
                        if tt < 60
                            rotation(n) = round((z/N_above)*90);
                            shown(n) = targetAngle;
                            rot_oneTrial(new_trial_fr) = 90;
                            StimData.Angle(n)  = targetAngle;
                            if idx(n) ~= idx(n -1)
                                Screen('DrawTexture', window, gabortex, [], [], targetAngle,...
                                    [], [], [], [], kPsychDontDoRotation, propertiesMat');
                                Screen('Flip', window, 0, 1);
                            end
                        elseif tt == 60
                            tt = 0;
                            ypt = 0;
                            new_trial_fr  = FramesZ;
                            rotation(n) = round((z/N_above)*90);
                            shown(n) = -1;
                            StimData.Angle(n)  = -1;
                        end
                    elseif sum(rot_oneTrial < 0) > 0
                        tt = tt+1;
                        if tt < 60
                            rotation(n) = round((z/N_above)*90);
                            shown(n) = 0;
                            rot_oneTrial(new_trial_fr) = 0;
                            StimData.Angle(n)  = 0;
                            if idx(n) ~= idx(n -1)
                                Screen('DrawTexture', window, gabortex, [], [], 0,...
                                    [], [], [], [], kPsychDontDoRotation, propertiesMat');
                                Screen('Flip', window, 0, 1);
                            end
                            if tt == 1
                                ITI_f = 15;
                            end
                        elseif tt == 60
                            tt = 0;
                            ypt = 0;
                            new_trial_fr  = FramesZ;
                            rotation(n) = round((z/N_above)*90);
                            shown(n) = -1;
                            StimData.Angle(n)  = -1;
                        end
                    end
                elseif new_trial_fr > FramesZ && new_trial_fr < FramesZ + (ITI_f)
                    if new_trial_fr == FramesZ + 1 && ITI_f == 90
                        Screen('FillRect', window, [0 0 0]);
                        Screen('Flip', window, 0, 1);
                    elseif new_trial_fr == FramesZ + 1 && ITI_f == 15
                        Screen('FillRect', window, [0.5 0.5 0.5]);
                        Screen('Flip', window, 0, 1);
                    end
                    rotation(n) = round((z/N_above)*90);
                    shown(n) = -1;
                    StimData.Angle(n)  = -1;
                    if new_trial_fr == FramesZ + 30 && ITI_f == 90
                        white_noise(2);
                        StimData.Aversive_frames(n) = 1;
                    end
                elseif new_trial_fr == FramesZ + (ITI_f)
                    tt = 0;
                    ypt = 0;
                    xpt = 0;
                    new_trial_fr = 0;
                    rot_oneTrial = zeros(1, FramesZ);
                    rotation(n) = round((z/N_above)*90);
                    shown(n) = -1;
                    StimData.Angle(n)  = -1;
                    avv = 0; scc = 0; ITI_f = 90;
                end
            
                set(hp, 'CData',im2uint16(rawData_tmp));
            
                if n < 90
                    pltRot = fliplr(shown(1:n));
                    pltRot(end+1 : 90) = nan;
                    pltRot = fliplr(pltRot);
            
                    pltLick = fliplr(licks(1:n));
                    pltLick(end+1 : 90) = nan;
                    pltLick = fliplr(pltLick);
            
                    pltRew = fliplr(StimData.Reward_frames(1:n));
                    pltRew(end+1 : 90) = nan;
                    pltRew = fliplr(pltRew);
                else
                    pltRot  = shown(n-89 : n);
                    pltLick = licks(n-89 : n);
                    pltRew = StimData.Reward_frames(n-89 : n);
                end
                set(h,'YData',pltRot);
                set(h3(1), 'YData', pltLick);
                set(h3(2), 'YData', pltRew);
            
                try
                    drawnow limitrate nocallbacks
                catch
                    warning('Could not plot')
                end
            end
            
            StimData.Rotation = rotation;
            StimData.RoiPerFrame  = roi_one_frame;
            StimData.z_score  = x;
            StimData.numCycles = numCycles;
            StimData.threshold_used = {rng; thresh};
            StimData.b_line_all = baseline.F0;
            StimData.plane_num = plane_num;
            StimData.plane_tot = plane_tot;
            StimData.RefImage = tts;
            
            if saveTrials == 1
                Folder = '\BCI_Trials';
                filename = ['\' num2str(1)];
                if ~exist([directory Folder], 'dir')
                    mkdir([directory Folder])
                end
                save([directory Folder filename], 'StimData', '-v6')
            end
            fprintf('That is it. Saving data')
            
            % Clean up
            fclose(tcpipServer);
            delete(tcpipServer);
            delete(ServerSend);
            ServerSend = [];
            clear all
            sca;
        end

        % Callback function
        function pushbutton5_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to pushbutton5 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            saveTrials = get(handles.checkbox6,'Value');
            if saveTrials == 0
                asw = questdlg('You are set not to save this trial. Are you sure this is right?');
                switch asw
                    case 'Yes'
                        disp ('As you command. I am only here to serve to my true master')
                    case 'No'
                        return
                    case 'cancel'
                        return
                end
            end
            directory = get(handles.edit11, 'String');
            if saveTrials == 1
                if strcmp('D:\Experiments\SessionName',directory)
                    answer = questdlg('You are using the default directory to save your data. Are you sure this is right?');
                    switch answer
                        case 'Yes'
                            disp ('As you command. I am only here to serve to my true master')
                        case 'No'
                            return
                        case 'cancel'
                            return
                    end
                end
            end
            % Here we call some default settings for setting up Psychtoolbox
            PsychDefaultSetup(2);
            
            % Seed the random number generator. Here we use the an older way to be
            % compatible with older systems. Newer syntax would be rng('shuffle'). Look
            % at the help function of rand "help rand" for more information
            rand('seed', sum(100 * clock));
            
            % Screen Number
            screenNumber = 2;%min(Screen('Screens'));
            
            % Define black, white and grey
            white = WhiteIndex(screenNumber);
            grey = white / 2;
            black = BlackIndex(screenNumber);
            
            % Open the screen
            [window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey, [], 32, 2,...
                [], [],  kPsychNeed32BPCFloat);
            
            % Flip to clear
            Screen('Flip', window);
            
            % Query the frame duration
            ifi = Screen('GetFlipInterval', window);
            
            % Maximum priority level
            topPriorityLevel = MaxPriority(window);
            
            % Get the centre coordinate of the window
            [xCenter, yCenter] = RectCenter(windowRect);
            
            %--------------------
            % Gabor information
            %--------------------
            
            % Dimensions
            % gaborDimPix = 500;
            gaborDimPixX = windowRect(3);
            gaborDimPixY = windowRect(4);
            
            % Sigma of Gaussian
            sigma = gaborDimPixX ;
            
            % Obvious Parameters
            orientation = 90;
            contrast = 1;
            aspectRatio = 5;
            
            % Spatial Frequency (Cycles Per Pixel)
            % One Cycle = Grey-Black-Grey-White-Grey i.e. One Black and One White Lobe
            if ~isempty(get(handles.edit5,'string'))
                numCycles = str2double(get(handles.edit5,'string'));
            else
                warning('Spatial Frequency automatically set to 5')
                numCycles = 5;
            end
            freq = numCycles / gaborDimPixX;
            
            % Build a procedural gabor texture
            gabortex = CreateProceduralGabor(window, gaborDimPixX, gaborDimPixY,...
                [], [1 1 1 0.5], 1, 100);
            
            % Positions of the Gabors
            dim = 8;
            [x, y] = meshgrid(-dim:dim, -dim:dim);
            
            % Calculate the distance in "Gabor numbers" of each gabor from the center
            % of the array
            dist = sqrt(x.^2 + y.^2);
            
            % Cut out an inner annulus
            innerDist = 3.5;
            x(dist <= innerDist) = nan;
            y(dist <= innerDist) = nan;
            
            % Cut out an outer annulus
            outerDist = 10;
            x(dist >= outerDist) = nan;
            y(dist >= outerDist) = nan;
            
            % Select only the finite values
            x = x(isfinite(x));
            y = y(isfinite(y));
            
            % Center the annulus coordinates in the centre of the screen
            xPos = x .* gaborDimPixX + xCenter;
            yPos = y .* gaborDimPixY + yCenter;
            
            % Count how many Gabors there are
            nGabors = numel(xPos);
            
            % Make the destination rectangles for all the Gabors in the array
            baseRect = [0 0 gaborDimPixX gaborDimPixY];
            allRects = nan(4, nGabors);
            for i = 1:nGabors
                allRects(:, i) = CenterRectOnPointd(baseRect, xPos(i), yPos(i));
            end
            
            % Drift speed for the 2D global motion
            if ~isempty(get(handles.edit6,'string'))
                tFreq = str2double(get(handles.edit6,'string'));
            else
                warning('Temporal Frequency automatically set to 2')
                tFreq = 2;
            end
            degPerSec = 360 * tFreq;
            degPerFrame =  degPerSec * ifi;
            
            % Randomise the Gabor orientations and determine the drift speeds of each gabor.
            % This is given by multiplying the global motion speed by the cosine
            % difference between the global motion direction and the global motion.
            % Here the global motion direction is 0. So it is just the cosine of the
            % angle we use. We re-orientate the array when drawing
            gaborAngles = 0; % rand(1, nGabors) .* 180 - 90;
            degPerFrameGabors = cosd(gaborAngles) .* degPerFrame;
            
            % Randomise the phase of the Gabors and make a properties matrix. We could
            % if we want have each Gabor with different properties in all dimensions.
            % Not just orientation and drift rate as we are doing here.
            % This is the power of using procedural textures
            phaseLine = 0.5;% rand(1, nGabors) .* 360;
            propertiesMat = repmat([NaN, freq, sigma, contrast, aspectRatio, 0, 0, 0],...
                nGabors, 1);
            propertiesMat(:, 1) = phaseLine';
            
            vbl = Screen('Flip', window);
            
            waitframes = 1;
            
            AngleNumb = str2double(get(handles.uibuttongroup1.SelectedObject,'String'));
            if AngleNumb == 8
                angles = 0:45:315;
            elseif AngleNumb == 12
                angles = 0:30:330;
            end
            
            
            
            if ~isempty(get(handles.edit2,'string'))
                tgrey = str2double(get(handles.edit2,'string'));
            else
                tgrey = 0;
            end
            if ~isempty(get(handles.edit3,'string'))
                tblack = str2double(get(handles.edit3,'string'));
            else
                tblack = 0;
            end
            if ~isempty(get(handles.edit4,'string'))
                tdrift = str2double(get(handles.edit4,'string'));
            else
                tdrift = 0;
            end
            
            
            presTimeDrift = repmat(tdrift, 1, numel(angles));
            % Animation loop
            % while ~KbCheck
            greyShow = get(handles.checkbox3,'Value');
            blackShow = get(handles.checkbox4,'Value');
            randomShow = get(handles.checkbox2,'Value');
            
            
            
            if randomShow == 1
                angles = angles(randperm(length(angles)));
            end
            
            if greyShow == 1
                GreyStim = repmat(5, 1, numel(angles)+1);
                presTimeGrey = repmat(tgrey, 1, numel(angles)+1);
            
                vStim = nan(1,length(angles) + length(GreyStim));
                vStim(1:2:end) = GreyStim;
                vStim(2:2:end) = angles;
            
                vStimTimes = nan(1,length(presTimeDrift) + length(presTimeGrey));
                vStimTimes(1:2:end) = presTimeGrey;
                vStimTimes(2:2:end) = presTimeDrift;
            else
                vStim = angles;
                vStimTimes = presTimeDrift;
            end
            
            BlTime = tblack;
            if blackShow == 1
                vStim = [3 vStim 3];
                vStimTimes = [BlTime vStimTimes BlTime];
            end
            
            trig = get(handles.checkbox9,'Value');
            if trig == 1
                ServerSend = tcpip('10.93.6.2',50000,'NetworkRole','server', ...
                    'OutputBufferSize', 12, 'Terminator', 'CR/LF');
                fopen(ServerSend);
                fwrite(ServerSend, uint32([0 1 1]), 'uint32');
                fclose(ServerSend);
                delete(ServerSend);
                ServerSend = [];
            end
            
            tduration = sum(vStimTimes,2);
            StimData(1).('TStamp') = nan(1, round(length(tduration/ifi)*1.5));
            StimData(1).('Angle')  = nan(1, round(length(tduration/ifi)*1.5));
            n = 0;
            timerGlobal = tic;
            
            repetions = 100;
            rep = 0;
            while rep < repetions
                rep = rep+1;
                for i = 1:numel(vStim)
                    timerVal = tic;
                    if vStim(i) ~= 3 && vStim(i) ~= 5
                        while toc(timerVal) < vStimTimes(i)
                            % Set the right blend function for drawing the gabors
                            Screen('BlendFunction', window, 'GL_ONE', 'GL_ZERO');
                            % Draw texture
                            Screen('DrawTexture', window, gabortex, [], [], vStim(i),...
                                [], [], [], [], kPsychDontDoRotation, propertiesMat');
                            % Flip our drawing to the screen
                            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
            
                            % Increment the phase of our Gabors
                            phaseLine = phaseLine + degPerFrameGabors;
                            propertiesMat(:, 1) = phaseLine';
            
                            n = n+1;
                            StimData.TStampGlobal(n) = toc(timerGlobal);
                            StimData.TStampLocal(n) = toc(timerVal);
                            StimData.Angle(n)  = vStim(i);
                        end
                    elseif vStim(i) == 5
                        while toc(timerVal) < vStimTimes(i)
                            % Color the screen a random color
                            Screen('FillRect', window, [0.5 0.5 0.5]);
            
                            % Flip to the screen
                            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
            
                            n = n+1;
                            StimData.TStampGlobal(n) = toc(timerGlobal);
                            StimData.TStampLocal(n) = toc(timerVal);
                            StimData.Angle(n)  = vStim(i);
                        end
                    elseif vStim(i) == 3
                        while toc(timerVal) < vStimTimes(i)
                            Screen('FillRect', window, [0 0 0]);
            
                            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
            
                            n = n+1;
                            StimData.TStampGlobal(n) = toc(timerGlobal);
                            StimData.TStampLocal(n) = toc(timerVal);
                            StimData.Angle(n)  = vStim(i);
                        end
                    end
                end
            end
            Fld = get(handles.uibuttongroup2.SelectedObject,'String');
            if saveTrials == 1
                currentCounterValue = str2double(get(handles.edit10, 'String'));
                newString = sprintf('%d', int32(currentCounterValue +1));
                set(handles.edit10, 'String', newString );
                trialNumber = currentCounterValue;
            
                if strcmp(Fld,'Initial')
                    Folder = '\Gratings_Start_Trials';
                elseif strcmp(Fld,'Final')
                    Folder = '\Gratings_End_Trials';
                end
                filename = ['\' num2str(trialNumber)];
                if ~exist([directory Folder], 'dir')
                    mkdir([directory Folder])
                end
                save([directory Folder filename], 'StimData')
            end
            sca;
        end

        % Button pushed function: pushbutton6
        function pushbutton6_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to pushbutton6 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            global tts vx vy Planes_imgs
            
            mc = get(handles.checkbox12,'Value');
            h = handles.axes4;
            % hold( h  , 'off')
            AllignOnly = get(handles.checkbox13,'Value');
            
            if ~isempty(get(handles.edit17,'string'))
                num = str2double(get(handles.edit17,'string'));
            else
                num = -1;
                warning('You did not selet how many frame you want to load. Defaults to all frames')
            end
            handles.num = num;
            
            plane_num = str2num(get(handles.edit18,'string'));
            plane_tot = str2double(get(handles.edit19,'string'));
            
            answer = questdlg('How do you want to load your image', ...
                'Options', ...
                'Acquire new one','Load from directory','Acquire new one');
            
            switch answer
                case 'Load from directory'
                    if ~isempty(get(handles.edit14,'string'))
                        fname = get(handles.edit14,'string');
                    else
                        error('Select a file on which you want to draw your ROIs')
                    end
                    rawData = sbxgrabframe(fname,ind,num);
                    rawData = squeeze(double(rawData));
                    rawData = rawData(:,:, plane_num:plane_tot:size(rawData,3));
                case 'Acquire new one'
                    pixelsX = str2double(get(handles.edit21,'string'));
                    LinesY  = str2double(get(handles.edit20,'string'));
                    FramesZ = str2double(get(handles.edit17,'string'));
                    fold = str2double(get(handles.edit30,'string'));
                    if fold ~= 0
                        fold_fct = LinesY./fold;
                    end

                    frames_sampled_tmp = cell(1, numel(plane_num));
                    for iplane = 1:numel(frames_sampled_tmp)
                        frames_sampled_tmp{iplane} = plane_num(iplane):plane_tot:(FramesZ*plane_tot);
                    end
                    C = frames_sampled_tmp;
                    maxLen = max(cellfun(@length, C));
                    padded = cellfun(@(v) [v, nan(1, maxLen - numel(v))], C, 'UniformOutput', false);
                    frames_sampled = vertcat(padded{:});

                    rawData_t = nan((pixelsX*LinesY),FramesZ);
                    plane_of_int= nan(1, FramesZ);
                    if fold ~= 0
                        rawData = uint16(zeros(fold, pixelsX,  FramesZ));
                    else
                        rawData = uint16(zeros(LinesY, pixelsX, FramesZ));
                    end

                    ServerSend = tcpip('10.93.6.2',50000,'NetworkRole','server', ...
                        'OutputBufferSize', 12, 'Terminator', 'CR/LF');
                    fopen(ServerSend);
                    fwrite(ServerSend, uint32([(FramesZ*plane_tot)  1  plane_tot]), 'uint32');
                    fclose(ServerSend);
            

                    tcpipServer = tcpip('10.93.6.2',30000,'NetworkRole','client',...
                        'Terminator', 'CR/LF');
                    tcpipServer.InputBufferSize = (pixelsX*LinesY*200);
                    fopen(tcpipServer);
            
                    pl = 0;
                    n = 0;
                    while pl < FramesZ*plane_tot
            
                        pl = pl + 1;
            
                        if sum(ismember(pl, frames_sampled)) == 0
                            fread(tcpipServer,(pixelsX*LinesY),'uint16');
                            continue
                        else
                            n = n+1;
                            rawData_t(:,n) = fread(tcpipServer,(pixelsX*LinesY),'uint16');
                            [plane_of_int(1, n), ~] = find(frames_sampled == pl, 1);
                        end
            
                        rawData_tmp = uint16(reshape(rawData_t(:,n), pixelsX, LinesY));
                        rawData_tmp = permute(rawData_tmp, [2 1]);
                        rawData_tmp = imcomplement(rawData_tmp);
            
                        if fold ~= 0
                            r_tmp = uint16(zeros(fold, pixelsX, fold_fct));
                            for ifc = 1:fold_fct %#ok<BDSCI>
                                fuck_index = (fold*(ifc-1))+ 1 : fold*ifc;
                                r_tmp(:, :, ifc) = rawData_tmp(fuck_index, :, 1);
                            end
                            rawData_tmp = uint16(mean(r_tmp, 3));
                        end
                        rawData(:, :, n) = rawData_tmp;
            
                        if n == 1
                            h = imshow(rawData_tmp, [], 'InitialMagnification', 800, 'Parent', h, 'border', 'tight');
                        end
                        if plane_of_int(1, n) == 1
                            set(h, 'CData', rawData_tmp);
                            drawnow limitrate nocallbacks
                        end
                    end
                    fclose(tcpipServer);
                    delete(tcpipServer);
                    delete(ServerSend);
                    ServerSend = [];
            
                    if num == -1
                        num = FramesZ;
                    end
            end
            
            n_mc = 3;
            n_mc_fld = 1;
            
            nfr = get(handles.uibuttongroup8.SelectedObject,'String');
            Planes_imgs = uint16(zeros(size(rawData, 1), size(rawData, 2), numel(plane_num)));
            if strcmp (nfr, 'Mean')
                for ipla = 1:numel(plane_num)
                    Planes_imgs(:, :, ipla) = uint16(mean(rawData(:, :, plane_of_int == ipla), 3));
                end

                tts = uint16(mean(rawData(:, :, plane_of_int == 1), 3));
                h = handles.axes4;
                himage = imshow((tts), [], 'InitialMagnification', 800, 'Parent', h, 'border', 'tight');
                set(handles.slider1, 'Min', 1, 'Max',max(plane_of_int), 'Value', 1, ...
                        'SliderStep', [0 0]);
                h.Toolbar.Visible = 'on';
                if mc == 1
                    tts = uint16(tts);
            
                    if AllignOnly == 0
                        answer = questdlg('Use this image for Motion Correction?', ...
                            'Options', ...
                            'Yes','No', 'Yes');
                        switch answer
                            case 'Yes'
                                if fold == 0
                                    [xcoord, ycoord] = ginput(n_mc);
                                    vx = nan(n_mc, 1);
                                    vy = nan(n_mc, 1);
                                    for ixx = 1:size(xcoord)
                                        vx(ixx, [1, 2]) = [round(xcoord(ixx) - 99), round(xcoord(ixx) + 99)];
                                        vy(ixx, [1, 2]) = [round(ycoord(ixx) - 64), round(ycoord(ixx) + 63)];
                                    end
            
                                    if any(vx(:, 1) < 100)
                                        idx = find(vx(:, 1) < 100);
                                        vx(idx, 1) = 100; vx(idx, 2) = 298;
                                    end
                                    if any(vx(:, 2) > 794)
                                        idx = find(vx(:, 2) > 794);
                                        vx(idx, 1) = 596; vx(idx, 2) = 794;
                                    end
                                    if any(vy(:, 1) < 1)
                                        idx = find(vy(:, 1) < 1);
                                        vy(idx, 1) = 1; vy(idx, 2) = 128;
                                    end
                                    if any(vy(:, 2) > 510)
                                        idx = find(vy(:, 2) > 510);
                                        vy(idx, 1) = 383; vy(idx, 2) = 510;
                                    end
                                else
                                    [xcoord, ~] = ginput(n_mc);
                                    vx = nan(n_mc, 1);
                                    vy = nan(n_mc, 1);
            
                                    for ixx = 1:size(xcoord)
                                        vx(ixx, [1, 2]) = [178, 535];
                                        vy(ixx, [1, 2]) = [1, 60];
                                    end
                                end
                            case 'No'
                                vx = nan;
                                vy = nan;
                        end
                    elseif AllignOnly == 1
                        vx = nan;
                        vy = nan;
                    end
                end
            
            elseif strcmp (nfr, 'Max')
                tts = uint16(max(rawData,[], 3));
                h = handles.axes4;
                himage = imshow(tts, [], 'InitialMagnification', 800, 'Parent', h, 'border', 'tight');
                set(handles.slider1, 'Min', 1, 'Max', 1, 'Value', 1, ...
                    'SliderStep', [0 0]);
            %     h.Toolbar.Visible = 'on';
                if mc == 1
                    tts = uint16(tts);
                    vy = size(tts,2)/4;
                    vx = size(tts,1)/4;
                end
            
            
            elseif strcmp (nfr, 'Std')
                tts = uint16(std(rawData,0,3));
                h = handles.axes4;
                himage = imshow(tts, [], 'InitialMagnification', 800, 'Parent', h, 'border', 'tight');
                set(handles.slider1, 'Min', 1, 'Max', 1, 'Value', 1, ...
                        'SliderStep', [0 0]);
            %     h.Toolbar.Visible = 'on';
                if mc == 1
                    tts = uint16(tts);
                    vy = size(tts,2)/4;
                    vx = size(tts,1)/4;
                end
            
            
            elseif strcmp (nfr, 'All frames')
                tts= rawData;
                h = handles.axes4;
                himage = imshow(tts(:,:,1), [], 'InitialMagnification', 800, 'Parent', h, 'border', 'tight');
                set(himage, 'CData',tts(:,:,1));
                h.Toolbar.Visible = 'on';
                if num == -1
                    set(handles.slider1, 'Min', 1, 'Max', size(tts,3), ...
                        'SliderStep', [1 1]/(size(tts,3)-1), 'Value', 1);
                else
                    set(handles.slider1, 'Min', 1, 'Max', num, ...
                        'SliderStep', [1 1]/(num-1),  'Value', 1);
                end
            end
            hold(handles.axes4, 'on')
            handles.frameSlider = tts;
            handles.nR = 0;
            handles.stopROI = 0;
            if isfield(handles, 'ROI')
                delete(handles.ROI);
            end
            guidata(hObject, handles);
        end

        % Value changed function: slider1
        function slider1_Callback(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to slider1 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: get(hObject,'Value') returns position of slider
            %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
            % get(handles.slider1,'Value')
            %
            % set(handles.slider1, 'Min', 1, 'Max', 1, ...
            %         'SliderStep', [1 1]/1 - 1, 'Value', 1);
            
            global Planes_imgs

            sl_v = get(hObject, 'Value');
            h = handles.axes4;
            imshow(Planes_imgs(:, :, round(sl_v)), [], 'InitialMagnification', 800, 'Parent', h, 'border', 'tight');
            h.Toolbar.Visible = 'on';

        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create figure1 and hide until all components are created
            app.figure1 = uifigure('Visible', 'off');
            app.figure1.Position = [761 155 1449 994];
            app.figure1.Name = 'Run_GUI';
            app.figure1.HandleVisibility = 'callback';
            app.figure1.Tag = 'figure1';

            % Create axes4
            app.axes4 = uiaxes(app.figure1);
            app.axes4.FontSize = 13.3333333333333;
            app.axes4.NextPlot = 'replace';
            app.axes4.Tag = 'axes4';
            app.axes4.Position = [396 395 1025 591];

            % Create uipanel2
            app.uipanel2 = uipanel(app.figure1);
            app.uipanel2.Title = 'BCI Panel';
            app.uipanel2.Tag = 'uipanel2';
            app.uipanel2.FontSize = 18.6666666666667;
            app.uipanel2.Position = [31 528 354 302];

            % Create axes1
            app.axes1 = uiaxes(app.uipanel2);
            app.axes1.FontSize = 11.3333333333333;
            app.axes1.NextPlot = 'replace';
            app.axes1.Tag = 'axes1';
            app.axes1.Position = [10 115 337 133];

            % Create pushbutton2
            app.pushbutton2 = uibutton(app.uipanel2, 'push');
            app.pushbutton2.ButtonPushedFcn = createCallbackFcn(app, @pushbutton2_Callback, true);
            app.pushbutton2.Tag = 'pushbutton2';
            app.pushbutton2.FontSize = 13.3333333333333;
            app.pushbutton2.Position = [230 5 80.8615384615384 30.8027027027027];
            app.pushbutton2.Text = 'Run BCI';

            % Create text2
            app.text2 = uilabel(app.uipanel2);
            app.text2.Tag = 'text2';
            app.text2.HorizontalAlignment = 'center';
            app.text2.VerticalAlignment = 'top';
            app.text2.FontSize = 10.6666666666667;
            app.text2.Position = [120 95 97 22];
            app.text2.Text = 'Feedback readout';

            % Create edit9
            app.edit9 = uieditfield(app.uipanel2, 'text');
            app.edit9.Tag = 'edit9';
            app.edit9.HorizontalAlignment = 'center';
            app.edit9.FontSize = 10.6666666666667;
            app.edit9.Position = [81 72 30.3230769230769 21.2432432432432];
            app.edit9.Value = '50';

            % Create text10
            app.text10 = uilabel(app.uipanel2);
            app.text10.Tag = 'text10';
            app.text10.HorizontalAlignment = 'center';
            app.text10.VerticalAlignment = 'top';
            app.text10.FontSize = 10.6666666666667;
            app.text10.Position = [25 73 49.4153846153846 15.9324324324324];
            app.text10.Text = 'Trials';

            % Create edit22
            app.edit22 = uieditfield(app.uipanel2, 'text');
            app.edit22.Tag = 'edit22';
            app.edit22.HorizontalAlignment = 'center';
            app.edit22.FontSize = 10.6666666666667;
            app.edit22.Tooltip = 'Target Activation';
            app.edit22.Position = [81 49 32 20];
            app.edit22.Value = '870';

            % Create text22
            app.text22 = uilabel(app.uipanel2);
            app.text22.Tag = 'text22';
            app.text22.HorizontalAlignment = 'center';
            app.text22.VerticalAlignment = 'top';
            app.text22.FontSize = 10.6666666666667;
            app.text22.Tooltip = 'Target Activation';
            app.text22.Position = [9 49 70.7538461538461 19.1189189189189];
            app.text22.Text = 'Frames (z)';

            % Create pushbutton14
            app.pushbutton14 = uibutton(app.uipanel2, 'push');
            app.pushbutton14.ButtonPushedFcn = createCallbackFcn(app, @pushbutton14_Callback, true);
            app.pushbutton14.Tag = 'pushbutton14';
            app.pushbutton14.FontSize = 13.3333333333333;
            app.pushbutton14.Position = [24 5 82 31];
            app.pushbutton14.Text = 'Test';

            % Create pushbutton15
            app.pushbutton15 = uibutton(app.uipanel2, 'push');
            app.pushbutton15.ButtonPushedFcn = createCallbackFcn(app, @pushbutton15_Callback, true);
            app.pushbutton15.Tag = 'pushbutton15';
            app.pushbutton15.FontSize = 13.3333333333333;
            app.pushbutton15.Position = [132 5 73 30.8027027027027];
            app.pushbutton15.Text = 'Baseline';

            % Create edit23
            app.edit23 = uieditfield(app.uipanel2, 'text');
            app.edit23.Tag = 'edit23';
            app.edit23.HorizontalAlignment = 'center';
            app.edit23.FontSize = 10.6666666666667;
            app.edit23.Position = [174 70 39.3076923076923 22.3054054054054];
            app.edit23.Value = '3';

            % Create text28
            app.text28 = uilabel(app.uipanel2);
            app.text28.Tag = 'text28';
            app.text28.HorizontalAlignment = 'center';
            app.text28.VerticalAlignment = 'top';
            app.text28.FontSize = 10.6666666666667;
            app.text28.Position = [121 71 50.5384615384615 19.1189189189189];
            app.text28.Text = 'Target';

            % Create checkbox11
            app.checkbox11 = uicheckbox(app.uipanel2);
            app.checkbox11.Tag = 'checkbox11';
            app.checkbox11.Text = 'Stim Baseline';
            app.checkbox11.FontSize = 10.6666666666667;
            app.checkbox11.Position = [225 250 97.6 24.8];

            % Create text29
            app.text29 = uilabel(app.uipanel2);
            app.text29.Tag = 'text29';
            app.text29.HorizontalAlignment = 'center';
            app.text29.VerticalAlignment = 'top';
            app.text29.FontSize = 10.6666666666667;
            app.text29.Position = [128 51 31.4461538461538 15.9324324324324];
            app.text29.Text = 'I.T.I.';

            % Create edit28
            app.edit28 = uieditfield(app.uipanel2, 'text');
            app.edit28.Tag = 'edit28';
            app.edit28.HorizontalAlignment = 'center';
            app.edit28.FontSize = 10.6666666666667;
            app.edit28.Position = [174 49 39.3076923076923 23.3675675675676];
            app.edit28.Value = '150';

            % Create checkbox12
            app.checkbox12 = uicheckbox(app.uipanel2);
            app.checkbox12.Tag = 'checkbox12';
            app.checkbox12.Text = {'Motion Correct'; ''};
            app.checkbox12.FontSize = 10.6666666666667;
            app.checkbox12.Position = [39 250 124.32 24.4297297297297];
            app.checkbox12.Value = true;

            % Create text31
            app.text31 = uilabel(app.uipanel2);
            app.text31.Tag = 'text31';
            app.text31.HorizontalAlignment = 'center';
            app.text31.VerticalAlignment = 'top';
            app.text31.FontSize = 10.6666666666667;
            app.text31.Position = [261 50 31.2 16];
            app.text31.Text = 'Fold';

            % Create edit30
            app.edit30 = uieditfield(app.uipanel2, 'text');
            app.edit30.Tag = 'edit30';
            app.edit30.HorizontalAlignment = 'center';
            app.edit30.FontSize = 10.6666666666667;
            app.edit30.Position = [297 48 39.2 23.2];
            app.edit30.Value = '0';

            % Create edit27
            app.edit27 = uieditfield(app.uipanel2, 'text');
            app.edit27.Tag = 'edit27';
            app.edit27.HorizontalAlignment = 'center';
            app.edit27.FontSize = 10.6666666666667;
            app.edit27.Position = [298 73 40 21];
            app.edit27.Value = '15000';

            % Create text27
            app.text27 = uilabel(app.uipanel2);
            app.text27.Tag = 'text27';
            app.text27.HorizontalAlignment = 'center';
            app.text27.VerticalAlignment = 'top';
            app.text27.FontSize = 10.6666666666667;
            app.text27.Position = [234 67 62 22];
            app.text27.Text = 'Baseline (z)';

            % Create uipanel4
            app.uipanel4 = uipanel(app.figure1);
            app.uipanel4.Title = 'General';
            app.uipanel4.Tag = 'uipanel4';
            app.uipanel4.FontSize = 18.6666666666667;
            app.uipanel4.Position = [38 847 353 138];

            % Create edit10
            app.edit10 = uieditfield(app.uipanel4, 'text');
            app.edit10.Tag = 'edit10';
            app.edit10.HorizontalAlignment = 'center';
            app.edit10.FontSize = 13.3333333333333;
            app.edit10.Position = [109 73 42.6186495176849 24.1009009009009];
            app.edit10.Value = '1';

            % Create text11
            app.text11 = uilabel(app.uipanel4);
            app.text11.Tag = 'text11';
            app.text11.HorizontalAlignment = 'center';
            app.text11.VerticalAlignment = 'top';
            app.text11.FontSize = 13.3333333333333;
            app.text11.Position = [9 73 97.5742765273312 19.7189189189189];
            app.text11.Text = 'Trial Number';

            % Create checkbox6
            app.checkbox6 = uicheckbox(app.uipanel4);
            app.checkbox6.Tag = 'checkbox6';
            app.checkbox6.Text = {'Save'; ''};
            app.checkbox6.FontSize = 13.3333333333333;
            app.checkbox6.Position = [18 42 61.6848874598071 29.5783783783784];
            app.checkbox6.Value = true;

            % Create edit11
            app.edit11 = uieditfield(app.uipanel4, 'text');
            app.edit11.Tag = 'edit11';
            app.edit11.HorizontalAlignment = 'center';
            app.edit11.FontSize = 13.3333333333333;
            app.edit11.Position = [28 4 299 25];
            app.edit11.Value = 'C:\Valerio\Testing_BCI';

            % Create text12
            app.text12 = uilabel(app.uipanel4);
            app.text12.Tag = 'text12';
            app.text12.HorizontalAlignment = 'center';
            app.text12.VerticalAlignment = 'top';
            app.text12.FontSize = 13.3333333333333;
            app.text12.Position = [130 44 61.6848874598071 19.7189189189189];
            app.text12.Text = 'Path';

            % Create edit25
            app.edit25 = uieditfield(app.uipanel4, 'text');
            app.edit25.Tag = 'edit25';
            app.edit25.HorizontalAlignment = 'center';
            app.edit25.FontSize = 13.3333333333333;
            app.edit25.Position = [277 85 42.6186495176849 24.1009009009009];
            app.edit25.Value = '0';

            % Create text25
            app.text25 = uilabel(app.uipanel4);
            app.text25.Tag = 'text25';
            app.text25.HorizontalAlignment = 'center';
            app.text25.VerticalAlignment = 'top';
            app.text25.FontSize = 13.3333333333333;
            app.text25.Position = [191 88 80.7511254019292 19.7189189189189];
            app.text25.Text = '% correct';

            % Create text26
            app.text26 = uilabel(app.uipanel4);
            app.text26.Tag = 'text26';
            app.text26.HorizontalAlignment = 'center';
            app.text26.VerticalAlignment = 'top';
            app.text26.FontSize = 13.3333333333333;
            app.text26.Position = [186 61 80.7511254019293 19.7189189189189];
            app.text26.Text = '% avers';

            % Create edit26
            app.edit26 = uieditfield(app.uipanel4, 'text');
            app.edit26.Tag = 'edit26';
            app.edit26.HorizontalAlignment = 'center';
            app.edit26.FontSize = 13.3333333333333;
            app.edit26.Position = [277 58 42.6186495176849 24.1009009009009];
            app.edit26.Value = '0';

            % Create text30
            app.text30 = uilabel(app.uipanel4);
            app.text30.Tag = 'text30';
            app.text30.HorizontalAlignment = 'center';
            app.text30.VerticalAlignment = 'top';
            app.text30.FontSize = 13.3333333333333;
            app.text30.Position = [206 35 57.1987138263665 19.7189189189189];
            app.text30.Text = 'Day';

            % Create edit29
            app.edit29 = uieditfield(app.uipanel4, 'text');
            app.edit29.Tag = 'edit29';
            app.edit29.HorizontalAlignment = 'center';
            app.edit29.FontSize = 13.3333333333333;
            app.edit29.Position = [277 32 42.6186495176849 24.1009009009009];
            app.edit29.Value = '1';

            % Create uipanel6
            app.uipanel6 = uipanel(app.figure1);
            app.uipanel6.Title = 'ROI';
            app.uipanel6.Tag = 'uipanel6';
            app.uipanel6.FontSize = 18.6666666666667;
            app.uipanel6.Position = [1113 1 308 361];

            % Create pushbutton6
            app.pushbutton6 = uibutton(app.uipanel6, 'push');
            app.pushbutton6.ButtonPushedFcn = createCallbackFcn(app, @pushbutton6_Callback, true);
            app.pushbutton6.Tag = 'pushbutton6';
            app.pushbutton6.FontSize = 10.6666666666667;
            app.pushbutton6.Position = [28 73 114.24 36.6081475787855];
            app.pushbutton6.Text = 'Load Image';

            % Create edit14
            app.edit14 = uieditfield(app.uipanel6, 'text');
            app.edit14.Tag = 'edit14';
            app.edit14.HorizontalAlignment = 'center';
            app.edit14.FontSize = 13.3333333333333;
            app.edit14.Position = [28 238 249.033210332103 24.4541935483871];
            app.edit14.Value = 'C:\Valerio\M01_000_006';

            % Create text15
            app.text15 = uilabel(app.uipanel6);
            app.text15.Tag = 'text15';
            app.text15.HorizontalAlignment = 'center';
            app.text15.VerticalAlignment = 'top';
            app.text15.FontSize = 13.3333333333333;
            app.text15.Position = [74 272 141 16];
            app.text15.Text = 'Add filename here';

            % Create edit17
            app.edit17 = uieditfield(app.uipanel6, 'text');
            app.edit17.ValueChangedFcn = createCallbackFcn(app, @edit17_Callback, true);
            app.edit17.Tag = 'edit17';
            app.edit17.HorizontalAlignment = 'center';
            app.edit17.FontSize = 13.3333333333333;
            app.edit17.Position = [39 190 44.870848708487 23.3909677419355];
            app.edit17.Value = '100';

            % Create text17
            app.text17 = uilabel(app.uipanel6);
            app.text17.Tag = 'text17';
            app.text17.HorizontalAlignment = 'center';
            app.text17.VerticalAlignment = 'top';
            app.text17.FontSize = 13.3333333333333;
            app.text17.Position = [28 218 63.940959409594 14.8851612903226];
            app.text17.Text = '# frames';

            % Create edit18
            app.edit18 = uieditfield(app.uipanel6, 'text');
            app.edit18.Tag = 'edit18';
            app.edit18.HorizontalAlignment = 'center';
            app.edit18.FontSize = 13.3333333333333;
            app.edit18.Position = [130 188 45.9926199261992 23.3909677419355];
            app.edit18.Value = '1';

            % Create edit19
            app.edit19 = uieditfield(app.uipanel6, 'text');
            app.edit19.Tag = 'edit19';
            app.edit19.HorizontalAlignment = 'center';
            app.edit19.FontSize = 13.3333333333333;
            app.edit19.Position = [212 189 44.8708487084871 23.3909677419355];
            app.edit19.Value = '1';

            % Create text18
            app.text18 = uilabel(app.uipanel6);
            app.text18.Tag = 'text18';
            app.text18.HorizontalAlignment = 'center';
            app.text18.VerticalAlignment = 'top';
            app.text18.FontSize = 13.3333333333333;
            app.text18.Position = [121 218 66.1845018450184 14.8851612903226];
            app.text18.Text = 'Plane #';

            % Create text19
            app.text19 = uilabel(app.uipanel6);
            app.text19.Tag = 'text19';
            app.text19.HorizontalAlignment = 'center';
            app.text19.VerticalAlignment = 'top';
            app.text19.FontSize = 13.3333333333333;
            app.text19.Position = [201 212 76.280442804428 21.2645161290323];
            app.text19.Text = 'Plane tot.';

            % Create uibuttongroup8
            app.uibuttongroup8 = uibuttongroup(app.uipanel6);
            app.uibuttongroup8.Tag = 'uibuttongroup8';
            app.uibuttongroup8.FontSize = 10.6666666666667;
            app.uibuttongroup8.Position = [76 126 137 52];

            % Create radiobutton14
            app.radiobutton14 = uiradiobutton(app.uibuttongroup8);
            app.radiobutton14.Tag = 'radiobutton14';
            app.radiobutton14.Text = 'Max';
            app.radiobutton14.FontSize = 10.6666666666667;
            app.radiobutton14.Position = [4 5 70.9333333333334 24.2699892133024];

            % Create radiobutton13
            app.radiobutton13 = uiradiobutton(app.uibuttongroup8);
            app.radiobutton13.Tag = 'radiobutton13';
            app.radiobutton13.Text = 'All frames';
            app.radiobutton13.FontSize = 10.6666666666667;
            app.radiobutton13.Position = [63 5 86.4500000000001 24.2699892133024];

            % Create radiobutton12
            app.radiobutton12 = uiradiobutton(app.uibuttongroup8);
            app.radiobutton12.Tag = 'radiobutton12';
            app.radiobutton12.Text = 'Mean';
            app.radiobutton12.FontSize = 10.6666666666667;
            app.radiobutton12.Position = [4 25 84.2333333333334 24.2699892133024];
            app.radiobutton12.Value = true;

            % Create radiobutton15
            app.radiobutton15 = uiradiobutton(app.uibuttongroup8);
            app.radiobutton15.Tag = 'radiobutton15';
            app.radiobutton15.Text = 'Std';
            app.radiobutton15.FontSize = 10.6666666666667;
            app.radiobutton15.Position = [63 25 96.4250000000001 24.2699892133024];

            % Create pushbutton10
            app.pushbutton10 = uibutton(app.uipanel6, 'push');
            app.pushbutton10.ButtonPushedFcn = createCallbackFcn(app, @pushbutton10_Callback, true);
            app.pushbutton10.Tag = 'pushbutton10';
            app.pushbutton10.FontSize = 10.6666666666667;
            app.pushbutton10.Position = [30 32 114.24 37.684857801691];
            app.pushbutton10.Text = 'Load ROI';

            % Create pushbutton11
            app.pushbutton11 = uibutton(app.uipanel6, 'push');
            app.pushbutton11.ButtonPushedFcn = createCallbackFcn(app, @pushbutton11_Callback, true);
            app.pushbutton11.Tag = 'pushbutton11';
            app.pushbutton11.FontSize = 10.6666666666667;
            app.pushbutton11.Position = [165 73 114.24 37.684857801691];
            app.pushbutton11.Text = 'Draw ROI';

            % Create pushbutton12
            app.pushbutton12 = uibutton(app.uipanel6, 'push');
            app.pushbutton12.ButtonPushedFcn = createCallbackFcn(app, @pushbutton12_Callback, true);
            app.pushbutton12.Tag = 'pushbutton12';
            app.pushbutton12.FontSize = 10.6666666666667;
            app.pushbutton12.Position = [165 35 114.24 37.684857801691];
            app.pushbutton12.Text = 'Save ROI';

            % Create checkbox13
            app.checkbox13 = uicheckbox(app.uipanel6);
            app.checkbox13.Tag = 'checkbox13';
            app.checkbox13.Text = 'Allignment Only';
            app.checkbox13.FontSize = 12;
            app.checkbox13.Position = [28 289 137 25];

            % Create slider1
            app.slider1 = uislider(app.figure1);
            app.slider1.MajorTicks = [];
            app.slider1.ValueChangedFcn = createCallbackFcn(app, @slider1_Callback, true);
            app.slider1.MinorTicks = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100];
            app.slider1.FontSize = 10.6666666666667;
            app.slider1.Tag = 'slider1';
            app.slider1.Position = [391.88 382.846153846154 1029.28 3];
            app.slider1.Value = 1;

            % Create uipanel6_2
            app.uipanel6_2 = uipanel(app.figure1);
            app.uipanel6_2.Title = 'MC Panel';
            app.uipanel6_2.Tag = 'uipanel6';
            app.uipanel6_2.FontSize = 18.6666666666667;
            app.uipanel6_2.Position = [31 408 175 109];

            % Create edit9_2
            app.edit9_2 = uieditfield(app.uipanel6_2, 'text');
            app.edit9_2.Tag = 'edit9';
            app.edit9_2.HorizontalAlignment = 'center';
            app.edit9_2.FontSize = 10.6666666666667;
            app.edit9_2.Position = [48 50 30.3230769230769 21.2432432432432];
            app.edit9_2.Value = '3';

            % Create text10_2
            app.text10_2 = uilabel(app.uipanel6_2);
            app.text10_2.Tag = 'text10';
            app.text10_2.HorizontalAlignment = 'center';
            app.text10_2.VerticalAlignment = 'top';
            app.text10_2.FontSize = 10.6666666666667;
            app.text10_2.Position = [20 50 31 17];
            app.text10_2.Text = 'Nq';

            % Create edit9_3
            app.edit9_3 = uieditfield(app.uipanel6_2, 'text');
            app.edit9_3.Tag = 'edit9';
            app.edit9_3.HorizontalAlignment = 'center';
            app.edit9_3.FontSize = 10.6666666666667;
            app.edit9_3.Position = [48 18 30.3230769230769 21.2432432432432];
            app.edit9_3.Value = '200';

            % Create text10_3
            app.text10_3 = uilabel(app.uipanel6_2);
            app.text10_3.Tag = 'text10';
            app.text10_3.HorizontalAlignment = 'center';
            app.text10_3.VerticalAlignment = 'top';
            app.text10_3.FontSize = 10.6666666666667;
            app.text10_3.Position = [18 13 31 22];
            app.text10_3.Text = 'NqX';

            % Create edit9_4
            app.edit9_4 = uieditfield(app.uipanel6_2, 'text');
            app.edit9_4.Tag = 'edit9';
            app.edit9_4.HorizontalAlignment = 'center';
            app.edit9_4.FontSize = 10.6666666666667;
            app.edit9_4.Position = [123 18 30.3230769230769 21.2432432432432];
            app.edit9_4.Value = '128';

            % Create text10_4
            app.text10_4 = uilabel(app.uipanel6_2);
            app.text10_4.Tag = 'text10';
            app.text10_4.HorizontalAlignment = 'center';
            app.text10_4.VerticalAlignment = 'top';
            app.text10_4.FontSize = 10.6666666666667;
            app.text10_4.Position = [93 13 31 22];
            app.text10_4.Text = 'NqY';

            % Create edit9_5
            app.edit9_5 = uieditfield(app.uipanel6_2, 'text');
            app.edit9_5.Tag = 'edit9';
            app.edit9_5.HorizontalAlignment = 'center';
            app.edit9_5.FontSize = 10.6666666666667;
            app.edit9_5.Position = [123 50 30.3230769230769 21.2432432432432];
            app.edit9_5.Value = '0';

            % Create text10_5
            app.text10_5 = uilabel(app.uipanel6_2);
            app.text10_5.Tag = 'text10';
            app.text10_5.HorizontalAlignment = 'center';
            app.text10_5.VerticalAlignment = 'top';
            app.text10_5.FontSize = 10.6666666666667;
            app.text10_5.Position = [92 45 33 22];
            app.text10_5.Text = 'Dsmp';

            % Create uipanel6_3
            app.uipanel6_3 = uipanel(app.figure1);
            app.uipanel6_3.Title = 'Exp design';
            app.uipanel6_3.Tag = 'uipanel6';
            app.uipanel6_3.FontSize = 18.6666666666667;
            app.uipanel6_3.Position = [216 409 169 108];

            % Create uibuttongroup8_2
            app.uibuttongroup8_2 = uibuttongroup(app.uipanel6_3);
            app.uibuttongroup8_2.Tag = 'uibuttongroup8';
            app.uibuttongroup8_2.FontSize = 10.6666666666667;
            app.uibuttongroup8_2.Position = [19 8 130 66];

            % Create radiobutton14_2
            app.radiobutton14_2 = uiradiobutton(app.uibuttongroup8_2);
            app.radiobutton14_2.Tag = 'radiobutton14';
            app.radiobutton14_2.Text = '2 Pop';
            app.radiobutton14_2.FontSize = 10.6666666666667;
            app.radiobutton14_2.Position = [4 19 71 24];

            % Create radiobutton12_2
            app.radiobutton12_2 = uiradiobutton(app.uibuttongroup8_2);
            app.radiobutton12_2.Tag = 'radiobutton12';
            app.radiobutton12_2.Text = '1 Pop';
            app.radiobutton12_2.FontSize = 10.6666666666667;
            app.radiobutton12_2.Position = [4 39 84 24];
            app.radiobutton12_2.Value = true;

            % Create radiobutton15_2
            app.radiobutton15_2 = uiradiobutton(app.uibuttongroup8_2);
            app.radiobutton15_2.Tag = 'radiobutton15';
            app.radiobutton15_2.Text = 'Pop Dynamics';
            app.radiobutton15_2.FontSize = 10.6666666666667;
            app.radiobutton15_2.Position = [4 -2 97 24];

            % Create ROIActivityPanel
            app.ROIActivityPanel = uipanel(app.figure1);
            app.ROIActivityPanel.Title = 'ROI Activity';
            app.ROIActivityPanel.FontSize = 18;
            app.ROIActivityPanel.Position = [415 1 688 365];

            % Create axes3
            app.axes3 = uiaxes(app.ROIActivityPanel);
            app.axes3.FontSize = 13.3333333333333;
            app.axes3.NextPlot = 'replace';
            app.axes3.Tag = 'axes3';
            app.axes3.Position = [20 19 658 319];

            % Create uipanel6_4
            app.uipanel6_4 = uipanel(app.figure1);
            app.uipanel6_4.Title = 'Image Panel';
            app.uipanel6_4.Tag = 'uipanel6';
            app.uipanel6_4.FontSize = 18.6666666666667;
            app.uipanel6_4.Position = [31 290 175 109];

            % Create edit20
            app.edit20 = uieditfield(app.uipanel6_4, 'text');
            app.edit20.Tag = 'edit20';
            app.edit20.HorizontalAlignment = 'center';
            app.edit20.FontSize = 10.6666666666667;
            app.edit20.Tooltip = 'Target Activation';
            app.edit20.Position = [58 14 30 21];
            app.edit20.Value = '512';

            % Create edit21
            app.edit21 = uieditfield(app.uipanel6_4, 'text');
            app.edit21.Tag = 'edit21';
            app.edit21.HorizontalAlignment = 'center';
            app.edit21.FontSize = 10.6666666666667;
            app.edit21.Tooltip = 'Target Activation';
            app.edit21.Position = [56 51 30 21];
            app.edit21.Value = '796';

            % Create text20
            app.text20 = uilabel(app.uipanel6_4);
            app.text20.Tag = 'text20';
            app.text20.HorizontalAlignment = 'center';
            app.text20.VerticalAlignment = 'top';
            app.text20.FontSize = 10.6666666666667;
            app.text20.Tooltip = 'Target Activation';
            app.text20.Position = [1 51 58.4 19.1189189189189];
            app.text20.Text = 'Pixels (x)';

            % Create text21
            app.text21 = uilabel(app.uipanel6_4);
            app.text21.Tag = 'text21';
            app.text21.HorizontalAlignment = 'center';
            app.text21.VerticalAlignment = 'top';
            app.text21.FontSize = 10.6666666666667;
            app.text21.Tooltip = 'Target Activation';
            app.text21.Position = [4 14 58.4 19.1189189189189];
            app.text21.Text = 'Lines (y)';

            % Create text21_2
            app.text21_2 = uilabel(app.uipanel6_4);
            app.text21_2.Tag = 'text21';
            app.text21_2.HorizontalAlignment = 'center';
            app.text21_2.VerticalAlignment = 'top';
            app.text21_2.FontSize = 10.6666666666667;
            app.text21_2.Tooltip = 'Target Activation';
            app.text21_2.Position = [86 48 43 23];
            app.text21_2.Text = 'Quality';

            % Create edit20_2
            app.edit20_2 = uieditfield(app.uipanel6_4, 'text');
            app.edit20_2.Tag = 'edit20';
            app.edit20_2.HorizontalAlignment = 'center';
            app.edit20_2.FontSize = 10.6666666666667;
            app.edit20_2.Tooltip = 'Target Activation';
            app.edit20_2.Position = [128 52 42 21];
            app.edit20_2.Value = 'uint16';

            % Create edit20_3
            app.edit20_3 = uieditfield(app.uipanel6_4, 'text');
            app.edit20_3.HorizontalAlignment = 'center';
            app.edit20_3.FontSize = 10.6666666666667;
            app.edit20_3.Tooltip = 'Target Activation';
            app.edit20_3.Position = [141 16 30 21];
            app.edit20_3.Value = '1';

            % Create text21_3
            app.text21_3 = uilabel(app.uipanel6_4);
            app.text21_3.HorizontalAlignment = 'center';
            app.text21_3.VerticalAlignment = 'top';
            app.text21_3.FontSize = 10.6666666666667;
            app.text21_3.Tooltip = 'Target Activation';
            app.text21_3.Position = [87 13 58 22];
            app.text21_3.Text = 'Dsmp';

            % Create uipanel6_5
            app.uipanel6_5 = uipanel(app.figure1);
            app.uipanel6_5.Title = 'Feedback';
            app.uipanel6_5.Tag = 'uipanel6';
            app.uipanel6_5.FontSize = 18.6666666666667;
            app.uipanel6_5.Position = [215 290 169 108];

            % Create edit7
            app.edit7 = uieditfield(app.uipanel6_5, 'text');
            app.edit7.Tag = 'edit7';
            app.edit7.HorizontalAlignment = 'center';
            app.edit7.FontSize = 10.6666666666667;
            app.edit7.Position = [78 22 30.3230769230769 21.2432432432432];
            app.edit7.Value = '6';

            % Create edit13
            app.edit13 = uieditfield(app.uipanel6_5, 'text');
            app.edit13.Tag = 'edit13';
            app.edit13.HorizontalAlignment = 'center';
            app.edit13.FontSize = 10.6666666666667;
            app.edit13.Position = [78 1 30.3230769230769 20.1810810810811];
            app.edit13.Value = '90';

            % Create text13
            app.text13 = uilabel(app.uipanel6_5);
            app.text13.Tag = 'text13';
            app.text13.HorizontalAlignment = 'center';
            app.text13.VerticalAlignment = 'top';
            app.text13.FontSize = 10.6666666666667;
            app.text13.Position = [5 3 68.5076923076923 15.9324324324324];
            app.text13.Text = 'Angle';

            % Create text8
            app.text8 = uilabel(app.uipanel6_5);
            app.text8.Tag = 'text8';
            app.text8.HorizontalAlignment = 'center';
            app.text8.VerticalAlignment = 'top';
            app.text8.FontSize = 10.6666666666667;
            app.text8.Position = [5 23 67.3846153846154 15.9324324324324];
            app.text8.Text = 'Spatial freq.';

            % Create uibuttongroup8_3
            app.uibuttongroup8_3 = uibuttongroup(app.uipanel6_5);
            app.uibuttongroup8_3.Tag = 'uibuttongroup8';
            app.uibuttongroup8_3.FontSize = 10.6666666666667;
            app.uibuttongroup8_3.Position = [5 46 158 30];

            % Create radiobutton14_3
            app.radiobutton14_3 = uiradiobutton(app.uibuttongroup8_3);
            app.radiobutton14_3.Tag = 'radiobutton14';
            app.radiobutton14_3.Text = 'Audio';
            app.radiobutton14_3.FontSize = 10.6666666666667;
            app.radiobutton14_3.Position = [104 2 49 24];

            % Create radiobutton12_3
            app.radiobutton12_3 = uiradiobutton(app.uibuttongroup8_3);
            app.radiobutton12_3.Tag = 'radiobutton12';
            app.radiobutton12_3.Text = 'Visual';
            app.radiobutton12_3.FontSize = 10.6666666666667;
            app.radiobutton12_3.Position = [4 3 53 24];
            app.radiobutton12_3.Value = true;

            % Create uipanel6_6
            app.uipanel6_6 = uipanel(app.figure1);
            app.uipanel6_6.Title = 'Signal';
            app.uipanel6_6.FontSize = 18;
            app.uipanel6_6.Position = [31 204 175 71];

            % Create uibuttongroup8_4
            app.uibuttongroup8_4 = uibuttongroup(app.uipanel6_6);
            app.uibuttongroup8_4.Position = [5 2 158 30];

            % Create radiobutton14_4
            app.radiobutton14_4 = uiradiobutton(app.uibuttongroup8_4);
            app.radiobutton14_4.Text = 'Spikes';
            app.radiobutton14_4.FontSize = 10.6666666666667;
            app.radiobutton14_4.Position = [104 2 54 24];

            % Create radiobutton12_4
            app.radiobutton12_4 = uiradiobutton(app.uibuttongroup8_4);
            app.radiobutton12_4.Text = 'DF/F0';
            app.radiobutton12_4.FontSize = 10.6666666666667;
            app.radiobutton12_4.Position = [4 3 53 24];
            app.radiobutton12_4.Value = true;

            % Show the figure after all components are created
            app.figure1.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Run_GUI_App(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.figure1)

            % Execute the startup function
            runStartupFcn(app, @(app)Run_GUI_OpeningFcn(app, varargin{:}))

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.figure1)
        end
    end
end