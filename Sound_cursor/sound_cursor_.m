function continuous_tone_cursor
    % Global variables to share state between callbacks.
    global currentFreq phase ap dt fs fig ax;
    
    % Initialize frequency and phase.
    currentFreq = 200; % Starting frequency (Hz)
    phase = 0;         % Initial phase for continuity
    
    % Audio parameters.
    fs = 44100;        % Sampling rate (Hz)
    dt = 0.05;         % Duration (seconds) of each audio block (10 ms)
    
    % Create a figure with an axes and set up the mouse motion callback.
    fig = figure('Name','Continuous Tone Control','NumberTitle','off',...
        'WindowButtonMotionFcn',@updateFrequency);
    ax = axes('Parent',fig);
    xlim(ax, [0 1]);   % x-range mapped to frequency changes
    ylim(ax, [0 1]);   % dummy y-range
    title(ax, 'Move mouse horizontally to change frequency');
    
    % Create the audio device writer (from Audio Toolbox)
    ap = audioDeviceWriter('SampleRate', fs);
    
    % Set up a timer to generate and play audio blocks at fixed intervals.
    tmr = timer('ExecutionMode','fixedRate', 'Period', dt, 'TimerFcn', @playChunk);
    start(tmr);
    
    % Ensure proper cleanup when the figure is closed.
    fig.CloseRequestFcn = @(src, event) closeFigure(src, event, tmr);
end

% Timer callback: generates a short sine-wave block and outputs audio.
function playChunk(~, ~)
    global currentFreq phase ap dt fs;
    
    % Use the current frequency to generate a sine wave block.
    t = (0:dt*fs-1)/fs;
    y = sin(2*pi*currentFreq*t + phase);
    
    % Update phase to maintain continuity.
    phase = phase + 2*pi*currentFreq*dt;
    
    % Output the audio block. (audioDeviceWriter expects a column vector.)
    step(ap, y');
end

% Callback for mouse motion: update frequency based on cursor's x-position.
function updateFrequency(~, ~)
    global currentFreq;
    
    % Get the current mouse position in axes coordinates.
    cp = get(gca, 'CurrentPoint');
    x = cp(1,1);
    
    % Get the current axes limits and clamp x.
    ax = gca;
    xLimits = get(ax, 'XLim');
    x = max(min(x, xLimits(2)), xLimits(1));
    
    % Map the x coordinate linearly to a frequency range [200, 2000] Hz.
    currentFreq = 200 + (2000-200)*((x - xLimits(1))/(xLimits(2)-xLimits(1)));
end

% Cleanup function: stops and deletes the timer and releases the audio device.
function closeFigure(src, ~, tmr)
    global ap;
    stop(tmr);
    delete(tmr);
    release(ap);
    delete(src);
end
