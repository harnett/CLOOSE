% Continuous Tone Controlled by Mouse Movement
% Move your cursor horizontally in the figure to change the tone frequency.

% Define a global variable to hold the current frequency.
global currentFreq;
currentFreq = 200;  % Start at 200 Hz

% Audio parameters
fs = 44100;       % Sampling rate in Hz
dt = 0.05;        % Duration (in seconds) of each audio chunk
phase = 0;        % Initial phase for continuous sine wave

% Create a figure with an axes.
fig = figure('Name', 'Continuous Tone Control', 'NumberTitle', 'off', ...
             'WindowButtonMotionFcn', @updateFrequency);
ax = axes('Parent', fig);
xlim(ax, [0 1]);  % x-range that will map to frequency changes
ylim(ax, [0 1]);  % dummy y-range
title(ax, 'Move mouse horizontally to change frequency');

% Create the audio player object (requires DSP System Toolbox)
ap = dsp.AudioPlayer('SampleRate', fs);

% Main loop: while the figure is open, continuously generate and play sound.
while ishandle(fig)
    % Get the current frequency (updated by the mouse callback).
    freq = currentFreq;
    
    % Generate a time vector and compute the sine wave.
    t = (0:dt*fs-1) / fs;
    y = sin(2*pi*freq*t + phase);
    
    % Update phase so that the sine wave is continuous.
    phase = phase + 2*pi*freq*dt;
    
    % Output the audio chunk.
    step(ap, y');
    
    % Allow MATLAB to process figure callbacks.
    drawnow;
end

% Clean up the audio player once the figure is closed.
release(ap);

% --- Callback function ---
function updateFrequency(~, ~)
    % This callback updates the global frequency based on the x-position 
    % of the mouse within the axes.
    global currentFreq;
    
    % Get the current mouse position in axes coordinates.
    cp = get(gca, 'CurrentPoint');
    x = cp(1,1);
    
    % Retrieve the current x-limits of the axes.
    ax = gca;
    xLimits = get(ax, 'XLim');
    
    % Clamp x between the x-limits.
    x = max(min(x, xLimits(2)), xLimits(1));
    
    % Map the x coordinate linearly to a frequency range [200, 2000] Hz.
    currentFreq = 200 + (2000-200) * ((x - xLimits(1)) / (xLimits(2) - xLimits(1)));
end
