clear all %#ok<CLALL>
close all

% Create figure and axes.
fig = figure('Name', 'Flexible dt Tone Control (No Callbacks)', 'NumberTitle', 'off');
ax = axes('Parent', fig);
xlim(ax, [0 1]);    % x-axis range that maps to frequency
ylim(ax, [0 1]);    % Dummy y-axis range
title(ax, 'Move mouse horizontally to change frequency');

% Audio parameters.
fs = 44100;         % Sampling rate in Hz
phase = 0;          % Initial phase for continuous sine wave

% Create the audio device writer object with variable-size input support.
ap = audioDeviceWriter('SampleRate', fs, 'SupportVariableSizeInput', true);

% Main loop: continuously generate and play sound until the figure is closed.
prevTime = tic;
while ishandle(fig)
    % Determine dt as the time elapsed since last loop iteration.
    dt = toc(prevTime);
    prevTime = tic;
    if dt < 0.01, dt = 0.01; end  % Avoid extremely short durations

    % Poll the current mouse position directly from the axes.
    cp = get(ax, 'CurrentPoint');
    x = cp(1,1);
    
    % Clamp x to the x-axis limits.
    xLimits = get(ax, 'XLim');
    x = max(min(x, xLimits(2)), xLimits(1));
    
    % Map x to a frequency range [200, 2000] Hz.
    freq = 200 + (2000 - 200) * ((x - xLimits(1)) / (xLimits(2) - xLimits(1)));
    
    % Generate a sine wave for dt seconds.
    numSamples = round(dt * fs);
    t = (0:numSamples-1) / fs;
    y = sin(2*pi*freq*t + phase);
    
    % Update phase to maintain continuity.
    phase = phase + 2*pi*freq*dt;
    
    % Play the generated audio chunk.
    ap(y');
    
    % Process pending graphics events.
    drawnow;
end

% Release the audio device writer when done.
release(ap);
