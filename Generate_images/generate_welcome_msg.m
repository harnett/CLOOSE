%% Parameters
nFrames   = 500;
imgHeight = 512;
imgWidth  = 796;
stack     = zeros(imgHeight, imgWidth, nFrames, 'uint16');  % Preallocate image stack

%% Create Final Text Image
% Create a blank image for text insertion.
baseImg = zeros(imgHeight, imgWidth, 'uint8');

% Define the text string.
textStr = 'WELCOME TO CLOOSE, BY THE HARNETT LAB';

% Center the text; note: insertText uses [x y] where x is horizontal.
pos = [imgWidth/2, imgHeight/2];

% Insert the text. 'BoxOpacity' is set to 0 to avoid a background box.
imgWithText = insertText(baseImg, pos, textStr, ...
    'FontSize', 30, 'BoxOpacity', 0, 'AnchorPoint', 'Center', 'TextColor', 'white');

% Convert the RGB image to grayscale.
grayImg = rgb2gray(imgWithText);

% Create a binary mask for the text (threshold = 128).
textMask = grayImg > 128;

%% Compute the Union Bounding Box of the Text
% Get coordinates of all text pixels.
[yCoords, xCoords] = find(textMask);
if isempty(xCoords)
    error('No text found.');
end
minX = min(xCoords);
maxX = max(xCoords);
x0 = minX;             % Leftmost column of the text
widthBB = maxX - minX; % Width of the text region

% Precompute a grid of column indices.
[X, ~] = meshgrid(1:imgWidth, 1:imgHeight);

%% Generate the Gradually Written Animation
% Over 10,000 frames, reveal the text from left (x0) to right (x0 + widthBB)
for i = 1:nFrames
    progress = i / nFrames;      % Progress from 0 to 1
    currentThreshold = x0 + widthBB * progress;
    % Reveal text pixels whose column index is less than or equal to the threshold.
    revealMask = (X <= currentThreshold) & textMask;
    % Create the frame: white for revealed text, black otherwise.
    frame = uint16(revealMask) * 65535;
    stack(:, :, i) = frame;
    
    figure(1)
    imshow(frame)
    drawnow limitrate
end
