clear all %#ok<CLALL>  % Clear all variables from the workspace
close all             % Close all open figure windows

%% INSERT YOUR RECEIVER IPv4 ADDRESS HERE
IPv4 = '127.0.0.1';

%% RUN THE REST
% Define Image Dimensions and Preallocate Storage
xpixels = 512;         % Image height (number of rows)
ypixels = 796;         % Image width (number of columns)
frames  = 500;         % Total number of frames expected to receive
% Preallocate a 3D array (uint16) to store the incoming image frames.
received_Data = uint16(zeros(xpixels, ypixels, frames));

% Set Up TCP/IP Connection Parameters
port = 30001;          % Port number used for communication
buffsz = 20;           % Buffer size multiplier to ensure the input buffer is large enough

% Create a TCP/IP client object to receive data from the local host.
% '127.0.0.1' is the loopback IP address (localhost), which is used for same-machine communication.
% 'InputBufferSize' is set to accommodate the expected number of pixels per frame times the buffer multiplier.
% 'Terminator' specifies the end-of-line characters used in the communication protocol.
tcpipServer = tcpip(IPv4, port, 'NetworkRole', 'server', ...
    'InputBufferSize', (xpixels*ypixels*buffsz), 'Terminator', 'CR/LF');

% Open the TCP/IP connection.
fopen(tcpipServer);

% Initialize Figure for Displaying Incoming Images
% Create a figure window for displaying the received images.
hFig = figure(1);
% Create an axes object in the figure.
hAx  = axes('Parent', hFig);
% Initialize the display with a blank image.
% The image object handle (hImg) will be used to update the displayed image data.
hImg = imshow(zeros(xpixels, ypixels, 'uint16'), 'Parent', hAx);

% Receive and Display Frames in a Loop
for iframe = 1:frames 
    % Read image data from the TCP/IP connection.
    % The number of elements to read is equal to the total number of pixels in the image.
    vect_img = uint16(fread(tcpipServer, xpixels*ypixels, 'uint16'));
    
    % Reshape the received vector into a 2D image frame.
    imgFrame = reshape(vect_img, xpixels, ypixels);
    
    % Store the current image frame in the preallocated 3D array.
    received_Data(:, :, iframe) = imgFrame;

    % Update the image display by modifying the 'CData' property of the existing image object.
    set(hImg, 'CData', imgFrame);
    % Refresh the display, limiting the rate of updates for efficiency.
    drawnow limitrate;
end

% Clean Up the TCP/IP Connection
% Close the TCP/IP connection once all frames have been received.
fclose(tcpipServer);
% Clear the TCP/IP object from the workspace.
clear tcpipServer
