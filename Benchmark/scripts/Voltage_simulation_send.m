clear all %#ok<CLALL>  % Clear all workspace variables
close all


filename = 'C:\Users\Harnettlab\OneDrive\Desktop\Voltage_Data\Sq_camera.bin';
nrow = 96;
ncol = 200;
stack = readBinMov(filename, nrow, ncol);

IPv4 = '127.0.0.1';

% Define TCP/IP communication parameters.
port = 30001;  % Set the port number for communication.
buffsz = 20;   % Buffer size multiplier for the output buffer.

% Create a TCP/IP server object to send data.
% '127.0.0.1' is the loopback address (localhost), meaning the connection is on the same machine.
% 'NetworkRole' is set to 'server' for this instance.
% 'OutputBufferSize' is determined by the image dimensions multiplied by the buffer size factor.
% 'Terminator' is set to 'CR/LF' (Carriage Return/Line Feed).
disp('ready')

send_device = tcpip(IPv4, port, 'NetworkRole', 'server', ...
    'OutputBufferSize', (size(stack, 1) * size(stack, 2) * buffsz), 'Terminator', 'CR/LF');

% Open the TCP/IP connection.
fopen(send_device);

% Loop through each frame in the image stack.
for iframe = 1:size(stack, 3)
    % Reshape the 2D image frame into a column vector.
    % This prepares the data for transmission over the TCP/IP connection.
    vect_img = reshape(stack(:, :, iframe), [], 1);
    
    % Write the image data to the TCP/IP connection as uint16.
    fwrite(send_device, vect_img, 'uint16');
end

% Close the TCP/IP connection after all frames have been sent.
fclose(send_device);

% Clear the TCP/IP object from the workspace.
clear send_device
