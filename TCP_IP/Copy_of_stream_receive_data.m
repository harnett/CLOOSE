% tcpipServer = tcpip('10.93.15.144',30000,'NetworkRole','client');
% tcpipServer = tcpip('10.93.14.143',30000,'NetworkRole','client');
close all
clear all

path = 'C:\Users\Harnettlab\OneDrive\Desktop\testing_bci_bench_marktest\';
name = 'fullimage_uint16_external.mat';

xpixels = 512; 
ypixels = 796;
z = 10000;

tcpipServer = tcpip('10.93.6.2',30001,'NetworkRole','client', 'Terminator', 'CR/LF');
tcpipServer.InputBufferSize = (xpixels*ypixels*20);

fopen(tcpipServer);
tt = nan(1, z);
for i = 1:z 
    tic
    received_Data = uint16(fread(tcpipServer, ypixels*xpixels, 'uint16'));
    tt(i) = toc;
end

fclose(tcpipServer);

save([path name], 'tt')