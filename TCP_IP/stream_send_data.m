clear all %#ok<CLALL> 

send_device = tcpip('10.93.6.226', 30001, 'NetworkRole', 'server', 'OutputBufferSize',(512*796*20), 'Terminator', 'CR/LF');
fopen(send_device);

xpx = 512;
ypx = 796;
img = uint16(zeros(ypx*xpx, 1));

for im = 1:10000
    fwrite(send_device, img, 'uint8')
end

fclose(send_device);
clear send_device