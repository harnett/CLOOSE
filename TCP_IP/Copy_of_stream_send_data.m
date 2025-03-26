clear snd
snd = tcpip('10.93.6.226', 30001, 'NetworkRole', 'server', 'OutputBufferSize',(512*796*20), 'Terminator', 'CR/LF');
fopen(snd);
xpx = 512;
img = uint16(zeros(796*xpx, 1));
for im = 1:10000
    fwrite(snd, img, 'uint8')
    
end
fclose(snd);
clear snd