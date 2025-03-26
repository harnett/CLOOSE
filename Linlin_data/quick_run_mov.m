%% load mov

DaqRate = 10000; sz=get(0,'screensize');
Info=textscan(fopen('experimental_parameters.txt'),'%s');
nrow = str2num(Info{1,1}{6,1}); ncol = str2num(Info{1,1}{3,1});
binName='movReg.bin'; [mov nframes] = readBinMov(binName, ncol, nrow);

% 
[ROI,Fclicky]=clicky_faster(mov); 