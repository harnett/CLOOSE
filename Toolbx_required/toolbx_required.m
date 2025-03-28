cd('C:\Users\MysticSocialist\Documents\GitHub\CLOOSE_MS');

% Point to the main script or entry point of your project
[~, requiredToolboxes] = matlab.codetools.requiredFilesAndProducts('stream_send_data.m');

% Display required toolboxes
for k = 1:length(requiredToolboxes)
    disp(requiredToolboxes(k).Name)
end