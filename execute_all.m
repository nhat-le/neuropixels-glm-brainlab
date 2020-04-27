for i = 12
%     try
    files = dir('*multiNeurons*.m');

    eval(files(i).name(1:end-2));
%     catch
%         errors = [errors i];
%     end
    
%     close all;
%     clear all;
end

