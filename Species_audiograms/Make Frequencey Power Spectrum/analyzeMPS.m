%Moises R. (4/24/2019): Runs the plotMPS_MR for all stim files (.mat from plotMPS_MR.m) selected (MUST be more than one file)

handles.stimDir = 'C:\Users\Woolley Lab-MB\Desktop\TEST'; %where to select your wav files
handles.figDir = 'C:\Users\Woolley Lab-MB\Desktop\TEST\figs'; %where MPS will be saved

conMPS = struct;
%Load mpsRel file into compiled structure conMPS
[fileNames, pathName] = uigetfile([handles.stimDir,'\*.mat'],'Choose files', 'MultiSelect','on');
nFiles = length(fileNames);

for s = 1:nFiles

    handles.input = fileNames(s);
    fileChar = char(handles.input);
    fileChar2 = char([pathName fileChar]);
    conMPS(s).File = fileChar(1:end-4); %populates the file names in the structure
    
    if pathName==0, error('None selected!');
        else conMPS(s).relMPS = load(fileChar2,'mpsRel'); %loads the relative MPS matrix of each file respectively into the structure
    end    
end
save([handles.figDir 'compiledMPS.mat']); %saves the compiled MPS .mat file (conMPS)


vecMPS = struct;
%converts MPS rel values matrices into vectors for each species/file

for v = 1:nFiles
    
    handles.vinput = fileNames(v);
    fileCharv = char(handles.vinput);
    fileCharv2 = char([pathName fileCharv]);
    vecMPS(v).File = fileCharv(1:end-4); %populates the file names in the structure
   
    if pathName==0, error('None selected!');
        else vecMPS(v).relVecMPS = reshape(conMPS(v).relMPS,1,[]); %just copies conMPS into vecMPS

    end  
end



% vecMPS = struct;
% %converts MPS rel values matrices into vectors for each species/file
% nFile = length(fileNames);
% 
% for v = 1:nFile
%     
%     handles.vinput = fileNames(v);
%     fileCharv = char(handles.vinput);
%     fileCharv2 = char([pathName fileCharv]);
%     vecMPS(v).File = fileCharv(1:end-4);
%     
%     if pathName==0, error('None selected!');
%         else vecMPS.relVecMPS = reshape(conMPS(v).relMPS,1,[]);
%     end  
% end


