clc; clear all 
addpath readReal

phyLength = [];
numReads = 1000;

for i = 1:numReads
    readFileName = strcat('fragment_', num2str(i - 1), '.txt');
    fid = fopen(readFileName);
    original_data = textscan(fid, '%s', 'delimiter', '\n');
    tempSum = 0;
    for k = 1:numel(original_data{1, 1})
        tempSum = tempSum + str2num(original_data{1, 1}{k});
    end
    phyLength = [phyLength; tempSum];
    fclose(fid)
end
save('phyLength', 'phyLength')