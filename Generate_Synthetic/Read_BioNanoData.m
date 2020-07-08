clc; clear all
%% In this script we extract the information for each cutter 
% The first cutter
fid = fopen('results_Nb.BbvCI.txt','r');
i= 1;
tline = fgetl(fid);
dictionary_genome = cell(1);
A{1} = tline;
dictionary_genome{1,1} = A{1};
cnt = 2;
frag_len_1 = [];
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
    if ~isempty(strfind(A{i}, '>chr'))
        dictionary_genome{1,cnt} = A{i};
        frag_len_1 = [frag_len_1,0];
        cnt = cnt +1;
    elseif ischar(A{i})
        if str2num(A{i}) > 0
            frag_len_1 = [frag_len_1,str2num(A{i})];
        end
    end
    if i/10000 == floor(i/10000)
        i
    end
end
fclose(fid);
clear A

% the second cutter
fid = fopen('results_NbBssSI.txt','r');
i= 1;
tline = fgetl(fid);
A{1} = tline;
dictionary_genome{2,1} = A{1};
cnt = 2;
frag_len_2 = [];
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
    if ~isempty(strfind(A{i}, '>chr'))
        dictionary_genome{2,cnt} = A{i};
        frag_len_2 = [frag_len_2,0];
        cnt = cnt +1;
    elseif ischar(A{i})
        if str2num(A{i}) > 0
            frag_len_2 = [frag_len_2,str2num(A{i})];
        end
    end
    if i/10000 == floor(i/10000)
        i
    end
end
fclose(fid);
clear A

% the third cutter
fid = fopen('results_nt.BspQI.txt','r');
i= 1;
tline = fgetl(fid);
A{1} = tline;
dictionary_genome{3,1} = A{1};
cnt = 2;
frag_len_3 = [];
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
    if ~isempty(strfind(A{i}, '>chr'))
        dictionary_genome{3,cnt} = A{i};
        frag_len_3 = [frag_len_3,0];
        cnt = cnt +1;
    elseif ischar(A{i}) 
        if str2num(A{i}) > 0
            frag_len_3 = [frag_len_3,str2num(A{i})];
        end
    end
    if i/10000 == floor(i/10000)
        i
    end
end
fclose(fid);

cd ../Synthetic_Data
save('fragment_len_1','frag_len_1')
save('fragment_len_2','frag_len_2')
save('fragment_len_3','frag_len_3')
save('dictionary_genome','dictionary_genome')
cd .. 