%% Import and setup data structure for thesis data

clear; clc;

subj = struct();

events = {'bsb','rc','wt'};
bag_type = {'hb','lb','hs','ls'};
sensors = {'anterior_thigh_right','anterior_thigh_left','dorsal_foot_left',...
    'dorsal_foot_right','medial_chest','proximal_lateral_shank_left',...
    'proximal_lateral_shank_right','sacrum'};
s_short = {'atr','atl','dfl','dfr','mdc','lsl','lsr','scm'};

%first check if the saved data file exists
if exist('..\Data\subj_data.mat','file')==2
   load('..\Data\subj_data.mat'); %if yes, load it
else %if not import all the data
    for i= [1,2,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]
        %base file name
        base = strcat('./thesis_data/asym_official_',num2str(i),'/ASYM_OFFICIAL/',num2str(i),'/');
        %annotation file path
        a_path = strcat(base,'annotations.csv');
        ann = readtable(a_path); %read the annotation file
        
        %setup data structure 
        for j=1:length(events)
            subj(i).(events{j}) = struct();
            for k=1:length(bag_type)
                subj(i).(events{j}).(bag_type{k}) = struct();
                for l=1:length(s_short)
                    subj(i).(events{j}).(bag_type{k}).(s_short{l}) = struct();
                end
            end
        end
        
        %loops to import data
        for j1=0:2 %3 different events
            ev = events{j1+1};

            for j2=1:3:12 %4 different bags, but need way to index 12 events
                bt = bag_type{ceil(j2/3)}; %back to 1:4
                j = j1+j2; %index for reading start,stop
                %pull start and stop values from annotation data
                start = str2double(ann{j,5});
                stop = str2double(ann{j,6});

                for k=1:length(sensors)
                    loc = s_short{k}; %sensor location in shorthand
                    fpath = strcat(base,sensors{k},'/'); %file path
                    a = importdata(strcat(fpath,'accel.csv')); %accel file path
                    g = importdata(strcat(fpath,'gyro.csv')); %gyro file path

                    %find the start and stop indices
                    [is,ie] = FindStartEnd(a.data(1:length(a.data),1),start,stop);
                    %extract relevant data
                    subj(i).(ev).(bt).(loc).a = a.data(is:ie,1:4);
                    [is,ie] = FindStartEnd(g.data(1:length(g.data)),start,stop);
                    subj(i).(ev).(bt).(loc).g = g.data(is:ie,1:4);
                end
            end
        end
    end
    save('subj_data.mat','subj');
end

    
function [is,ie] = FindStartEnd(t,start,stop)
[trash, is] = min(abs(start-t));
[trash, ie] = min(abs(stop-t));
end
    
    
    
    
    