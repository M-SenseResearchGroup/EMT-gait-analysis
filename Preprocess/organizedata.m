%% Create Struct
%%
issave = 1;
save_path = 'C:\Users\Laraw\OneDrive - University of Vermont\UVM\Research\McGinnis\Rescue_Climb\Data\Preprocessed';
%%
% Load data 
locations ={'anterior_thigh_right','proximal_lateral_shank_left','dorsal_foot_left','proximal_lateral_shank_right','dorsal_foot_right','sacrum','anterior_thigh_left','medial_chest'}; 
subject = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12','S13','S14','S15','S16','S17','S18','S19','S20','S21','S22','S23','S24','S25'};
d_path = 'C:\Users\Laraw\OneDrive - University of Vermont\UVM\Research\McGinnis\Rescue_Climb\Data\Raw Data\ASYM_OFFICIAL';
for i=1%[1,2,4:23]
    T = readtable(fullfile(d_path,sprintf('%d',i),'annotations.csv'));
    if size(T,1)<12
        fprintf('Error:Not Enough Annotations Subject %d\n',i)
    elseif size(T,1)>12
        if length(find(strcmp(T.EventType, 'Blind Standing Balance')))>4
            nums = find(strcmp(T.EventType, 'Blind Standing Balance'));
            for k = 2:length(nums)
                if nums(k)-1==nums(k-1)
                    T(nums(k-1),:)=[];
                end
            end
        end
        if length(find(strcmp(T.EventType, 'Walk and Turn')))>4
            nums = find(strcmp(T.EventType, 'Walk and Turn'));
            for k = 2:length(nums)
                if nums(k)-1==nums(k-1)
                    T(nums(k-1),:)=[];
                end
            end     
        end
        if length(find(strcmp(T.EventType, 'Rescue Climb')))>4
        nums = find(strcmp(T.EventType, 'Rescue Climb'));
            for k = 2:length(nums)
                if nums(k)-1==nums(k-1)
                    T(nums(k-1),:)=[];
                end
            end  
        end
        fprintf('Error:Annotations Subject %d\n',i)
    end
    fprintf('%s\n',subject{i})
    data.(subject{i}).annotations=T;
        
    for j=1:length(locations)
        fd1 = dir (fullfile(d_path,sprintf('%d',i),locations{j}));
        fd2 = dir (fullfile(d_path,sprintf('%d',i),locations{j},fd1(3).name));
        fn = fullfile(d_path,sprintf('%d',i),locations{j},fd1(3).name,fd2(3).name);
        
        data.(subject{i}).(locations{j}).accel = dlmread(fullfile(fn,'accel.csv'),',',1,1);
        data.(subject{i}).(locations{j}).time = dlmread(fullfile(fn,'accel.csv'),',',[1 0 0 0]);
        data.(subject{i}).(locations{j}).gyro = dlmread(fullfile(fn,'gyro.csv'),',',1,1);
    end
end

if issave
   save(fullfile(save_path,'EMTdataTest.mat'),'data') 
end

