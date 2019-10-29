%% Segement activity data
%%
issave = 1;
save_path = 'C:\Users\Laraw\OneDrive - University of Vermont\UVM\Research\McGinnis\Rescue_Climb\Data\Activity Segmentation';
%%
load('C:\Users\Laraw\OneDrive - University of Vermont\UVM\Research\McGinnis\Rescue_Climb\data\Preprocessed\EMTdata.mat')

locations ={'anterior_thigh_right','proximal_lateral_shank_left','dorsal_foot_left','proximal_lateral_shank_right','dorsal_foot_right','sacrum','anterior_thigh_left','medial_chest'}; 
subject = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12','S13','S14','S15','S16','S17','S18','S19','S20','S21','S22','S23','S24','S25'};
foot={'dorsal_foot_right','dorsal_foot_left'};
foot_short = {'R','L'};
bag={'HB','LB','HS','LS'};
task = {'WAT', 'RC'};
SegmentTimes_WAT = [];
SegmentTimes_RC = [];
Subject = {};
Load = [];
for i=[1,2,4:23]%subject
    gSall = data.(subject{i}).sacrum.gyro; % Rates of turn in sensor frame.   
    tSall = data.(subject{i}).sacrum.time/1000;
    gSHlall = data.(subject{i}).proximal_lateral_shank_left.gyro;
    tSHlall = data.(subject{i}).proximal_lateral_shank_left.time/1000;
    gSHrall = data.(subject{i}).proximal_lateral_shank_right.gyro;
    tSHrall = data.(subject{i}).proximal_lateral_shank_right.time/1000;
    startTimes = str2double(data.(subject{i}).annotations.StartTimestamp_ms_)/1000;
    endTimes = str2double(data.(subject{i}).annotations.StopTimestamp_ms_)/1000;
    AnnotInd_stand = find(strcmp(data.(subject{i}).annotations.EventType,'Blind Standing Balance'));
    AnnotInd_RC = find(strcmp(data.(subject{i}).annotations.EventType,'Rescue Climb'));
    AnnotInd_WAT = find(strcmp(data.(subject{i}).annotations.EventType,'Walk and Turn'));
    
    if length(startTimes)~=12
        fprintf('Annotation Error: %s',subject{i});
    end
    
    if length(AnnotInd_stand)~=4
        fprintf('Annotation Stand Error: %s',subject{i});
    end
    
    if length(AnnotInd_RC)~=4
        fprintf('Annotation RC Error: %s',subject{i});
    end
    
    for k = 1:4 %BagType
       for j=1:2 %Task: 1=walk and turn, 2 = rescue climb
            if j == 1
                %only look at BSB to WAT
                gS=gSall(tSall>=startTimes(AnnotInd_stand(k)) & tSall<=endTimes(AnnotInd_WAT(k)),:);%sacrum stand to end of recue climb
                tS=tSall(tSall>=startTimes(AnnotInd_stand(k)) & tSall<=endTimes(AnnotInd_WAT(k)));%sacrum
                gSHl=gSHlall(tSHlall>=startTimes(AnnotInd_stand(k)) & tSHlall<=endTimes(AnnotInd_WAT(k)),:);%shin
                tSHl=tSHlall(tSHlall>=startTimes(AnnotInd_stand(k)) & tSHlall<=endTimes(AnnotInd_WAT(k)));%shin
                gSHr=gSHrall(tSHrall>=startTimes(AnnotInd_stand(k)) & tSHrall<=endTimes(AnnotInd_WAT(k)),:);%shin
                tSHr=tSHrall(tSHrall>=startTimes(AnnotInd_stand(k)) & tSHrall<=endTimes(AnnotInd_WAT(k)));%shin

                %plot
                figure;
                plot(tS(tS>=startTimes(AnnotInd_WAT(k))),sqrt(sum(gS(tS>=startTimes(AnnotInd_WAT(k)),:).^2,2)),'r') %Sacrum gyro mag
                plot(tS(tS>=startTimes(AnnotInd_WAT(k))),gS(tS>=startTimes(AnnotInd_WAT(k)),2)) %Sacrum gyro yaw
                hold on
                plot(tSHl(tSHl>=startTimes(AnnotInd_WAT(k))),gSHl(tSHl>=startTimes(AnnotInd_WAT(k)),3))
                plot(tSHr(tSHr>=startTimes(AnnotInd_WAT(k))),-gSHr(tSHr>=startTimes(AnnotInd_WAT(k)),3))
                title(sprintf('%s %s sacrum gyro %s',subject{i},bag{k},'Walk and Turn'))
                axis tight
                turntimes_WAT = ginput(2);
                Subject = [Subject;subject{i}];
                Load = [Load; bag{k}];
                SegmentTimes_WAT = [SegmentTimes_WAT;startTimes(AnnotInd_WAT(k)),turntimes_WAT(:,1)',endTimes(AnnotInd_WAT(k))];
                      
% %             elseif j == 2
                gS=gSall(tSall>=startTimes(AnnotInd_stand(k)) & tSall<=endTimes(AnnotInd_RC(k)),:);%sacrum stand to end of recue climb
                tS=tSall(tSall>=startTimes(AnnotInd_stand(k)) & tSall<=endTimes(AnnotInd_RC(k)));%sacrum
                gSHl=gSHlall(tSHlall>=startTimes(AnnotInd_stand(k)) & tSHlall<=endTimes(AnnotInd_RC(k)),:);%shin
                tSHl=tSHlall(tSHlall>=startTimes(AnnotInd_stand(k)) & tSHlall<=endTimes(AnnotInd_RC(k)));%shin
                gSHr=gSHrall(tSHrall>=startTimes(AnnotInd_stand(k)) & tSHrall<=endTimes(AnnotInd_RC(k)),:);%shin
                tSHr=tSHrall(tSHrall>=startTimes(AnnotInd_stand(k)) & tSHrall<=endTimes(AnnotInd_RC(k)));%shin

                figure;
                plot(tS(tS>=startTimes(AnnotInd_RC(k))),gS(tS>=startTimes(AnnotInd_RC(k)),2))
                hold on
                plot(tSHl(tSHl>=startTimes(AnnotInd_RC(k))),gSHl(tSHl>=startTimes(AnnotInd_RC(k)),3))
                plot(tSHr(tSHr>=startTimes(AnnotInd_RC(k))),-gSHr(tSHr>=startTimes(AnnotInd_RC(k)),3))
                title(sprintf('%s %s sacrum gyro %s',subject{i},bag{k},'Rescue Climb'))
                axis tight
                turntimes_RC = ginput(6);
                SegmentTimes_RC = [SegmentTimes_RC;startTimes(AnnotInd_RC(k)),turntimes_RC(:,1)',endTimes(AnnotInd_RC(k))]; 
            end
       end
    end
end
 T = table(Subject,Load,SegmentTimes_WAT,SegmentTimes_RC);
 
if issave
    save(fullfile(save_path,'SegmentTimes.mat'),'T')
end