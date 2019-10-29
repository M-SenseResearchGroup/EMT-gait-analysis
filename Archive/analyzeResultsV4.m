%% AnalyzeResultsV4
% corresponds to stairclimbv4
% Lara Weed
% 6 Dec 2019
%%
%load tempStep
%load tempStepBig
%% Create indices
A=stepT.Type_Activity=='A';
D=stepT.Type_Activity=='D';
Landing=stepT.Type_Activity=='L';
R=stepT.fside=='R';
L=stepT.fside=='L';
subject = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12','S13','S14','S15','S16','S17','S18','S19','S20','S21','S22','S23','S24','S25'};
%%
for i=[1,2,4:23]
    sub_ind = strcmp(stepT.Subject,subject(:,i));
    figure;
    p1 = subplot(4,1,1);
    title(subject(:,i))
    plot(stepT.start_time(sub_ind & A),stepT.step_length(sub_ind & A),'r*')
    hold on
    plot(stepT.start_time(sub_ind & D),stepT.step_length(sub_ind & D),'b*')
    hold on
    plot(stepT.start_time(sub_ind & Landing),stepT.step_length(sub_ind & Landing),'g*')
  
    p2 = subplot(4,1,2);
    plot(stepT.start_time(sub_ind & A),stepT.lateral_dev(sub_ind & A),'r*')
    hold on
    plot(stepT.start_time(sub_ind & D),stepT.lateral_dev(sub_ind & D),'b*')
    hold on
    plot(stepT.start_time(sub_ind & Landing),stepT.lateral_dev(sub_ind & Landing),'g*')
  
    p3 = subplot(4,1,3);
    plot(stepT.start_time(sub_ind & A),stepT.step_height(sub_ind & A),'r*')
    hold on
    plot(stepT.start_time(sub_ind & D),stepT.step_height(sub_ind & D),'b*')
    hold on
    plot(stepT.start_time(sub_ind & Landing),stepT.step_height(sub_ind & Landing),'g*')
    
    p4 = subplot(4,1,4);
    plot(stepT.start_time(sub_ind & A),stepT.max_swing_vel(sub_ind & A),'r*')
    hold on
    plot(stepT.start_time(sub_ind & D),stepT.max_swing_vel(sub_ind & D),'b*')
    hold on
    plot(stepT.start_time(sub_ind & Landing),stepT.max_swing_vel(sub_ind & Landing),'g*')
    linkaxes([p1 p2 p3 p4],'x')
    
    figure;
    p5 = subplot(4,1,1);
    title(subject(:,i))
    plot(stepT.start_time(sub_ind & A),stepT.foot_attack_angle(sub_ind & A),'r*')
    hold on
    plot(stepT.start_time(sub_ind & D),stepT.foot_attack_angle(sub_ind & D),'b*')
    hold on
    plot(stepT.start_time(sub_ind & Landing),stepT.foot_attack_angle(sub_ind & Landing),'g*')
    
    p6 = subplot(4,1,2);
    plot(stepT.start_time(sub_ind & A),stepT.contact_time(sub_ind & A),'r*')
    hold on
    plot(stepT.start_time(sub_ind & D),stepT.contact_time(sub_ind & D),'b*')
    hold on
    plot(stepT.start_time(sub_ind & Landing),stepT.contact_time(sub_ind & Landing),'g*')
    
    p7 = subplot(4,1,3);
    plot(stepT.start_time(sub_ind & A),stepT.step_time(sub_ind & A),'r*')
    hold on
    plot(stepT.start_time(sub_ind & D),stepT.step_time(sub_ind & D),'b*')
    hold on
    plot(stepT.start_time(sub_ind & Landing),stepT.step_time(sub_ind & Landing),'g*')
  
    p8 = subplot(4,1,4);
    plot(stepT.start_time(sub_ind & A),stepT.cadence(sub_ind & A),'r*')
    hold on
    plot(stepT.start_time(sub_ind & D),stepT.cadence(sub_ind & D),'b*')
    hold on
    plot(stepT.start_time(sub_ind & Landing),stepT.cadence(sub_ind & Landing),'g*')
    linkaxes([p5 p6 p7 p8],'x')
end



