function analysis(dateStrings1, dateStrings2)
tic
numDates = numel(dateStrings1);
% input example: dateStrings1 = ['20230511', '20230518', etc], dateStrings2 = ['2023-05-11', '2023-05-18', etc] 
hit_eventPhase_cell = {};
miss_eventPhase_cell = {};
falsealarm_eventPhase_cell = {};
correctrejection_eventPhase_cell = {};
eventPhase_cell = {};
not_eventPhase_cell = {};
all_hit_eventMeanPhases = [];
all_miss_eventMeanPhases = [];
all_correctrejection_eventMeanPhases = [];
all_falsealarm_eventMeanPhases = [];
all_eventMeanPhases = [];
all_not_eventMeanPhases = [];
all_circ_otest1_results = [];
all_circ_otest2_results = [];
all_circ_otest3_results = [];
all_circ_otest4_results = [];
all_circ_otest5_results = [];
all_circ_otest6_results = [];
all_hit_eventMeanPhases_deg = [];
all_miss_eventMeanPhases_deg = [];
all_correctrejection_eventMeanPhases_deg = [];
all_falsealarm_eventMeanPhases_deg = [];
all_eventMeanPhases_deg = [];
all_not_eventMeanPhases_deg = [];
eventPhase_alldates = [];
not_eventPhase_alldates = [];
hit_eventPhase_alldates = [];
miss_eventPhase_alldates = [];
falsealarm_eventPhase_alldates = [];
correctrejection_eventPhase_alldates = [];
% loading files
a = 'Y:\Projects\Pulv_bodysignal\Magnus_SDT\';
    for num = 1:numDates 
dateString1 = dateStrings1{num};
ecgFilename = [a dateString1 '_ecg.mat'];
ecg = load(ecgFilename);
dateString2 = dateStrings2{num};
filePattern = [a 'Magcombined' dateString2 '*.mat'];
fileList = dir(filePattern);
numFiles = numel(fileList);
Magcombined = cell(1, numFiles);
for i = 1:numFiles
    names{i} = fileList(i).name;
    loadedData = load([fileList(i).name]);
    Magcombined{i} = loadedData;
    blockId = arrayfun(@(x) str2double(x.name(32:33)), fileList);
end
% storing trials of each type in cells and converting to structures and arrays
for i = 1:numFiles
    numTrials = numel(Magcombined{i}.trial);
    idx_hit_trials{i} = false(1, numTrials);
    idx_miss_trials{i} = false(1, numTrials);
    idx_falsealarm_trials{i} = false(1, numTrials);
    idx_correctrejection_trials{i} = false(1, numTrials);
    for h = 1:numTrials
        trialType = Magcombined{i}.trial(h).SDT_trial_type;
        idx_hit_trials{i}(h) = trialType == 1;
        idx_miss_trials{i}(h) = trialType == 2;
        idx_falsealarm_trials{i}(h) = trialType == 3;
        idx_correctrejection_trials{i}(h) = trialType == 4;
    end
    hit_trials = Magcombined{i}.trial(idx_hit_trials{i});
    hit_trials_cell{i} = Magcombined{i}.trial(idx_hit_trials{i});
    for jh = 1:numel(hit_trials)
    hit_times{i}(jh) = hit_trials(jh).TDT_state_onsets_aligned_to_1st_INI(4,1);
    end  
    miss_trials = Magcombined{i}.trial(idx_miss_trials{i});
    miss_trials_cell{i} = Magcombined{i}.trial(idx_miss_trials{i});
    for jm = 1:numel(miss_trials)
    miss_times{i}(jm) = miss_trials(jm).TDT_state_onsets_aligned_to_1st_INI(4,1);
    end 
    falsealarm_trials = Magcombined{i}.trial(idx_falsealarm_trials{i});
    falsealarm_trials_cell{i} = Magcombined{i}.trial(idx_falsealarm_trials{i});
    for jf = 1:numel(falsealarm_trials)   
    falsealarm_times{i}(jf) = falsealarm_trials(jf).TDT_state_onsets_aligned_to_1st_INI(4,1);
    end 
    correctrejection_trials = Magcombined{i}.trial(idx_correctrejection_trials{i});
    correctrejection_trials_cell{i} = Magcombined{i}.trial(idx_correctrejection_trials{i});
    for jc = 1:numel(correctrejection_trials)
    correctrejection_times{i}(jc) = correctrejection_trials(jc).TDT_state_onsets_aligned_to_1st_INI(4,1);
    end 
    idx_completed_trials = [Magcombined{i}.trial.completed];
    valid_indices2 = idx_completed_trials >= 1 & idx_completed_trials <= numTrials;
    completed_trials{i} = Magcombined{i}.trial(valid_indices2);
    completed_trials1 = Magcombined{i}.trial(valid_indices2);
    for jcom = 1:numel(completed_trials1)
    completed_times{i}(jcom) = completed_trials1(jcom).TDT_state_onsets_aligned_to_1st_INI(4,1);
    end 
        idx_noncompleted_trials = ~[Magcombined{i}.trial.completed];
        noncompleted_trials1 = Magcombined{i}.trial(idx_noncompleted_trials);
        noncompleted_trials{i} = Magcombined{i}.trial(idx_noncompleted_trials);
        structures_with_4s = {};
        
        for idx1 = 1:numel(noncompleted_trials1)
            if any(noncompleted_trials1(idx1).TDT_states == 4)
                structures_with_4s = noncompleted_trials1;
                for jnonn = 1:numel(structures_with_4s)
                not_event_times{i}(jnonn) = structures_with_4s(jnonn).TDT_state_onsets_aligned_to_1st_INI(4,1);
                end
            end
        end
end

for l = 1:size(blockId,1)
    n(l) = blockId(l);
    idx = n(l);
eventTimesNorm = cell(size(blockId, 1), 1);
not_eventTimesNorm = cell(size(blockId, 1), 1);
hit_eventTimesNorm = cell(size(blockId, 1), 1);
miss_eventTimesNorm = cell(size(blockId, 1), 1);
falsealarm_eventTimesNorm = cell(size(blockId, 1), 1);
correctrejection_eventTimesNorm = cell(size(blockId, 1), 1);
ecgCycleDurations = [ecg.out(idx).R2R_t(ecg.out(idx).idx_valid_R2R_consec)];
R2R_t = [ecg.out(idx).R2R_t(ecg.out(idx).idx_valid_R2R_consec)];
R2R_valid = [ecg.out(idx).R2R_valid(ecg.out(idx).idx_valid_R2R_consec)];
numCycles = length(ecgCycleDurations);
for m = 1:numCycles
    cycleStart = R2R_t(m)- R2R_valid(m);
    cycleEnd = R2R_t(m);
    cycleDuration = cycleEnd - cycleStart;
    eventTimesCycle = completed_times{1,l}((completed_times{1,l} >= cycleStart) & (completed_times{1,l} < cycleEnd));
    eventTimesNorm{l}((completed_times{1,l} >= cycleStart) & (completed_times{1,l} < cycleEnd)) = (eventTimesCycle - cycleStart) / cycleDuration;
    not_eventTimesCycle = not_event_times{l}((not_event_times{l} >= cycleStart) & (not_event_times{l} < cycleEnd));
    not_eventTimesNorm{l}((not_event_times{l} >= cycleStart) & (not_event_times{l} < cycleEnd)) = (not_eventTimesCycle - cycleStart) / cycleDuration;
    hit_eventTimesCycle =  hit_times{l}(( hit_times{l} >= cycleStart) & ( hit_times{l} < cycleEnd));
    hit_eventTimesNorm{l}((hit_times{l} >= cycleStart) & (hit_times{l} < cycleEnd)) = (hit_eventTimesCycle - cycleStart) / cycleDuration;
    miss_eventTimesCycle = miss_times{l}((miss_times{l} >= cycleStart) & (miss_times{l} < cycleEnd));
    miss_eventTimesNorm{l}((miss_times{l} >= cycleStart) & (miss_times{l} < cycleEnd)) = (miss_eventTimesCycle - cycleStart) / cycleDuration;
    falsealarm_eventTimesCycle = falsealarm_times{l}((falsealarm_times{l} >= cycleStart) & (falsealarm_times{l} < cycleEnd));
    falsealarm_eventTimesNorm{l}((falsealarm_times{l} >= cycleStart) & (falsealarm_times{l} < cycleEnd)) = (falsealarm_eventTimesCycle - cycleStart) / cycleDuration;
    correctrejection_eventTimesCycle = correctrejection_times{l}((correctrejection_times{l} >= cycleStart) & (correctrejection_times{l} < cycleEnd));
    correctrejection_eventTimesNorm{l}((correctrejection_times{l} >= cycleStart) & (correctrejection_times{l} < cycleEnd)) = (correctrejection_eventTimesCycle - cycleStart) / cycleDuration;

end
hit_eventPhase{l} = 2*pi*hit_eventTimesNorm{l};
miss_eventPhase{l} = 2*pi*miss_eventTimesNorm{l};
falsealarm_eventPhase{l} = 2*pi*falsealarm_eventTimesNorm{l};
correctrejection_eventPhase{l} = 2*pi*correctrejection_eventTimesNorm{l};
eventPhase{l} = 2*pi*eventTimesNorm{l};
not_eventPhase{l} = 2*pi*not_eventTimesNorm{l};

eventPhase_all = horzcat(eventPhase{:});
not_eventPhase_all = horzcat(not_eventPhase{:});   
hit_eventPhase_all = horzcat(hit_eventPhase{:});
miss_eventPhase_all = horzcat(miss_eventPhase{:});
falsealarm_eventPhase_all = horzcat(falsealarm_eventPhase{:});
correctrejection_eventPhase_all = horzcat(correctrejection_eventPhase{:});


hit_eventMeanPhase_session = circ_mean(hit_eventPhase_all');
miss_eventMeanPhase_session = circ_mean(miss_eventPhase_all');
falsealarm_eventMeanPhase_session = circ_mean(falsealarm_eventPhase_all');
correctrejection_eventMeanPhase_session = circ_mean(correctrejection_eventPhase_all');
eventMeanPhase_session = circ_mean(eventPhase_all');
not_eventMeanPhase_session = circ_mean(not_eventPhase_all');

hit_eventMeanPhase_session(hit_eventMeanPhase_session < 0) = hit_eventMeanPhase_session(hit_eventMeanPhase_session < 0) + 2*pi;
miss_eventMeanPhase_session(miss_eventMeanPhase_session < 0) = miss_eventMeanPhase_session(miss_eventMeanPhase_session < 0) + 2*pi;
falsealarm_eventMeanPhase_session(falsealarm_eventMeanPhase_session < 0) = falsealarm_eventMeanPhase_session(falsealarm_eventMeanPhase_session < 0) + 2*pi;
correctrejection_eventMeanPhase_session(correctrejection_eventMeanPhase_session < 0) = correctrejection_eventMeanPhase_session(correctrejection_eventMeanPhase_session < 0) + 2*pi;
eventMeanPhase_session(eventMeanPhase_session < 0) = eventMeanPhase_session(eventMeanPhase_session < 0) + 2*pi;
not_eventMeanPhase_session(not_eventMeanPhase_session < 0) = not_eventMeanPhase_session(not_eventMeanPhase_session < 0) + 2*pi;

hit_eventMeanPhase_deg_session = hit_eventMeanPhase_session * (180/pi);
miss_eventMeanPhase_deg_session = miss_eventMeanPhase_session * (180/pi);
falsealarm_eventMeanPhase_deg_session = falsealarm_eventMeanPhase_session * (180/pi);
correctrejection_eventMeanPhase_deg_session = correctrejection_eventMeanPhase_session * (180/pi);
eventMeanPhase_deg_session = eventMeanPhase_session * (180/pi);
not_eventMeanPhase_deg_session = not_eventMeanPhase_session * (180/pi);

end

eventPhase_alldates{num} = eventPhase_all;
not_eventPhase_alldates{num} = not_eventPhase_all;
hit_eventPhase_alldates{num} = hit_eventPhase_all;
miss_eventPhase_alldates{num} = miss_eventPhase_all;
falsealarm_eventPhase_alldates{num} = falsealarm_eventPhase_all;
correctrejection_eventPhase_alldates{num} = correctrejection_eventPhase_all;

all_hit_eventMeanPhases = [all_hit_eventMeanPhases; hit_eventMeanPhase_session];
all_miss_eventMeanPhases = [all_miss_eventMeanPhases; miss_eventMeanPhase_session];
all_correctrejection_eventMeanPhases = [all_correctrejection_eventMeanPhases; correctrejection_eventMeanPhase_session];
all_falsealarm_eventMeanPhases = [all_falsealarm_eventMeanPhases; falsealarm_eventMeanPhase_session];
all_eventMeanPhases = [all_eventMeanPhases; eventMeanPhase_session];
all_not_eventMeanPhases = [all_not_eventMeanPhases; not_eventMeanPhase_session];
    
all_hit_eventMeanPhases_deg = [all_hit_eventMeanPhases_deg; hit_eventMeanPhase_deg_session];
all_miss_eventMeanPhases_deg = [all_miss_eventMeanPhases_deg; miss_eventMeanPhase_deg_session];
all_correctrejection_eventMeanPhases_deg = [all_correctrejection_eventMeanPhases_deg; correctrejection_eventMeanPhase_deg_session];
all_falsealarm_eventMeanPhases_deg = [all_falsealarm_eventMeanPhases_deg; falsealarm_eventMeanPhase_deg_session];
all_eventMeanPhases_deg = [all_eventMeanPhases_deg; eventMeanPhase_deg_session];
all_not_eventMeanPhases_deg = [all_not_eventMeanPhases_deg; not_eventMeanPhase_deg_session];
        end    
circ_otest1_sessions = circ_rtest(all_hit_eventMeanPhases);
circ_otest2_sessions = circ_rtest(all_miss_eventMeanPhases);
circ_otest3_sessions = circ_rtest(all_correctrejection_eventMeanPhases);
circ_otest4_sessions = circ_rtest(all_falsealarm_eventMeanPhases);
circ_otest5_sessions = circ_rtest(all_eventMeanPhases);
circ_otest6_sessions = circ_rtest(all_not_eventMeanPhases); 


figure;
subplot(2, 3, 1);
if ~isempty(all_hit_eventMeanPhases)
    polar(all_hit_eventMeanPhases, ones(size(all_hit_eventMeanPhases)), 'ro'); 
    title('Hit');
    hold on
axis equal;
hold off
end
subplot(2, 3, 2);
if ~isempty(all_miss_eventMeanPhases)
polar(all_miss_eventMeanPhases, ones(size(all_miss_eventMeanPhases)), 'ro'); 
    title('Miss');
    hold on
axis equal;
hold off
end

subplot(2, 3, 3);
if ~isempty(all_correctrejection_eventMeanPhases)
polar(all_correctrejection_eventMeanPhases, ones(size(all_correctrejection_eventMeanPhases)), 'ro'); 
    title('Correct rejection');
    hold on
axis equal;
hold off
end

subplot(2, 3, 4);
if ~isempty(all_falsealarm_eventMeanPhases)
polar(all_falsealarm_eventMeanPhases, ones(size(all_falsealarm_eventMeanPhases)), 'ro'); 
    title('False alarms');
    hold on
axis equal;
hold off
end

subplot(2, 3, 5);
if ~isempty(all_eventMeanPhases)
polar(all_eventMeanPhases, ones(size(all_eventMeanPhases)), 'ro'); 
    title('Completed');
    hold on
axis equal;
hold off
end

subplot(2, 3, 6);
if ~isempty(all_not_eventMeanPhases)
polar(all_not_eventMeanPhases, ones(size(all_not_eventMeanPhases)), 'ro'); 
    title('Not completed');
    hold on
    axis equal;
    hold off
end

 total = sum(all_hit_eventMeanPhases_deg);
    mean1 = total / length(all_hit_eventMeanPhases_deg);
disp('Hit Event Mean Phases for all sessions (dates):');
disp(mean1);

 total = sum(all_miss_eventMeanPhases_deg);
    mean2 = total / length(all_miss_eventMeanPhases_deg);
disp('Miss Event Mean Phases for all sessions:');
disp(mean2);

 total = sum(all_correctrejection_eventMeanPhases_deg);
    mean3 = total / length(all_correctrejection_eventMeanPhases_deg);
disp('Correct rejections Event Mean Phases for all sessions:');
disp(mean3);

total = sum(all_falsealarm_eventMeanPhases_deg);
    mean4 = total / length(all_falsealarm_eventMeanPhases_deg);
disp('False alarm Event Mean Phases for all sessions:');
disp(mean4);

total = sum(all_eventMeanPhases_deg);
    mean6 = total / length(all_eventMeanPhases_deg);
disp('All Completed Event Mean Phases for sessions:');
disp(mean6);

total = sum(all_not_eventMeanPhases_deg);
    mean7 = total / length(all_not_eventMeanPhases_deg);
disp('All Not Completed Event Mean Phases for sessions:');
disp(mean7);

disp('Hit R-tests for sessions:');
disp(circ_otest1_sessions);

disp('All Miss R-tests for sessions:');
disp(circ_otest2_sessions);

disp('All Correct Rejection R-tests for sessions:');
disp(circ_otest3_sessions);

disp('All False alarms R-tests for sessions:');
disp(circ_otest4_sessions);

disp('All Completed R-tests for sessions:');
disp(circ_otest5_sessions);

disp('All Not Completed R-tests for sessions:');
disp(circ_otest6_sessions);

eventPhase_alldates1 = horzcat(eventPhase_alldates{:});
not_eventPhase_alldates1 = horzcat(not_eventPhase_alldates{:});   
hit_eventPhase_alldates1 = horzcat(hit_eventPhase_alldates{:});
miss_eventPhase_alldates1 = horzcat(miss_eventPhase_alldates{:});
falsealarm_eventPhase_alldates1 = horzcat(falsealarm_eventPhase_alldates{:});
correctrejection_eventPhase_alldates1 = horzcat(correctrejection_eventPhase_alldates{:});

hit_eventMeanPhase = circ_mean(hit_eventPhase_alldates1');
hit_mean_length = circ_r(hit_eventPhase_alldates1')*100;
hit_x_end = hit_mean_length * cos(hit_eventMeanPhase);
hit_y_end = hit_mean_length * sin(hit_eventMeanPhase);

miss_eventMeanPhase = circ_mean(miss_eventPhase_alldates1');
miss_mean_length = circ_r(miss_eventPhase_alldates1')*100;
miss_x_end = miss_mean_length * cos(miss_eventMeanPhase);
miss_y_end = miss_mean_length * sin(miss_eventMeanPhase);

falsealarm_eventMeanPhase = circ_mean(falsealarm_eventPhase_alldates1');
falsealarm_mean_length = circ_r(falsealarm_eventPhase_alldates1')*100;
falsealarm_x_end = falsealarm_mean_length * cos(falsealarm_eventMeanPhase);
falsealarm_y_end = falsealarm_mean_length * sin(falsealarm_eventMeanPhase);

correctrejection_eventMeanPhase = circ_mean(correctrejection_eventPhase_alldates1');
correctrejection_mean_length = circ_r(correctrejection_eventPhase_alldates1')*100;
correctrejection_x_end = correctrejection_mean_length * cos(correctrejection_eventMeanPhase);
correctrejection_y_end = correctrejection_mean_length * sin(correctrejection_eventMeanPhase);

eventMeanPhase = circ_mean(eventPhase_alldates1');
mean_length = circ_r(eventPhase_alldates1')*100;
x_end = mean_length * cos(eventMeanPhase);
y_end = mean_length * sin(eventMeanPhase);

not_eventMeanPhase = circ_mean(not_eventPhase_alldates1');
not_mean_length = circ_r(not_eventPhase_alldates1')*100;
not_x_end = not_mean_length * cos(not_eventMeanPhase);
not_y_end = not_mean_length * sin(not_eventMeanPhase);

hit_eventMeanPhase(hit_eventMeanPhase < 0) = hit_eventMeanPhase(hit_eventMeanPhase < 0) + 2*pi;
miss_eventMeanPhase(miss_eventMeanPhase < 0) = miss_eventMeanPhase(miss_eventMeanPhase < 0) + 2*pi;
falsealarm_eventMeanPhase(falsealarm_eventMeanPhase < 0) = falsealarm_eventMeanPhase(falsealarm_eventMeanPhase < 0) + 2*pi;
correctrejection_eventMeanPhase(correctrejection_eventMeanPhase < 0) = correctrejection_eventMeanPhase(correctrejection_eventMeanPhase < 0) + 2*pi;
eventMeanPhase(eventMeanPhase < 0) = eventMeanPhase(eventMeanPhase < 0) + 2*pi;
not_eventMeanPhase(not_eventMeanPhase < 0) = not_eventMeanPhase(not_eventMeanPhase < 0) + 2*pi;

hit_total_count = numel(hit_eventPhase_alldates1);
miss_total_count = numel(miss_eventPhase_alldates1);
falsealarm_total_count = numel(falsealarm_eventPhase_alldates1);
correctrejection_total_count = numel(correctrejection_eventPhase_alldates1);
total_count = numel(eventPhase_alldates1);
not_total_count = numel(not_eventPhase_alldates1);

figure;
subplot(2, 2, 1);
if ~isempty(hit_eventPhase_alldates1)
    ig_rose(hit_eventPhase_alldates1, 20, true); 
    title('Hit');
    hold on
text(-1.5, 1.5, ['Total: ', num2str(hit_total_count)], 'Color', 'black', 'FontWeight', 'bold');
plot([0, hit_x_end], [0, hit_y_end], 'r', 'LineWidth', 2);
axis equal;
    hold off;
end
subplot(2, 2, 2);
if ~isempty(miss_eventPhase_alldates1)
    ig_rose(miss_eventPhase_alldates1, 20, true);
    title('Miss');
    hold on;
    plot([0, miss_x_end], [0, miss_y_end], 'r', 'LineWidth', 2);
    text(-1.5, 1.5, ['Total: ', num2str(miss_total_count)], 'Color', 'black', 'FontWeight', 'bold');
    axis equal;
    hold off;
end

subplot(2, 2, 3);
if ~isempty(falsealarm_eventPhase_alldates1)
    ig_rose(falsealarm_eventPhase_alldates1, 20, true);
    title('False Alarm');
    hold on;
    plot([0, falsealarm_x_end], [0, falsealarm_y_end], 'r', 'LineWidth', 2);
    text(-1.5, 1.5, ['Total: ', num2str(falsealarm_total_count)], 'Color', 'black', 'FontWeight', 'bold');
    axis equal;
    hold off;
end

subplot(2, 2, 4);
if ~isempty(correctrejection_eventPhase_alldates1)
    ig_rose(correctrejection_eventPhase_alldates1, 20, true);
    title('Correct Rejection');
    hold on;
    plot([0, correctrejection_x_end], [0, correctrejection_y_end], 'r', 'LineWidth', 2);
    text(-1.5, 1.5, ['Total: ', num2str(correctrejection_total_count)], 'Color', 'black', 'FontWeight', 'bold');
    axis equal;
    hold off;
end
spacing = 1;

figure;
subplot(1, 2, 1);
if ~isempty(eventPhase_alldates1)
    ig_rose(eventPhase_alldates1, 20, true); 
    title('completed');
    hold on
text(-1.5, 1.5, ['Total: ', num2str(total_count)], 'Color', 'black', 'FontWeight', 'bold');
plot([0, x_end], [0, y_end], 'r', 'LineWidth', 2);
axis equal;
    hold off;
end
subplot(1, 2, 2);
if ~isempty(not_eventPhase_alldates1)
    ig_rose(not_eventPhase_alldates1, 20, true);
    title('Not completed');
    hold on;
    plot([0, not_x_end], [0, not_y_end], 'r', 'LineWidth', 2);
    text(-1.5, 1.5, ['Total: ', num2str(not_total_count)], 'Color', 'black', 'FontWeight', 'bold');
    axis equal;
    hold off;
end

hit_eventMeanPhase_deg = hit_eventMeanPhase * (180/pi);
miss_eventMeanPhase_deg = miss_eventMeanPhase * (180/pi);
falsealarm_eventMeanPhase_deg = falsealarm_eventMeanPhase * (180/pi);
correctrejection_eventMeanPhase_deg = correctrejection_eventMeanPhase * (180/pi);
eventMeanPhase_deg = eventMeanPhase * (180/pi);
not_eventMeanPhase_deg = not_eventMeanPhase * (180/pi);

disp('Circular Mean Phase Angles (degrees):');
disp(['Completed: ', num2str(eventMeanPhase_deg)]);
disp(['Not completed: ', num2str(not_eventMeanPhase_deg)]);
disp(['Hits: ', num2str(hit_eventMeanPhase_deg)]);
disp(['Misses: ', num2str(miss_eventMeanPhase_deg)]);
disp(['False Alarms: ', num2str(falsealarm_eventMeanPhase_deg)]);
disp(['Correct Rejections: ', num2str(correctrejection_eventMeanPhase_deg)]);
circ_otest1 = circ_otest(hit_eventPhase_alldates1,1);
circ_otest2 = circ_otest(miss_eventPhase_alldates1,1);
circ_otest3 = circ_otest(correctrejection_eventPhase_alldates1,1);
circ_otest4 = circ_otest(falsealarm_eventPhase_alldates1,1);
circ_otest5 = circ_otest(eventPhase_alldates1,1);
circ_otest6 = circ_otest(not_eventPhase_alldates1,1);
disp('O-test p-values counted for all trials together:');
disp(['Completed: ', num2str(circ_otest5)]);
disp(['Not completed: ', num2str(circ_otest6)]);
disp(['Hits: ', num2str(circ_otest1)]);
disp(['Misses: ', num2str(circ_otest2)]);
disp(['False Alarms: ', num2str(circ_otest3)]);
disp(['Correct Rejections: ', num2str(circ_otest4)]);

disp('Rayleigh Test p-value for all trials together:');
[h, p] = circ_rtest(hit_eventPhase_alldates1);
fprintf('Hits: %.4f\n', p);
[h1, p1] = circ_rtest(eventPhase_alldates1);
fprintf('Completed: %.4f\n', h1);
[h2, p2] = circ_rtest(miss_eventPhase_alldates1);
fprintf('Misses: %.4f\n', h2);
[h3, p3] = circ_rtest(not_eventPhase_alldates1);
fprintf('Not completed: %.4f\n', h3);
[h4, p4] = circ_rtest(correctrejection_eventPhase_alldates1);
fprintf('Correct rejections: %.4f\n', h4);
[h5, p5] = circ_rtest(falsealarm_eventPhase_alldates1);
fprintf('False alarms: %.4f\n', h5);
total_trials = numel(eventPhase_alldates1);
fprintf('Number of completed trials: %d\n', total_trials);
pHit = numel(hit_eventPhase_alldates1) / total_trials;
pFA = numel(falsealarm_eventPhase_alldates1) / total_trials;
testsim_dprime(pHit, pFA);
  toc
  

