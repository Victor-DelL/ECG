clear all; close all; clc
%Загрузка сигнала
fileID = fopen('Signal\S10.txt','r');
signal = fscanf(fileID,'%f ');
fclose(fileID);
fs = 250;% Выбор частоты дискртизации
%% График сигнала
%График сигнала
 figure('Name','Исходный ЭКГ сигнал','NumberTitle','off');
 plot(signal);
 title('Исходный ЭКГ сигнал');
 
%% CWT 
%Непрерывное вейвлет преобразование
Fa = 13;
waveletName = 'bior1.5';
waveletFrequency = centfrq(waveletName);
scale = round((waveletFrequency*fs)/Fa);
ecg_d = cwt(signal,scale,waveletName);
ecg_s = ecg_d.^2;
%% скользящее среднее
ecg_m = conv(ecg_s ,ones(1 ,round(0.150*fs))/round(0.150*fs));%conv-Свертка и умножение полиномов
%% Result
delay = round(0.15*fs/2);
delay_signal = zeros(numel(ecg_m)-delay,1);
for i=1:numel(ecg_m)-delay
 delay_signal(i) = ecg_m(i+delay-1);
end
%% поиск R пиков
peaks = [];
window_size = round(2.5*fs);
threshold = 0.5;
i = 1;
step = window_size*0.3;
while i<numel(delay_signal)
 window_start_point = i;
 if (i+window_size > numel(delay_signal)) % Обработка выхода за пределы
 window_start_point = numel(delay_signal)-window_size;
 step = 2*window_size;
 end
 
 % Поиск максимума в окне
 windowMax = max(delay_signal(window_start_point:window_start_point+window_size));
 % Определение порогового значения
 windowThreshold = threshold*windowMax;
 
 % поиск интервала, в котором лежит R пик
 maxInterval = []; % Интервал, в котором лежит R пик
 maxFlag = 0; 
 for j=round(window_start_point):(window_start_point+window_size)
    if delay_signal(j) > windowThreshold
        maxInterval(end+1) = j; % Увеличение интервала
        maxFlag = 1;
    else
        if (maxFlag == 1)
         % Поиск максимума в найденном интервале в исходном сигнале
        [value, location] = max(signal(maxInterval(1):(maxInterval(end)-round(numel(maxInterval)/3))));
        % Добавление "R"-пика
        peaks(end+1) = j - numel(maxInterval) + location; 
        % Очистка найденного интервала
        maxInterval = []; 
        end
        maxFlag = 0;
     end
 end
 
 i = i+step;
end
 
peaks = peaks - 1;
y1=get(gca,'ylim');

figure('Name','R зубцы','NumberTitle','off');
hold on;
plot(signal, 'b');
title('R зубцы');
 for j=1:numel(peaks)
    hold on,plot([peaks(j) peaks(j)], y1, 'r');
 end
hold off;

%% 
peaks = unique(peaks);
clearPeaks = [];
clearPeaks(end+1) = 1; % Временный пик
minDistance = 0.2*fs;
i = 1;
while i < numel(peaks)
 if(peaks(i)-clearPeaks(end) < minDistance)
 firstPeak = signal(peaks(i));
 secondPeak = signal(clearPeaks(end));
 if firstPeak > secondPeak
    clearPeaks(end) = [];
    clearPeaks(end+1) = peaks(i); 
 end
 else
    clearPeaks(end+1) = peaks(i);
 end
 i = i + 1;
end
clearPeaks(1) = [];
y1=get(gca,'ylim'); 
%% Save 
fileID = fopen('./Rpeaks.txt','w');
fprintf(fileID,'%d ', clearPeaks);
fclose(fileID);
%% Find P and T peaks
numel(signal)
% Вырезаем qRS комплексы (интерполяция)
PT_ecg = signal;
PT_threshold = 0.3;
PT_ecg = signal;
PT_threshold = 0.3;
for i=1:numel(clearPeaks)
 peakCWT_Value = delay_signal(clearPeaks(i));
 % find right bound
 k = clearPeaks(i);
 while delay_signal(k) > (peakCWT_Value*PT_threshold) && k > 1
 k = k - 1;
 end
 rightBound = k;
 % find left bound
 k = clearPeaks(i);
 while delay_signal(k) > (peakCWT_Value*PT_threshold) && k < numel(signal)
 k = k + 1;
 end
 leftBound = k;
 
 % Interpolation
 delta_A = signal(leftBound) - signal(rightBound);
 delta_x = leftBound - rightBound;
 step = delta_A/delta_x;
 current_value = signal(leftBound);
 for k=rightBound:leftBound
 PT_ecg(k) = current_value;
 current_value = current_value + step;
 end
end

Fa = 7;
waveletFrequency = centfrq(waveletName);
scale = round((waveletFrequency*fs)/Fa);
PT_ecg_d = cwt(PT_ecg,scale,waveletName); 
PT_ecg_s = PT_ecg_d.^2;
PT_ecg_m = conv(PT_ecg_s ,ones(1 ,round(0.15*fs))/round(0.15*fs));
delay = round(0.15*fs/2);
PT_delay_signal = zeros(numel(PT_ecg_m)-delay,1);
for i=1:numel(PT_ecg_m)-delay
 PT_delay_signal(i) = PT_ecg_m(i+delay-1);
end

%% Find PT peaks 
PT_peaks = [];
window_size = round(2.5*fs);
threshold = 0.2;
i = 1;
step = window_size*0.2;
while i<numel(PT_delay_signal)
 window_start_point = i;
 if (i+window_size > numel(PT_delay_signal))
 window_start_point = numel(PT_delay_signal)-window_size;
 step = 2*window_size;
 end

windowMax = max(PT_delay_signal(window_start_point:window_start_point+window_size));
 windowThreshold = threshold*windowMax;
 
 maxInterval = [];
 maxFlag = 0;
 for j=round(window_start_point):(window_start_point+window_size)
 if PT_delay_signal(j) > windowThreshold
 maxInterval(end+1) = j;
 maxFlag = 1;
 else
 if (maxFlag == 1) 
 [value, location] = max(signal(maxInterval(1):(maxInterval(end)-round(numel(maxInterval)/3))));
 PT_peaks(end+1) = j - numel(maxInterval) + location; 
 maxInterval = []; 
 end
 maxFlag = 0;
 end
 end
 
 i = i+step;
end

PT_peaks = PT_peaks - 1;
y1=get(gca,'ylim'); 

figure(11);
hold on;
 plot(signal, 'b');
 for j=1:numel(PT_peaks)
 hold on,plot([PT_peaks(j) PT_peaks(j)], y1, 'r');
 end
hold off;

%%
PT_peaks = unique(PT_peaks);
clearPT_Peaks = [];
clearPT_Peaks(end+1) = 1;
minDistance = 0.3*fs;
i = 1;
while i < numel(PT_peaks)
 if(PT_peaks(i)-clearPT_Peaks(end) < minDistance)
 firstPeak = signal(PT_peaks(i));
 secondPeak = signal(clearPT_Peaks(end));
 if firstPeak > secondPeak
 clearPT_Peaks(end) = [];
 clearPT_Peaks(end+1) = PT_peaks(i); 
 end
 else
 clearPT_Peaks(end+1) = PT_peaks(i);
 end
 i = i + 1;
end
clearPT_Peaks(1) = [];
y1=get(gca,'ylim');
figure('Name','Обнаруженные PT-зубцы','NumberTitle','off');
hold on;
 plot(signal,'g'); 
 plot(signal, 'b');
 title('Обнаруженные P и T зубцы');
 for j=1:numel(clearPT_Peaks)
 hold on,plot([clearPT_Peaks(j) clearPT_Peaks(j)], y1, 'r');
 end
hold off;

 %% Save
fileID = fopen('./PTpeaks.txt','w');
fprintf(fileID,'%d ', clearPT_Peaks);
fclose(fileID);
 


