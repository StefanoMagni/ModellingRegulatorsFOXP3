% Script to compute optimal time-point set up for Feng's experiment
% Created on January the 11th, 2018
% Authors: stefano.magni@uni.lu, rucha.sawleker@uni.lu

close all; clear; clc;

MinTimePoint = 0;%100;
MaxTimePoint = 600;%420;
TimePointStep = 2;%20;

FOXP3WindowWIDTH = 80;

MinFOXP3WindowSTART = 100;
StepFOXP3WindowSTART = 2;%20;
MaxFOXP3WindowSTART = 320;

workingDir = pwd;

%%%%%%%%% 1 MEASUREMENT %%%%%%%%%

ProbabilitiesArrayMissing = [];
ProbabilitiesHittingArray = [];
TimesArray = [];

for MeasurementTime1 = MinTimePoint:TimePointStep:MaxTimePoint
    TimesArray = [TimesArray, MeasurementTime1];
    ProbabilityOfMissingCurrentWindowARRAY1 = [];
    ProbabilityOfHittingCurrentWindowARRAY1 = [];
    for CurrentFOXP3WindowSTART = MinFOXP3WindowSTART:StepFOXP3WindowSTART:MaxFOXP3WindowSTART
        CurrentFOXP3WindowEND = CurrentFOXP3WindowSTART + FOXP3WindowWIDTH;
        if MeasurementTime1 >= CurrentFOXP3WindowSTART & MeasurementTime1 <= CurrentFOXP3WindowEND 
            ProbabilityOfMissingCurrentWindow = 0;
        elseif MeasurementTime1 < CurrentFOXP3WindowSTART | MeasurementTime1 > CurrentFOXP3WindowEND
            ProbabilityOfMissingCurrentWindow = 1;
        else
            disp('ERROR');
        end
        ProbabilityOfMissingCurrentWindowARRAY1 = [ProbabilityOfMissingCurrentWindowARRAY1, ProbabilityOfMissingCurrentWindow];
        ProbabilityOfHittingCurrentWindowARRAY1 = [ProbabilityOfHittingCurrentWindowARRAY1, (1-ProbabilityOfMissingCurrentWindow)];
    end
    ProbabilityOfMissingALLWindows = sum(ProbabilityOfMissingCurrentWindowARRAY1)/size(ProbabilityOfMissingCurrentWindowARRAY1,2);
    ProbabilityOfHittingALLWindows = sum(ProbabilityOfHittingCurrentWindowARRAY1)/size(ProbabilityOfHittingCurrentWindowARRAY1,2);
    ProbabilitiesArrayMissing =  [ProbabilitiesArrayMissing, ProbabilityOfMissingALLWindows];
    ProbabilitiesHittingArray =  [ProbabilitiesHittingArray, ProbabilityOfHittingALLWindows];
end

MyFigure = figure;
plot(TimesArray,ProbabilitiesArrayMissing);
ylim([0,1.1]);
xlabel('Time of Measurement (min)');
ylabel('Probability of missing FOXP3 window');

filename = [sprintf('ProbabilityMissingFOXP3_1_Measurements.jpg')];
fullname = fullfile(workingDir,filename);
saveas(MyFigure,fullname);  % here you save the figure

[MinProbOfMissingALL_1_measurement, IndexOfMin] = min(ProbabilitiesArrayMissing);
CorrespondingProbOfHittingALL_1_measurement = ProbabilitiesHittingArray(IndexOfMin);

% %%%%%%%%% 2 MEASUREMENTS %%%%%%%%%
% % Since we can assume that "Measurement 1 hitting FOXP3 window" and 
% % "Measurement 2 hitting FOXP3 window" are independent random variables,
% % we can use the above computed Probability distribution for each of them,
% % and simply combine these P(A) and P(B) using the basic rules for
% % probability
% 
% % P(A and B) = P(A)xP(B)
% % P(A or B) = P(A) + P(B) - P(A and B)
% 
% % A = "missing with M1"
% % B = "missing with M2"
% % P(missing the window) = P("missing with M1" and "missing with M1")
% 
% Z = ProbabilitiesArray'*ProbabilitiesArray;
% 
% % Z = zeros(size(TimesArray,2),size(TimesArray,2))
% % Z(6,3) = 1
% % Z(6,4) = 1
% % Z(7,3) = 1
% % Z(7,4) = 1
% 
% figure;
% surf(TimesArray,TimesArray,Z)
% %plot(TimesArray,ProbabilitiesArray);
% 
% % Not reliable result, thus perhaps do it with for loop to double check.

%%%%%%%%% 2 MEASUREMENT - WITH FOR LOOPS %%%%%%%%%

ProbabilityOfMissingCurrentWindow_ALLtimepoints_MATRIX = [];

for MeasurementTime2 = MinTimePoint:TimePointStep:MaxTimePoint
    ProbabilitiesArrayMissing = []; 
    for MeasurementTime1 = MinTimePoint:TimePointStep:MaxTimePoint
        ProbabilityOfMissingCurrentWindow_ALLtimepoints_ARRAY = [];
        for CurrentFOXP3WindowSTART = MinFOXP3WindowSTART:StepFOXP3WindowSTART:MaxFOXP3WindowSTART
            CurrentFOXP3WindowEND = CurrentFOXP3WindowSTART + FOXP3WindowWIDTH;
            % Measurement 1
            if MeasurementTime1 >= CurrentFOXP3WindowSTART && MeasurementTime1 <= CurrentFOXP3WindowEND 
                ProbabilityOfMissingCurrentWindow1 = 0;
            elseif MeasurementTime1 < CurrentFOXP3WindowSTART || MeasurementTime1 > CurrentFOXP3WindowEND
                ProbabilityOfMissingCurrentWindow1 = 1;
            else
                disp('ERROR');
            end
            % Measurement 2
            if MeasurementTime2 >= CurrentFOXP3WindowSTART && MeasurementTime2 <= CurrentFOXP3WindowEND 
                ProbabilityOfMissingCurrentWindow2 = 0;
            elseif MeasurementTime2 < CurrentFOXP3WindowSTART || MeasurementTime2 > CurrentFOXP3WindowEND
                ProbabilityOfMissingCurrentWindow2 = 1;
            else
                disp('ERROR');
            end
            % Missing with both Measurement 1 and 2
            if ProbabilityOfMissingCurrentWindow1 == 1 && ProbabilityOfMissingCurrentWindow2 == 1 
                ProbabilityOfMissingCurrentWindow_ALLtimepoints = 1;
            elseif ProbabilityOfMissingCurrentWindow1 == 0 || ProbabilityOfMissingCurrentWindow2 == 0
                ProbabilityOfMissingCurrentWindow_ALLtimepoints = 0;
            else
                disp('Error!');
            end
            ProbabilityOfMissingCurrentWindow_ALLtimepoints_ARRAY = ...
                [ProbabilityOfMissingCurrentWindow_ALLtimepoints_ARRAY, ...
                ProbabilityOfMissingCurrentWindow_ALLtimepoints];
        end
        
        ProbabilityOfMissingALLWindows = sum(ProbabilityOfMissingCurrentWindow_ALLtimepoints_ARRAY)/size(ProbabilityOfMissingCurrentWindow_ALLtimepoints_ARRAY,2);
        ProbabilitiesArrayMissing =  [ProbabilitiesArrayMissing, ProbabilityOfMissingALLWindows];
    end
    ProbabilityOfMissingCurrentWindow_ALLtimepoints_MATRIX = ...
    [ProbabilityOfMissingCurrentWindow_ALLtimepoints_MATRIX; ProbabilitiesArrayMissing];
end

Z = ProbabilityOfMissingCurrentWindow_ALLtimepoints_MATRIX; % CHANGE THIS
MinProbOfMissingALL_2_measurement = min(min(Z));

figure;
%subplot(1,2,1);
colormap(flipud(parula));
surf(TimesArray,TimesArray,Z);
caxis([0 1]);
title(['Minimum of Probability = ', num2str(MinProbOfMissingALL_2_measurement)]);
xlabel('Time of Measurement A (min)');
ylabel('Time of Measurement B (min)');
zlabel('Probability of missing FOXP3 with both measurements A and B');
h = colorbar;
ylabel(h, 'Probability of missing FOXP3 with both measurements A and B')

%subplot(1,2,2);
MyFigure = figure;
colormap(flipud(parula));
contourf(TimesArray,TimesArray,Z,'LineStyle','none','LevelStep',0.001);
caxis([0 1]);
title(['Minimum of Probability = ', num2str(MinProbOfMissingALL_2_measurement)]);
%colormap([0,1,0;1,0,0]);
%colormap('winter');
xlabel('Time of Measurement A (min)');
ylabel('Time of Measurement B (min)');
h = colorbar;
ylabel(h, 'Probability of missing FOXP3 with both measurements A and B')

filename = [sprintf('ProbabilityMissingFOXP3_2_Measurements.jpg')];
fullname = fullfile(workingDir,filename);
saveas(MyFigure,fullname);  % here you save the figure


%%%%%%%%% 3 MEASUREMENT - WITH FOR LOOPS %%%%%%%%%
mkdir(workingDir,'images');

ProbabilityOfMissingCurrentWindow_ALLtimepoints_MATRIX_3 = {};

IndexTensor = 0;
for MeasurementTime3 = MinTimePoint:TimePointStep:MaxTimePoint
    ProbabilityOfMissingCurrentWindow_ALLtimepoints_MATRIX = [];
    for MeasurementTime2 = MinTimePoint:TimePointStep:MaxTimePoint
        ProbabilitiesArrayMissing = []; 
        for MeasurementTime1 = MinTimePoint:TimePointStep:MaxTimePoint
            ProbabilityOfMissingCurrentWindow_ALLtimepoints_ARRAY = [];
            for CurrentFOXP3WindowSTART = MinFOXP3WindowSTART:StepFOXP3WindowSTART:MaxFOXP3WindowSTART
                CurrentFOXP3WindowEND = CurrentFOXP3WindowSTART + FOXP3WindowWIDTH;
                % Measurement 1
                if MeasurementTime1 >= CurrentFOXP3WindowSTART && MeasurementTime1 <= CurrentFOXP3WindowEND 
                    ProbabilityOfMissingCurrentWindow1 = 0;
                elseif MeasurementTime1 < CurrentFOXP3WindowSTART || MeasurementTime1 > CurrentFOXP3WindowEND
                    ProbabilityOfMissingCurrentWindow1 = 1;
                else
                    disp('ERROR');
                end
                % Measurement 2
                if MeasurementTime2 >= CurrentFOXP3WindowSTART && MeasurementTime2 <= CurrentFOXP3WindowEND 
                    ProbabilityOfMissingCurrentWindow2 = 0;
                elseif MeasurementTime2 < CurrentFOXP3WindowSTART || MeasurementTime2 > CurrentFOXP3WindowEND
                    ProbabilityOfMissingCurrentWindow2 = 1;
                else
                    disp('ERROR');
                end
                % Measurement 3
                if MeasurementTime3 >= CurrentFOXP3WindowSTART && MeasurementTime3 <= CurrentFOXP3WindowEND 
                    ProbabilityOfMissingCurrentWindow3 = 0;
                elseif MeasurementTime3 < CurrentFOXP3WindowSTART || MeasurementTime3 > CurrentFOXP3WindowEND
                    ProbabilityOfMissingCurrentWindow3 = 1;
                else
                    disp('ERROR');
                end
                % Missing with both Measurement 1 and 2
                if ProbabilityOfMissingCurrentWindow1 == 1 && ProbabilityOfMissingCurrentWindow2 == 1 && ProbabilityOfMissingCurrentWindow3 == 1 
                    ProbabilityOfMissingCurrentWindow_ALLtimepoints = 1;
                elseif ProbabilityOfMissingCurrentWindow1 == 0 || ProbabilityOfMissingCurrentWindow2 == 0 || ProbabilityOfMissingCurrentWindow3 == 0
                    ProbabilityOfMissingCurrentWindow_ALLtimepoints = 0;
                else
                    disp('Error!');
                end
                ProbabilityOfMissingCurrentWindow_ALLtimepoints_ARRAY = ...
                    [ProbabilityOfMissingCurrentWindow_ALLtimepoints_ARRAY, ...
                    ProbabilityOfMissingCurrentWindow_ALLtimepoints];
            end

            ProbabilityOfMissingALLWindows = sum(ProbabilityOfMissingCurrentWindow_ALLtimepoints_ARRAY)/size(ProbabilityOfMissingCurrentWindow_ALLtimepoints_ARRAY,2);
            ProbabilitiesArrayMissing =  [ProbabilitiesArrayMissing, ProbabilityOfMissingALLWindows];
        end
        ProbabilityOfMissingCurrentWindow_ALLtimepoints_MATRIX = ...
        [ProbabilityOfMissingCurrentWindow_ALLtimepoints_MATRIX; ProbabilitiesArrayMissing];
    end
%     ProbabilityOfMissingCurrentWindow_ALLtimepoints_MATRIX
%     disp(' ');
    IndexTensor = IndexTensor + 1;
    ProbabilityOfMissingCurrentWindow_ALLtimepoints_MATRIX_3{IndexTensor} = ...
    ProbabilityOfMissingCurrentWindow_ALLtimepoints_MATRIX;

    Z = ProbabilityOfMissingCurrentWindow_ALLtimepoints_MATRIX; % CHANGE THIS
    MinProbOfMissingALL_CurrentWindow = min(min(Z));
    
    MyFigure = figure;
    set(MyFigure, 'Visible', 'off');
    colormap(flipud(parula));
    contourf(TimesArray,TimesArray,Z,'LineStyle','none','LevelStep',0.001);
    caxis([0 1]);
    title(['Time of Measurement C = ', num2str(IndexTensor*TimePointStep-TimePointStep) ,' min, '...
        'Minimum of Probability = ', num2str(MinProbOfMissingALL_CurrentWindow)]);
    xlabel('Time of Measurement A (min)');
    ylabel('Time of Measurement B (min)');
    h = colorbar;
    ylabel(h, 'Probability of missing FOXP3 with both measurements A, B and C');
    
    %saveas(MyFigure,['Image', num2str(IndexTensor),'.jpg']);  % here you save the figure
    filename = [sprintf('%03d',IndexTensor) '.jpg'];
    fullname = fullfile(workingDir,'images',filename);
    saveas(MyFigure,fullname);  % here you save the figure
end

%%% Make video %%%
shuttleVideo = VideoReader('shuttle.avi');
imageNames = dir(fullfile(workingDir,'images','*.jpg'));
imageNames = {imageNames.name}';
outputVideo = VideoWriter(fullfile(workingDir,'ProbabilityMissingFOXP3_3_Measurements.avi'));
SlowDown = 1.0/(2*TimePointStep);
outputVideo.FrameRate = SlowDown*shuttleVideo.FrameRate;
open(outputVideo)
for ii = 1:length(imageNames)
   img = imread(fullfile(workingDir,'images',imageNames{ii}));
   writeVideo(outputVideo,img)
end
close(outputVideo)
shuttleAviMovieAvi = VideoReader(fullfile(workingDir,'ProbabilityMissingFOXP3_3_Measurements.avi'));
ii = 1;
while hasFrame(shuttleAviMovieAvi)
   mov(ii) = im2frame(readFrame(shuttleAviMovieAvi));
   ii = ii+1;
end
%%% End video %%%

MinARRAY = [];
for i = 1:IndexTensor
    Matrix = ProbabilityOfMissingCurrentWindow_ALLtimepoints_MATRIX_3{i};
    Min2D = min(min(Matrix));
    MinARRAY = [MinARRAY, Min2D];
end
MinProbOfMissingALL_3_measurement = min(MinARRAY);


%%%%%%%%%% Probability Function of Timepoints %%%%%%%%%%

ProbabilityMISSING_ALLtimepoints_1patient = [1, MinProbOfMissingALL_1_measurement, ...
    MinProbOfMissingALL_2_measurement, MinProbOfMissingALL_3_measurement];
ProbabilityHITTING_AtLeast1Timepoint_1patient = [1,1,1,1] - ProbabilityMISSING_ALLtimepoints_1patient;
% CorrespondingProbOfHittingALL_1_measurement, ...

ProbabilityMISSING_ALLtimepoints_ALLpatients = ProbabilityMISSING_ALLtimepoints_1patient.^3;
ProbabilityHITTING_AtLeast1Timepoint_ALLpatients = ProbabilityHITTING_AtLeast1Timepoint_1patient.^3;

MyFigure = figure;
MeasurementsNumber = [0,1,2,3];

subplot(2,1,1);
plot(MeasurementsNumber,ProbabilityMISSING_ALLtimepoints_1patient,'-o', ...
    'Color','r','MarkerSize',8,'MarkerEdgeColor','r','MarkerFaceColor','r'); hold on;
plot(MeasurementsNumber,ProbabilityHITTING_AtLeast1Timepoint_1patient,'-o', ...
    'Color','g','MarkerSize',8,'MarkerEdgeColor','g','MarkerFaceColor','g'); hold off;
title('1 Patient')
xticks(MeasurementsNumber);
xlabel('Number of intermediate measurements (between 0h and 6h)');
ylabel('Probability');
legend('MINIMUM of Probability of MISSING with ALL measurements', ...
    'MAXIMUM of Probability of HITTING with AT LEAST 1 measurement');

subplot(2,1,2);
plot(MeasurementsNumber,ProbabilityMISSING_ALLtimepoints_ALLpatients,'-o', ...
    'Color','r','MarkerSize',8,'MarkerEdgeColor','r','MarkerFaceColor','r'); hold on;
plot(MeasurementsNumber,ProbabilityHITTING_AtLeast1Timepoint_ALLpatients,'-o', ...
    'Color','g','MarkerSize',8,'MarkerEdgeColor','g','MarkerFaceColor','g'); hold on;
title('3 Patients')
xticks(MeasurementsNumber);
xlabel('Number of intermediate measurements (between 0h and 6h)');
ylabel('Probability');
legend('MINIMUM of Probability of MISSING with ALL measurements for ALL 3 patients', ...
    'MAXIMUM of Probability of HITTING with AT LEAST 1 measurement for ALL 3 patients');

suptitle('Probabilities of MISSING/HITTING FOXP3 window');

filename = [sprintf('MinimumProbMissingALLFOXP3.jpg')];
fullname = fullfile(workingDir,filename);
saveas(MyFigure,fullname);  % here you save the figure

%%% Play video %%%
figure
imshow(mov(1).cdata, 'Border', 'tight')
movie(mov,100,shuttleAviMovieAvi.FrameRate)



