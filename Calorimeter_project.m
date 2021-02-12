%%Header 
%Linus Schmitz 
%Final Version 10/25/2020

%% Givens 
clear 
clc
close 
%Equation 1
%Import givens for constant values in the equation
% C_sav; % specific heat of sample
C_cav = 0.895; %specific heat of the calorimeter in J/gC
% T_0; %initial temperature of the calorimeter
% T_1; %initial temperature of the sample
% T_2; %final temperature of the calorimeter and sample at equilibrium
m_c = 510; %mass of the calorimeter in g
m_s = 39.306; %mass of the sample in g 

%% Reallocating the data into indivdiual vectors
%Seperate the data into vectors for easier handeling
data = importdata('SampleB');
field = 'textdata';
data = rmfield(data, field);
data = struct2array(data);
data(:,1)=[];
time = data(:,1);
CalTemp1 = data(:,2);
SampleTemp = data(:,3);
AirTemp = data(:,4);
CalTemp2 = data(:,5);

%Combing the two calorimeter temperature vectors
CalTemp = (CalTemp1 + CalTemp2)/2;

%% Plotting initial data for understanding 
%Plot all of the data on a single plot to make sense of it. Turn various
%parts on and off to get a better view
figure(1)
hold on
title('Time vs Temperature without linear fit')
time1 = seconds(time);
time1.Format = 'mm:ss';
plot(time1, CalTemp)
plot(time1, SampleTemp)
plot(time1, AirTemp)
legend('Calorimeter Temperature','Sample Temperature', 'Air Temperature', 'Location','East' )
ylabel('Temperature in C')
xlabel('Time in Minutes')
hold off

%% Important times 
%Looking at the total plot and the raw data I made best estimates for when
%significant times occured during the test process
time_start = 1;     % This is the first value
time_sample = 310;  % This is the time where the sample was placed into the calorimeter
time_stop = 410;    % This is the time where the sample stopped warming the calorimeter

 %% Best Fit line first linear (1) tempurature 
 %Using the function solve_line_fit plot linear regressions for each
 %significant time segments of the temperature of the calorimeter. Note:
 %not all of these will be used however durin the coding process I wanted
 %to make sure that the function was working and I wanted to see if any of
 %the other linear regressions might be useful so plotted the extra
 %segments
 
 %In addition total error was calculated for each of the line segments.
 %This is the 4th line of each code section 
 
 %Initial time segment. From start to right before the sample was added
t1 = time(time_start:time_sample,:);
y1 = CalTemp(time_start:time_sample,:);
[A1,B1, A_1_error, B_1_error] = solve_line_fit(2, t1,y1,0);
T1_error = sqrt((A_1_error)^2+(B_1_error)^2);     %Total Error calculation
G1 = @(t0) A1 + B1 * t0;

%Second time segment. This is an extra line, from when the sample was added
%to when the calorimeter and sample were at equilibrium 
t2 = time(time_sample:(time_stop),:);
y2 = CalTemp(time_sample:(time_stop),:);
[A2,B2, A_2_error, B_2_error] = solve_line_fit(3,t2,y2,-1);
T2_error = sqrt((A_2_error)^2+(B_2_error)^2);   %Total Error calculation
G2 = @(t2) A2 + B2 * t2;

%Third time segment. From the start of when the sample and calorimeter are
%at equilibrium to the end of the experiment 
t3 = time(time_stop:end,:);
y3 = CalTemp(time_stop:end,:);
[A3,B3, A_3_error, B_3_error] = solve_line_fit(4,t3,y3,0);
T3_error = sqrt((A_3_error)^2+(B_3_error)^2);   %Total Error calculation
G3 = @(t0) A3 + B3 * t0;

%This is an extra segment and is the sample temperature line of best fit
%for the entire experiment
t4 = time(1:end,:);
y4 = SampleTemp(1:end,:);
[A4,B4, A_4_error, B_4_error] = solve_line_fit(5,t4,y4,0);
T4_error = sqrt((A_4_error)^2+(B_4_error)^2);   %Total Error calculation
G4 = @(t0) A4 + B4 * t0;

%% Start solving for the unknowns in the equation
%This is the time that all the variable values will be calculated at. 
time_testing = ceil((413+312)/2);

%Finding the Temperature values based off the the functions that were
%written with the solve_line_fit outputs for each line segment. All are in
%a y = mx + b format 
T_0 = G1(time_testing);
T_1 = SampleTemp(time_testing);
T_2 = G3(time_testing);

%Solving for C_sav
C_sav = (m_c*C_cav*(T_2-T_0))/(m_s*(T_1-T_2));

%Declaring error in the two masses based on significant figures in the mass
T5_error = .5; 
T6_error = 0.0005; 

%Partial derivatives for each of the variables and their coresponding error
%calculations 
error_2 = (((m_c * C_cav * (T_1 - T_0))/(m_s * (T_1-T_2)^2)) * T3_error);
error_1 = (-((m_c * C_cav * (T_2 - T_0))/(m_s * (T_1-T_2)^2)) * T4_error);
error_0 = (-((m_c * C_cav)/(m_s * (T_1 -T_2)) * T1_error));
error_3 = ((( C_cav * (T_2 - T_0))/(m_s * (T_1-T_2))) *T5_error );
error_4 = (-((m_c * C_cav * (T_2 - T_0))/(m_s^2 * (T_1-T_2))) * T6_error );

format long
%Total error calculation
Error_C_sav = sqrt((error_2)^2 + (error_1)^2 +(error_0)^2 +(error_3)^2 +(error_4)^2);
fprintf("The C_sav value calculated was %.3f (J/gC) and the associated error in that value was +- %.3f (J/gC)\n", C_sav, Error_C_sav)


% solve for the missing temp values and then plug them into the equation
% from beggining to get C_cav or whatever

%% Put everything on the same plot
figure(7)
hold on
title('Time vs Temperature with linear fits')
plot(time1, CalTemp)
plot(time1(time_start:time_sample+80,:), G1(time_start:time_sample+80), 'LineWidth', 2);
%plot(time(time_sample:time_stop,:), G2(time(time_sample:time_stop,:)), 'LineWidth', 2);
plot(time1(time_stop-80:end,:), G3(time(time_stop-80:end,:)), 'LineWidth', 2);
xline(time1(time_testing), 'c')
plot(time1(time_testing),G1(time_testing), 'b*');

%The subtraction of 0.1 is for visualization purposes. I can not get the
%plot to display properly despite the functions working correctly 
plot(time1(time_testing),G3(time_testing)-0.1, 'm*');
avgTime = time1(time_testing);


lgd = legend('Calorimeter Temp', 'Calorimeter w/o sample', 'Calorimeter & sample @ equilibrium', 'Average Time (used for calcs)', ...
'Temp of Calorimeter (T_0)','Temp of Cal and Sample (T_2)', 'Location', 'East');
lgd.FontSize = 7;
ylabel('Temperature in C')
xlabel('Time in Minutes')


hold off







