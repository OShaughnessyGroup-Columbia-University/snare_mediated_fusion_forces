clear all

folder_prefix = './';
nrod_arr = [3,4,5,6,7,9];
nsn = length(nrod_arr);
col_arr = [0.94,0,0.82; % magenta
    0.93,0.69,0.13; % mustard yellow
    0.05,0.62,0.31; % green
    0,0,0.85; % blue
    0.92,0.39,0.07; % orange
    0.05,0.18,0.61  % navy blue
    ];
ncorr = 1;
SNAREs = zeros(nsn, 1);
r_ini = zeros(nsn, 1);
n_measure = zeros(nsn, 1);
n_measure_tot = zeros(nsn, 1);
TMD_radius = zeros(nsn, 2);
Zippering_force = zeros(nsn, 2);
Total_squeezing_force = zeros(nsn, 2);
Squeezing_force_per_LD = zeros(nsn, 2);
Radial_linker_force = zeros(nsn, 2);
Theta_linker_force = zeros(nsn, 2);
Radial_entropic_force = zeros(nsn, 2);
Theta_entropic_force = zeros(nsn, 2);
Z_entropic_force = zeros(nsn, 2);


for i = 1:length(nrod_arr)
    
    nrod = nrod_arr(i);
    
    fname = [folder_prefix 'force_data_nrod_' int2str(nrod) '_rtmd_ini_5.91.txt'];
    fileID = fopen(fname, 'r');

    % Initialize an empty array to store numbers
    numberArray = [];
        
    % Read the file line by line
    while ~feof(fileID)
        line = fgetl(fileID);
            
        % Convert the line to a number, returns NaN if it's not a number
        num = str2double(line);
            
        % Check if the line is a number (not NaN)
        if ~isnan(num)
            numberArray = [numberArray; num];
        end
    end
        
    % Close the file
    fclose(fileID);
        
    % Extract the data
    SNAREs(i) = numberArray(1);
    r_ini(i) = numberArray(2);
    n_measure(i) = numberArray(3);
    n_measure_tot(i) = numberArray(3)*nrod;
    if n_measure(i) ~= 0
        TMD_radius(i, :) = numberArray(4:5);
        Zippering_force(i, :) = numberArray(6:7);
        Total_squeezing_force(i, :) = numberArray(8:9);
        Squeezing_force_per_LD(i, :) = numberArray(10:11);
        Radial_linker_force(i, 1) = -numberArray(12);
        Radial_linker_force(i, 2) = numberArray(13);
        Theta_linker_force(i, :) = numberArray(14:15);
        Radial_entropic_force(i, :) = numberArray(16:17);
        Theta_entropic_force(i, :) = numberArray(18:19);
        Z_entropic_force(i, :) = numberArray(20:21);
    else
        TMD_radius(i, :) = -1;
        Zippering_force(i, :) = -1;
        Total_squeezing_force(i, :) = -1;
        Squeezing_force_per_LD(i, :) = -1;
        Radial_linker_force(i, :) = -1;
        Theta_linker_force(i, :) = -1;
        Radial_entropic_force(i, :) = -1;
        Theta_entropic_force(i, :) = -1;
        Z_entropic_force(i, :) = -1;
    end
end

x_arr = categorical(nrod_arr);

% Equil. ring radius 
col = col_arr(2,:);
figure;
hold on
bar(x_arr,TMD_radius(:,1), 0.8,'FaceColor', col)
errorbar(x_arr,TMD_radius(:,1),TMD_radius(:,2),'black','LineWidth',5,'Capsize',12,'LineStyle', 'none')
%errorbar(x_arr,TMD_radius(:,1),TMD_radius(:,2)./n_measure,'black','LineWidth',5,'Capsize',12,'LineStyle', 'none')

xlabel('# SNAREs');
ylabel('Ring radius (nm)');
ax = gca;
ax.FontSize = 16;
ax.FontWeight="bold";
ax.XMinorTick="off";
ax.YMinorTick="on";
ax.Box = 'on';
ax.TickLength = [0.03 0.025];
ax.LineWidth = 1;
%ax.XLim = ([1,6]);
%ax.YLim = ([0,14]);

% f_squeeze per SNARE 
col = col_arr(1,:);
figure;
hold on
bar(x_arr,Squeezing_force_per_LD(:,1), 0.8,'FaceColor', col)
errorbar(x_arr,Squeezing_force_per_LD(:,1),Squeezing_force_per_LD(:,2),'black','LineWidth',5,'Capsize',12,'LineStyle', 'none')
%errorbar(x_arr,Squeezing_force_per_LD(:,1),Squeezing_force_per_LD(:,2)./n_measure_tot,'black','LineWidth',5,'Capsize',12,'LineStyle', 'none')

xlabel('# SNAREs');
ylabel('Squeezing force per SNARE (pN)');
ax = gca;
ax.FontSize = 16;
ax.FontWeight="bold";
ax.XMinorTick="off";
ax.YMinorTick="on";
ax.Box = 'on';
ax.TickLength = [0.03 0.025];
ax.LineWidth = 1;
%ax.XLim = ([1,6]);
ax.YLim = ([0,30]);

% F_squeeze total
col = col_arr(6,:);
figure;
hold on
bar(x_arr,Total_squeezing_force(:,1), 0.8,'FaceColor', col)
errorbar(x_arr,Total_squeezing_force(:,1),Total_squeezing_force(:,2),'black','LineWidth',5,'Capsize',12,'LineStyle', 'none')
%errorbar(x_arr,Total_squeezing_force(:,1),Total_squeezing_force(:,2)./n_measure,'black','LineWidth',5,'Capsize',12,'LineStyle', 'none')

xlabel('# SNAREs');
ylabel('Total squeezing force (pN)');
ax = gca;
ax.FontSize = 16;
ax.FontWeight="bold";
ax.XMinorTick="off";
ax.YMinorTick="on";
ax.Box = 'on';
ax.TickLength = [0.03 0.025];
ax.LineWidth = 1;
%ax.XLim = ([1,6]);
%ax.YLim = ([0,14]);


% Radial forces
figure;
hold on;
% Combine the data for the bar plot
data = [Radial_linker_force(:,1), Radial_entropic_force(:,1)];
% Define the error data
%error_data = [Radial_linker_force(:,2), Radial_entropic_force(:,2)];
error_data = [Radial_linker_force(:,2)./sqrt(n_measure_tot/ncorr), Radial_entropic_force(:,2)./sqrt(n_measure_tot/ncorr)];
% Create a grouped bar plot
bar_handle = bar(x_arr, data, 'grouped','BarWidth', 1);
bar_handle(1).FaceColor = col_arr(3,:);
bar_handle(2).FaceColor = col_arr(5,:);

% Set the x-axis to have categorical labels
set(gca, 'XTickLabel', x_arr);

% Plot error bars
for i = 1:size(data, 2)
    x = bar_handle(i).XEndPoints; % Get the x positions for error bars
    errorbar(x, data(:, i), error_data(:, i),'black','LineWidth',5,'Capsize',0,'LineStyle', 'none');
end

ax = gca;
ax.FontSize = 16;
ax.FontWeight="bold";
ax.XMinorTick="off";
ax.YMinorTick="on";
ax.Box = 'on';
ax.TickLength = [0.03 0.025];
ax.LineWidth = 1;
%xlabel('# SNAREs');
%ylabel('Radial force per SNARE (pN)');
set(gca, 'XTickLabel', []);
xlabel('');
set(gca, 'YTickLabel', []);
ylabel('');
hold off