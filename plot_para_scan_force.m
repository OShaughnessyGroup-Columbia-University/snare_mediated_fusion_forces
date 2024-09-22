clear all

folder_prefix = './';
nrod = 7;
lp_arr = [0.3,0.5,0.7];
Nunz_arr = [7,10,14];
nsn = length(Nunz_arr);
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

figure;
hold on

for j = 1:length(lp_arr)
    lp = lp_arr(j);

    for i = 1:length(Nunz_arr)
        
        nunz = Nunz_arr(i);
        
        fname = [folder_prefix 'force_data_scan_nrod_' int2str(nrod) '_Nunzip_' num2str(nunz, '%.2f') '_lp_' num2str(lp, '%.2f') '.txt'];
        fileID = fopen(fname, 'r');
    
        % Initialize an empty array to store numbers
        numberArray = [];
            
        % Read the file line by line
        while ~feof(fileID)
            line = fgetl(fileID); % Read a line from the file
                
            % Convert the line to a number, returns NaN if it's not a number
            num = str2double(line);
                
            % Check if the line is a number (not NaN)
            if ~isnan(num)
                numberArray = [numberArray; num]; % Append the number to the array
            end
        end
            
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
    
    x_arr = categorical(Nunz_arr);
    
    % f_squeeze per SNARE 
    col = col_arr(j,:);
    
    %bar(x_arr,Zippering_force(:,1), 0.8,'FaceColor', col)
    errorbar(x_arr,Zippering_force(:,1),Zippering_force(:,2),'black','LineWidth',3,'Capsize',10,'LineStyle', 'none', 'Marker', 'o', 'Color', col)
end

xlabel('# Unzippered LD residues');
ylabel('Zippering force per LD (pN)');
ax = gca;
ax.FontSize = 16;
ax.FontWeight="bold";
ax.XMinorTick="off";
ax.YMinorTick="on";
ax.Box = 'on';
ax.TickLength = [0.03 0.025];
ax.LineWidth = 1;
%ax.XLim = ([1,6]);
ax.YLim = ([0,50]);
legend({'l_{p}=0.3', 'l_{p}=0.5', 'l_{p}=0.7'});
% Get the legend handle
lgd = legend;
set(lgd, 'Color', 'none');
set(lgd, 'EdgeColor', 'none');