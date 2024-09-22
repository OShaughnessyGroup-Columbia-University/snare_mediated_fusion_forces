clear all

filename_fl = './force_data_fixed_coord_nostaple_nrod_7_rtmd_2.27.txt';
%filename_cut = './force_data_fixed_coord_cutsnare_nostaple_nrod_7_rtmd_2.27.txt';
filename_cut = './force_data_fixed_coord_cutsnare_nrod_7_rtmd_2.27.txt';
nrod = 7;

fileID = fopen(filename_fl, 'r');
% Read the file line by line
while ~feof(fileID)
    % Read the current line
    currentLine = fgetl(fileID);
    
    if contains (currentLine, 'Number of measurements per SNARE')
        nextLine = fgetl(fileID);
        n_measure_fl = str2double(nextLine);
    end
    % Check if the current line contains the search phrase
    if contains(currentLine, 'Radial entropic force (pN)')
        nextLine = fgetl(fileID);
        mean_fl = str2double(nextLine);
        
        nextLine = fgetl(fileID);
        sd_fl = str2double(nextLine);
        break;
    end
end


fileID = fopen(filename_cut, 'r');
% Read the file line by line
while ~feof(fileID)
    % Read the current line
    currentLine = fgetl(fileID);
    
    if contains (currentLine, 'Number of measurements per SNARE')
        nextLine = fgetl(fileID);
        n_measure_cut = str2double(nextLine);
    end
    % Check if the current line contains the search phrase
    if contains(currentLine, 'Radial entropic force (pN)')
        nextLine = fgetl(fileID);
        mean_cut = str2double(nextLine);
        
        nextLine = fgetl(fileID);
        sd_cut = str2double(nextLine);
        break;
    end
end

x_arr = [1,2];
x_arr = categorical(x_arr);

m_f_sq = [mean_fl, mean_cut];
npt = [n_measure_fl*nrod, n_measure_cut*nrod];
err_f_sd = [sd_fl, sd_cut]; % SD
err_f_sem = err_f_sd./sqrt(npt); % SEM

col_arr = [0.94,0,0.82;
    0.05,0.62,0.31;
    0,0,0.85;
    0.92,0.39,0.07;
    0.93,0.69,0.13];

col_fl = col_arr(3,:);
col_cut = col_arr(5,:);

figure
hold on
%bar(x_arr,m_f_sq, 0.6,'FaceColor', col_fl)
%errorbar(x_arr,m_f_sq,err_f_sd,'black','LineWidth',5,'Capsize',12)
bar(x_arr(1),m_f_sq(1), 0.6,'FaceColor', col_fl)
bar(x_arr(2),m_f_sq(2), 0.6,'FaceColor', col_cut)
errorbar(x_arr(1),m_f_sq(1),err_f_sd(1),'black','LineWidth',5,'Capsize',12)
errorbar(x_arr(2),m_f_sq(2),err_f_sd(2),'black','LineWidth',5,'Capsize',12)
set(gca,'xticklabel',{'Full-length', 'Truncated'});

ax = gca;
ax.FontSize = 22;
ax.FontWeight="bold";
ax.XMinorTick="off";
ax.YMinorTick="on";
ax.Box = 'on';
ax.TickLength = [0.03 0.025];
ax.LineWidth = 1;
%ax.XLim([0.5,3.5])
ylabel('Radial entropic force per SNARE (pN)','Fontsize',22)
