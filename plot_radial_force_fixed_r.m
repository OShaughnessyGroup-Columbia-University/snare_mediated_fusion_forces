clear all

folder_prefix = './rtmd_fix_datafiles/';
nrod_arr = [3,5,7,9];
col_arr = [0.94,0,0.82;
    0.05,0.62,0.31;
    0,0,0.85;
    0.92,0.39,0.07;
    0.93,0.69,0.13];
figure;
hold on

for irod = 1:length(nrod_arr)
    
    nrod = nrod_arr(irod);
    
    fname1 = [folder_prefix 'force_data_fixed_coord_nrod_' int2str(nrod) '_rtmd_*.txt'];
    fname2 = [folder_prefix 'force_data_fixed_coord_nostaple_nrod_' int2str(nrod) '_rtmd_*.txt'];
    fname_pat1 = ['force_data_fixed_coord_nrod_' int2str(nrod) '_rtmd_%f.txt'];
    fname_pat2 = ['force_data_fixed_coord_nostaple_nrod_' int2str(nrod) '_rtmd_%f.txt'];

    col = col_arr(irod,:);
    
    % Get the list of files
    files1 = dir(fname1);
    files2 = dir(fname2);
    
    % Extract the numeric part from the filenames and sort them
    fileNames1 = {files1.name};
    fileNames2 = {files2.name};
    fileNames = [fileNames1, fileNames2];
    lastNumber = cellfun(@(x) sscanf(x, fname_pat1), fileNames1);
    lastNumber = [lastNumber, cellfun(@(x) sscanf(x, fname_pat2), fileNames2)];
    [~, sortedIndices] = sort(lastNumber);
    sortedFileNames = fileNames(sortedIndices);
    
    % Preallocate arrays to store the data
    numFiles = length(sortedFileNames);
    nmax = numFiles;%-3;
    SNAREs = zeros(numFiles, 1);
    r_ini = zeros(numFiles, 1);
    n_measure = zeros(numFiles, 1);
    TMD_radius = zeros(numFiles, 2);
    Zippering_force = zeros(numFiles, 2);
    Total_squeezing_force = zeros(numFiles, 2);
    Squeezing_force_per_LD = zeros(numFiles, 2);
    Radial_linker_force = zeros(numFiles, 2);
    Theta_linker_force = zeros(numFiles, 2);
    Radial_entropic_force = zeros(numFiles, 2);
    Theta_entropic_force = zeros(numFiles, 2);
    Z_entropic_force = zeros(numFiles, 2);
    
    % Loop through the files and read the data
    for i = 1:numFiles
        filename = [folder_prefix sortedFileNames{i}];
        fileID = fopen(filename, 'r');
        
        % Check if the file opened successfully
        if fileID == -1
            error('File cannot be opened');
        end
        
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
        
        fclose(fileID);
        
        % Extract the data
        SNAREs(i) = numberArray(1);
        r_ini(i) = numberArray(2);
        n_measure(i) = numberArray(3);
        if n_measure(i) ~= 0
            TMD_radius(i, :) = numberArray(4:5);
            Zippering_force(i, :) = numberArray(6:7);
            Total_squeezing_force(i, :) = numberArray(8:9);
            Squeezing_force_per_LD(i, :) = numberArray(10:11);
            Radial_linker_force(i, :) = -numberArray(12:13);
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
   
    errorbar(r_ini(1:nmax), Radial_entropic_force(1:nmax, 1), Radial_entropic_force(1:nmax, 2)./sqrt(n_measure(1:nmax)), 'o', ...
    'MarkerSize', 6, 'MarkerFaceColor', col, 'MarkerEdgeColor', col, 'LineWidth', 1, 'Color', col);
    errorbar(r_ini(1:nmax), Radial_linker_force(1:nmax, 1), Radial_linker_force(1:nmax, 2)./sqrt(n_measure(1:nmax)), 'o', ...
    'MarkerSize', 6, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', col, 'LineWidth', 1, 'Color', col);
end

xlabel('Ring radius (nm)');
ylabel('Radial force (pN)');
ax = gca;
ax.FontSize = 16;
ax.FontWeight="bold";
ax.XMinorTick="off";
ax.YMinorTick="on";
ax.Box = 'on';
ax.TickLength = [0.03 0.025];
ax.LineWidth = 1;
ax.XLim = ([0,6.5]);
ax.YLim = ([6,15]);
lgd = legend('3 SNAREs', '', '5 SNAREs', '', '7 SNAREs', '', '9 SNAREs', '');
set(lgd, 'Box', 'off');
set(lgd, 'Color', 'none');
