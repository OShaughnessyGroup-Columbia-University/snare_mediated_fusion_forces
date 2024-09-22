clear all

folder_prefix = '../paper_code/';
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

tfus_arr = [4.45e-02,9.16e-02,6.35e-02;
    1.28e-01,1.12e+00,3.98e+00;
    5.78e+00,-1,-1;]; % from fusion_times_process_ld_scan.ipynb
nfus_arr = [5,5,5;
    5,4,2;
    2,1000,1000;];

figure;
hold on

for j = 1:length(lp_arr)
    lp = lp_arr(j);

    x_arr = categorical(Nunz_arr);
    
    % f_squeeze per SNARE 
    col = col_arr(j,:);
    tfus = tfus_arr(:,j);
    nfus = nfus_arr(:,j);
    
    errorbar(x_arr,tfus,tfus./sqrt(nfus),'black','LineWidth',3,'Capsize',10,'LineStyle', 'none', 'Marker', 'o', 'Color', col)
end

xlabel('# Unzippered LD residues');
ylabel('Fusion time (ms)');
ax = gca;
ax.FontSize = 16;
ax.FontWeight="bold";
ax.XMinorTick="off";
ax.YMinorTick="on";
ax.Box = 'on';
ax.TickLength = [0.03 0.025];
ax.LineWidth = 1;
%ax.XLim = ([1,6]);
ax.YLim = ([0,5]);
legend({'l_{p}=0.3', 'l_{p}=0.5', 'l_{p}=0.7'});
% Get the legend handle
lgd = legend;
set(lgd, 'Color', 'none');
set(lgd, 'EdgeColor', 'none');