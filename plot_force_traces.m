clear all

col_arr = [0.94,0,0.82;
    0.05,0.62,0.31;
    0,0,0.85;
    0.92,0.39,0.07;
    0.93,0.69,0.13];

folder_prefix = './';
%fig 2
%fname = 'ves_ves_avg_nostaple_rtmd_2.27_deltah_2.10_Nunzip_10.00_lp_0.50_nrod_7_rod_r_0.38_dves_55.9_tension_0.05_fent_rad.mat'; 
fname = 'ves_ves_cutsnare_rtmd_2.27_deltah_2.10_Nunzip_10.00_lp_0.50_nrod_7_rod_r_0.38_dves_55.9_tension_0.05_fent_rad.mat';

%fig 3
%fname = 'ves_ves_avg_nostaple_rtmd_ini_2.27_deltah_2.10_Nunzip_10.00_lp_0.50_nrod_7_rod_r_0.38_dves_55.9_tension_0.05_fent_rad.mat';
%fname = 'ves_ves_cutsnare_nostaple_rtmd_ini_2.27_deltah_2.10_Nunzip_10.00_lp_0.50_nrod_7_rod_r_0.38_dves_55.9_tension_0.05_fent_rad.mat';
%fname = 'ves_ves_cutsnare_rtmd_ini_2.27_deltah_2.10_Nunzip_10.00_lp_0.50_nrod_7_rod_r_0.38_dves_55.9_tension_0.05_fent_rad.mat';

fname = [folder_prefix fname];
data = load(fname);
t = data.array(:,1);
fent = data.array(:,2);

col = col_arr(5,:);
%{
figure;
plot(t, fent,LineWidth=1,Color=col);
xlabel('Time (\mus)');
ylabel('Radial entropic force (pN)');
ax = gca;
ax.FontSize = 16;
ax.FontWeight="bold";
ax.XMinorTick="off";
ax.YMinorTick="on";
ax.Box = 'on';
ax.TickLength = [0.03 0.025];
ax.LineWidth = 1;
ax.XLim = ([0,5]);
ax.YLim = ([-20,60]);
%}


%figure;

binEdges = -100:2:100;
fent = fent(floor(length(fent) / 2)+1:end);
histogram(fent, binEdges, 'Normalization', 'pdf', 'FaceColor', col);

xlabel('Radial entropic force (pN)');
ylabel('Probability Density (pN^{-1))');
ax = gca;
ax.FontSize = 16;
ax.FontWeight="bold";
ax.XMinorTick="off";
ax.YMinorTick="on";
ax.Box = 'on';
ax.TickLength = [0.03 0.025];
ax.LineWidth = 1;
ax.XLim = ([-40,60]);
ax.YLim = ([0,0.08]);
%}
