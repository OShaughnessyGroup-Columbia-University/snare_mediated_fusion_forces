clear all

col_arr = [0.94,0,0.82;
    0.05,0.62,0.31;
    0,0,0.85;
    0.92,0.39,0.07;
    0.93,0.69,0.13];

folder_prefix = './';
%fname = 'ves_ves_avg_nostaple_rtmd_ini_2.27_deltah_2.10_Nunzip_10.00_lp_0.50_nrod_7_rod_r_0.38_dves_55.9_tension_0.05_ring_rad_movavg.mat';
%fname = 'ves_ves_cutsnare_nostaple_rtmd_ini_2.27_deltah_2.10_Nunzip_10.00_lp_0.50_nrod_7_rod_r_0.38_dves_55.9_tension_0.05_ring_rad_movavg.mat';
fname = 'ves_ves_cutsnare_rtmd_ini_2.27_deltah_2.10_Nunzip_10.00_lp_0.50_nrod_7_rod_r_0.38_dves_55.9_tension_0.05_ring_rad_movavg.mat';

fname = [folder_prefix fname];
data = load(fname);
t = data.array(:,1)/1000;
fent = data.array(:,2);
fent_sd = data.array(:,3);

col = col_arr(5,:);
%figure;
hold on
plot(t, fent,LineWidth=1,Color=col);
fill([t; flipud(t)], [fent - fent_sd; flipud(fent + fent_sd)], ...
    col, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
xlabel('Time (\mus)');
ylabel('Ring radius (nm)');
ax = gca;
ax.FontSize = 16;
ax.FontWeight="bold";
ax.XMinorTick="off";
ax.YMinorTick="on";
ax.Box = 'on';
ax.TickLength = [0.03 0.025];
ax.LineWidth = 1;
ax.XLim = ([0,2]);
ax.YLim = ([0,7]);
