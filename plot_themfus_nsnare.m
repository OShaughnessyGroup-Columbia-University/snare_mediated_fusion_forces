clear
close all

N_snare = [3,4,5,6,7,9];
source_cdf = {'statistics_survivals_n_snare_3_lp_0.5_27aug24.txt',...
    'statistics_survivals_n_snare_4_lp_0.5_27aug24.txt',...
    'statistics_survivals_n_snare_5_lp_0.5_27aug24.txt',...
    'statistics_survivals_n_snare_6_lp_0.5_27aug24.txt',...
    'statistics_survivals_n_snare_7_lp_0.5_27aug24.txt',...
    'statistics_survivals_n_snare_9_lp_0.5_27aug24.txt'};
source_cha_time = {'statistics_chara_time_n_snare_3_lp_0.5.txt',...
    'statistics_chara_time_n_snare_4_lp_0.5.txt',...
    'statistics_chara_time_n_snare_5_lp_0.5.txt',...
    'statistics_chara_time_n_snare_6_lp_0.5.txt',...
    'statistics_chara_time_n_snare_7_lp_0.5.txt',...
    'statistics_chara_time_n_snare_9_lp_0.5.txt'};

lwidth=5;
color_arr=['c','r','y','b','k','g'];
capsize=12;
bar_width=0.8;
N_fuse = [0,0,2,9,17,19];
N_hemi = [6,17,20,19,20,20];
fig1=figure(1);
hold on;
fig2=figure(2);
hold on;
fig3=figure(3);
hold on;
fig4=figure(4);
hold on
fig5=figure(5);
hold on
fig6=figure(6);
hold on
fig7=figure(7);
hold on
fig8=figure(8);
hold on

for i = 1:1:numel(N_snare)
    A = readtable(char(source_cdf(i)));
    Data = A.Variables;
    mean_hf = Data(1,1:end);
    low_hf = Data(5,1:end);
    high_hf = Data(9,1:end);
    mean_f = Data(2,1:end);
    low_f = Data(6,1:end);
    high_f = Data(10,1:end);
    mean_h2f = Data(3,1:end);
    low_h2f = Data(7,1:end);
    high_h2f = Data(11,1:end);
    mean_ih = Data(4,1:end);
    low_ih = Data(8,1:end);
    high_ih = Data(12,1:end);
    
    label = num2str(N_snare(i));
    t = [0:1:1499] * 1.36 / 1000;
    color = color_arr(i);
    len_i_char=categorical(cellstr(num2str(N_snare(i))));
    figure(1)
    p1=plot(t,1 - mean_hf,color,'LineWidth',lwidth);
    h=fill([t, fliplr(t)], [1 - high_hf, fliplr(1 - low_hf)],...
        color, 'facealpha', 0.1, 'edgecolor', 'none');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    % hold on;

    
    figure(2)
    p2=plot(t,1 - mean_h2f,color,'LineWidth',lwidth);
    h=fill([t, fliplr(t)], [1 - high_h2f, fliplr(1 - low_h2f)],...
        color, 'facealpha', 0.1, 'edgecolor', 'none');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    %hold on;
    

    figure(3)
    p3=plot(t,1 - mean_f,color,'LineWidth',lwidth);
    h=fill([t, fliplr(t)], [1 - high_f, fliplr(1 - low_f)],...
        color, 'facealpha', 0.1, 'edgecolor', 'none');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    % hold on;

    
    figure(4)
    p4=plot(t,1 - mean_ih,color,'LineWidth',lwidth);
    h=fill([t, fliplr(t)], [1 - high_ih, fliplr(1 - low_ih)],...
        color, 'facealpha', 0.1, 'edgecolor', 'none');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    

    B = readtable(char(source_cha_time(i)));
    Data2 = B.Variables;
    m_hf_t = Data2(1);
    s_hf_t = Data2(5);
    m_f_t = Data2(2);
    s_f_t = Data2(6);
    m_h2f_t = Data2(3);
    s_h2f_t = Data2(7);
    m_ih_t = Data2(4);
    s_ih_t = Data2(8);
    
    figure(5)
    bar(len_i_char,m_hf_t,bar_width,'Facecolor',[0.94,0.00,0.82])
    %errorbar(len_i_char,m_hf_t,s_hf_t/sqrt(N_hemi(i)),s_hf_t/sqrt(N_hemi(i)),'black','LineWidth',lwidth,'Capsize',capsize)
    errorbar(len_i_char,m_hf_t,m_hf_t/sqrt(N_hemi(i)),m_hf_t/sqrt(N_hemi(i)),'black','LineWidth',lwidth,'Capsize',capsize)

    figure(6)
    % if len_i==4
    %     hf2f_t_charac=0;
    %     hf2f_t_err_neg=0;
    %     hf2f_t_err_pos=0;
    % end
    bar(len_i_char,m_h2f_t,bar_width,'green')
    errorbar(len_i_char,m_h2f_t,m_h2f_t/sqrt(N_fuse(i)),m_h2f_t/sqrt(N_fuse(i)),'black','LineWidth',lwidth,'Capsize',capsize)

    figure(7)
    % if len_i==4
    %     fus_t_charac=0;
    %     fus_t_charac_err_neg=0;
    %     fus_t_charac_err_pos=0;
    % end
    bar(len_i_char,m_f_t,bar_width,'Facecolor',[0.05,0.62,0.31])
    errorbar(len_i_char,m_f_t,m_f_t/sqrt(N_fuse(i)),m_f_t/sqrt(N_fuse(i)),'black','LineWidth',lwidth,'Capsize',capsize)
    
    figure(8)
    % if len_i==4
    %     fus_t_charac=0;
    %     fus_t_charac_err_neg=0;
    %     fus_t_charac_err_pos=0;
    % end
    bar(len_i_char,m_ih_t,bar_width,'red')
    errorbar(len_i_char,m_ih_t,m_ih_t/sqrt(N_fuse(i)),m_ih_t/sqrt(N_fuse(i)),'black','LineWidth',lwidth,'Capsize',capsize)
    
end
plot_arr=[1,2,3,4,5,6,7,8];
%legend(fig1_arr)
%legend(fig2_arr)
%legend(fig3_arr)
for i=1:numel(plot_arr)
    figure(i)
    fig=figure(i);
    %set(fig,'Position',[500,200,800,600])
    ax=gca;
    ax.FontSize=22;
    if i==1
        ylim([-0.05,1.05])
        xlim([0 2.04])
        %ylabel('Hemifusion prob.')
        %xlabel('Time (ms)')
        %legend(gca,'show')
        %legend('3','5','6','7','9')
    elseif i==3
        ylim([-0.05,1.05])
        xlim([0 2.04])
        %ylabel('Fusion probability')
        %xlabel('Time (ms)')
        %legend(gca,'show')
        %legend('3','5','6','7','9')
    elseif i==2
        ylim([-0.05,1.05])
        xlim([0 2.04])
        %ylabel('Hemifus. to fus. prob.','Fontsize',20)
        %xlabel('Time (ms)')
        %legend(gca,'show')
        %legend('3','5','6','7','9')
    elseif i==4
        ylim([-0.05,1.05])
        xlim([0 2.04])
        %ylabel('Irrev. hemifus. prob.','Fontsize',20)
        %xlabel('Time (ms)')
        %legend(gca,'show')
        %legend('3','5','6','7','9')
    elseif i==5
        %ylabel('Hemifus. time (ms)','Fontsize',20)
        %xlabel('Number of SNAREs')
        %ax.XTick="off";
    elseif i==6
        %ylabel({'Hemifus. to', 'fus. time (ms)'},'Fontsize',20)
        %xlabel('Number of SNAREs')
        %ax.XTick="off";
    elseif i==7
        ylim([0,20])
        %ylabel('Fusion time (ms)')
        %xlabel('Number of SNAREs')
        %ax.XTick="off";
    elseif i==8
        ylim([0,20])
        %ylabel('Irrev. hemifus. time (ms)')
        %xlabel('Number of SNAREs')
    end
    %xlim([-0.5,10.9])
    %ylabel('Time (\mus)')
    %xlim([1,15])
    %xticks(snare_len)
    
    ax.FontWeight="bold";
    ax.XMinorTick="off";
    ax.YMinorTick="on";
    ax.Box = 'on';
    ax.TickLength = [0.03 0.025];
    ax.LineWidth = 1;
end

%{
savefig(fig1, 'cdf_hem_nolabel.fig')
savefig(fig2, 'cdf_hf_nolabel.fig')
savefig(fig3, 'cdf_fus_nolabel.fig')
savefig(fig4, 'cdf_ihem_nolabel.fig')
savefig(fig5, 't_hem_nolabel.fig')
savefig(fig6, 't_hf_nolabel.fig')
savefig(fig7, 't_fus_nolabel.fig')
savefig(fig8, 't_ihem_nolabel.fig')
%}