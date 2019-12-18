function [h]=func_Boxplot_multi_2_Save_new(VarS1, VarS2, VarM1, VarM2, textS1, textS2, textM1, textM2, Var_name,  dir_out_name_fig, y_limit, title_text)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Behzad Asadieh   , Ph.D. Candidate                  %%%
%%% Civil Engineering Department - Water Resources      %%%
%%% The City College of The City University of New York %%%
%%% basadie00@citymail.cuny.edu                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=8; % Number of boxplots per variable

Var_M=NaN(1,n*10); % Combining VarM1 and VarM2 into one variable
for i=1:n
    cc1=(i-1)*10+1;
    Var_M(cc1:cc1+5-1)=VarM1(:,i)';    
    cc2=(i-1)*10+1+5;
    Var_M(cc2:cc2+5-1)=VarM2(:,i)';
    
end
Group_M=NaN(1,n*10);
for i=1:n*2
    cc=(i-1)*5+1;
    Group_M(cc:cc+5-1)=i;
    
end
Position_M=NaN(1,n*2);
for i=1:n
    cc1=(i-1)*2+1;
    Position_M(cc1)=i;
    Position_M(cc1+1)=i+0.3;
end

positions_1=NaN(1,n); positions_1(1,1)=1;
for i=2:n
    positions_1(1,i)=positions_1(1,i-1)+1;
end

positions_2=NaN(1,8); positions_2(1,1)=1.3;
for i=2:n
    positions_2(1,i)=positions_2(1,i-1)+1;
end


%%% BoxPlot %%%
h=figure; hold on;

scatter(positions_1, VarS1, 100 , [0 0.6 0.6] , '*', 'LineWidth', 6);
scatter(positions_2, VarS2, 100 , [0.5 0 0.5] , '*', 'LineWidth', 6);

h1=boxplot(Var_M,Group_M, 'positions', Position_M);
set(h1,'LineWidth',3);
set(gca,'xtick',[mean(Position_M(1:2)) mean(Position_M(3:4)) mean(Position_M(5:6)) mean(Position_M(7:8)) mean(Position_M(9:10)) mean(Position_M(11:12)) mean(Position_M(13:14)) mean(Position_M(15:16))])
set(gca,'xticklabel',Var_name)

color = repmat(['w', 'y'], [1 8]);
h2 = findobj(gca,'Tag','Box');
for j=1:length(h2)
   patch(get(h2(j),'XData'),get(h2(j),'YData'),color(j),'FaceAlpha',.5);
end

legend(textS1, textS2, textM1, textM2,'location','northwest')

scatter(positions_1, VarS1, 100 , [0 0.6 0.6] , '*', 'LineWidth', 6);
scatter(positions_2, VarS2, 100 , [0.5 0 0.5] , '*', 'LineWidth', 6);

ylim(y_limit)
title(title_text)
set(gca,'FontSize',32, 'FontName', 'MS Sens Serif') % Axis Numbers and rages Font
set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window

saveas(h,dir_out_name_fig)
saveas(h,dir_out_name_fig,'png')
% close

end
