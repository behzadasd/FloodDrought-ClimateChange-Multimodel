function [h1]=func_GeoshowMap_Save(Variable, Lat_img, Lon_img, var_limit, x_limit, y_limit, dir_out_name_fig, title_text)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Behzad Asadieh   , Ph.D. Candidate                  %%%
%%% Civil Engineering Department - Water Resources      %%%
%%% The City College of The City University of New York %%%
%%% basadie00@citymail.cuny.edu                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
load coast
geoshow(flipud(lat), flipud(long),'DisplayType', 'polygon', 'FaceColor', [0.94 0.97 1])
%geoshow('landareas.shp', 'FaceColor', [0.97 0.97 0.97]);  %%% Enable this to Change the LandArea color
%states = shaperead('usastatehi', 'UseGeoCoords', true);    %%% Enable this to add USA state boundaries
%geoshow(states,'DefaultFaceColor', 'white', 'DefaultEdgeColor', 'blue', 'FaceColor', 'white');
xlim(x_limit); ylim(y_limit);

hold on

h1=imagesc(Lon_img, Lat_img, Variable , var_limit);
set(h1,'alphadata',~isnan(Variable)) % Sets NaN values no color (colors them white)
%set(gca, 'CLim', [-max(abs(max(max(Variable))),abs(min(min(Variable)))), max(abs(max(max(Variable))),abs(min(min(Variable))))]);
colormap jet
title(title_text)
xlabel('Longitude ºE')
ylabel('Latitude ºN')
set(gca,'YDir','normal') % Prevents fliping latitude directions because of the upper arrays being greater than lower arrays
colorbar('location','Eastoutside')
set(gca,'FontSize',32, 'FontName', 'MS Sens Serif') % Axis Numbers and rages Font
set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
set(gcf,'PaperPositionMode','auto') % The most important line for proper saving the figures ***
%set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 40 20])
saveas(h1,dir_out_name_fig)
print(dir_out_name_fig,'-dpng','-r300')
close


end

