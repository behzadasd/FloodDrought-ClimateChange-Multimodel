function [h1]=func_FloodDroughtRisk_Map_16c(Var_P05, Var_P95, Lat_img, Lon_img, x_limit, y_limit, dir_out_name_fig, title_text)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Behzad Asadieh   , Ph.D. Candidate                  %%%
%%% Civil Engineering Department - Water Resources      %%%
%%% The City College of The City University of New York %%%
%%% basadie00@citymail.cuny.edu                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

colormap_risk = [0.200,	0.000,	0.400;...
    0.600,	0.000,	0.298;...
    0.600,	0.000,	0.000;...
    0.800,	0.000,	0.000;...
    0.600,	0.200,	1.000;...
    1.000,	0.000,	1.000;...
    1.000,	0.000,	0.000;...
    1.000,	0.400,	0.400;...
    0.000,	0.000,	0.600;...
    0.000,	0.502,	1.000;...
    0.000,	1.000,	0.000;...
    0.800,	0.608,	0.600;...
    0.200,	0.200,	1.000;...
    0.600,	1.000,	1.000;...
    0.600,	1.000,	0.800;...
    0.000,	0.400,	0.000];

Lat_n=size(Var_P05,1);
Lon_n=size(Var_P05,2);

Color_Value=NaN(Lat_n, Lon_n);

for lt=1:Lat_n
    for ln=1:Lon_n
        if ~isnan(Var_P05(lt,ln)) && ~isnan(Var_P95(lt,ln))
            
            %%% Flood and Drought risk increase %%%
            if (Var_P95(lt,ln)>=0.5 && Var_P95(lt,ln)<= 1) && (Var_P05(lt,ln)<= -0.5 && Var_P05(lt,ln)>= -1) % Line 1 ColorMap
                Color_Value(lt,ln) = 0.5;
            elseif (Var_P95(lt,ln)>=0.5 && Var_P95(lt,ln)<= 1) && (Var_P05(lt,ln)<= 0 && Var_P05(lt,ln)> -0.5) % Line 2 ColorMap
                Color_Value(lt,ln) = 1.5;
            elseif (Var_P95(lt,ln)>=0 && Var_P95(lt,ln)< 0.5) && (Var_P05(lt,ln)<= -0.5 && Var_P05(lt,ln)>= -1) % Line 5 ColorMap
                Color_Value(lt,ln) = 4.5;
            elseif (Var_P95(lt,ln)>=0 && Var_P95(lt,ln)< 0.5) && (Var_P05(lt,ln)<= 0 && Var_P05(lt,ln)> -0.5) % Line 6 ColorMap
                Color_Value(lt,ln) = 5.5;
            end
            
            %%% Flood risk increase, Drought risk decrease %%%
            if (Var_P95(lt,ln)>=0.5 && Var_P95(lt,ln)<= 1) && (Var_P05(lt,ln)<= 0.5 && Var_P05(lt,ln)> 0) % Line 3 ColorMap
                Color_Value(lt,ln) = 2.5;
            elseif (Var_P95(lt,ln)>=0.5 && Var_P95(lt,ln)<= 1) && (Var_P05(lt,ln)<= 1 && Var_P05(lt,ln)> 0.5) % Line 4 ColorMap
                Color_Value(lt,ln) = 3.5;
            elseif (Var_P95(lt,ln)>=0 && Var_P95(lt,ln)< 0.5) && (Var_P05(lt,ln)<= 0.5 && Var_P05(lt,ln)> 0) % Line 7 ColorMap
                Color_Value(lt,ln) = 6.5;
            elseif (Var_P95(lt,ln)>=0 && Var_P95(lt,ln)< 0.5) && (Var_P05(lt,ln)<= 1 && Var_P05(lt,ln)> 0.5) % Line 8 ColorMap
                Color_Value(lt,ln) = 7.5;
            end
            
            %%% Drought risk increase, Flood risk decrease %%%
            if (Var_P95(lt,ln)>= -0.5 && Var_P95(lt,ln)< 0) && (Var_P05(lt,ln)<= -0.5 && Var_P05(lt,ln)>= -1) % Line 9 ColorMap
                Color_Value(lt,ln) = 8.5;
            elseif (Var_P95(lt,ln)>= -0.5 && Var_P95(lt,ln)< 0) && (Var_P05(lt,ln)<= 0 && Var_P05(lt,ln)> -0.5) % Line 10 ColorMap
                Color_Value(lt,ln) = 9.5;
            elseif (Var_P95(lt,ln)>= -1 && Var_P95(lt,ln)< -0.5) && (Var_P05(lt,ln)<= -0.5 && Var_P05(lt,ln)>= -1) % Line 13 ColorMap
                Color_Value(lt,ln) = 12.5;
            elseif (Var_P95(lt,ln)>= -1 && Var_P95(lt,ln)< -0.5) && (Var_P05(lt,ln)<= 0 && Var_P05(lt,ln)> -0.5) % Line 14 ColorMap
                Color_Value(lt,ln) = 13.5;
            end
            
            %%% Flood and Drought risk decrease %%%
            if (Var_P95(lt,ln)>= -0.5 && Var_P95(lt,ln)< 0) && (Var_P05(lt,ln)<= 0.5 && Var_P05(lt,ln)> 0) % Line 11 ColorMap
                Color_Value(lt,ln) = 10.5;
            elseif (Var_P95(lt,ln)>= -0.5 && Var_P95(lt,ln)< 0) && (Var_P05(lt,ln)<= 1 && Var_P05(lt,ln)> 0.5) % Line 12 ColorMap
                Color_Value(lt,ln) = 11.5;
            elseif (Var_P95(lt,ln)>= -1 && Var_P95(lt,ln)< -0.5) && (Var_P05(lt,ln)<= 0.5 && Var_P05(lt,ln)> 0) % Line 15 ColorMap
                Color_Value(lt,ln) = 14.5;
            elseif (Var_P95(lt,ln)>= -1 && Var_P95(lt,ln)< -0.5) && (Var_P05(lt,ln)<= 1 && Var_P05(lt,ln)> 0.5) % Line 16 ColorMap
                Color_Value(lt,ln) = 15.5;
            end
            
            
        end
    end
end


figure
load coast
geoshow(flipud(lat), flipud(long),'DisplayType', 'polygon', 'FaceColor', [0.94 0.97 1])
%geoshow('landareas.shp', 'FaceColor', [0.97 0.97 0.97]);  %%% Enable this to Change the LandArea color
%states = shaperead('usastatehi', 'UseGeoCoords', true);    %%% Enable this to add USA state boundaries
%geoshow(states,'DefaultFaceColor', 'white', 'DefaultEdgeColor', 'blue', 'FaceColor', 'white');
xlim(x_limit); ylim(y_limit);

hold on

h1=imagesc(Lon_img, Lat_img, Color_Value, [0 16]);
set(h1,'alphadata',~isnan(Color_Value)) % Sets NaN values no color (colors them white)
set(gca,'YDir','normal') % Prevents fliping latitude directions because of the upper arrays being greater than lower arrays
set(gcf, 'ColorMap', colormap_risk)
%colorbar('location','Eastoutside')
title(title_text)
xlabel('Longitude ºE')
ylabel('Latitude ºN')
set(gca,'FontSize',32, 'FontName', 'MS Sens Serif') % Axis Numbers and rages Font
set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
set(gcf,'PaperPositionMode','auto') % The most important line for proper saving the figures ***
%set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 40 20])
saveas(h1,dir_out_name_fig)
print(dir_out_name_fig,'-dpng','-r300')
%close


end

