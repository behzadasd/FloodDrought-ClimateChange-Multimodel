%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Discharge Percentile Trend  %%%
%%%         ISI-MIP - data        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Behzad Asadieh   , Ph.D. Candidate                  %%%
%%% Civil Engineering Department - Water Resources      %%%
%%% The City College of The City University of New York %%%
%%% basadie00@citymail.cuny.edu                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
clear;
clc;

load All_Results_hist_rcp8p5

load WetCells_WFD_GPCC_WBM
%feature('UseGenericOpengl', 1); % Sets OpenGL to GenericOpengl to prevent the texts in plots being upside-down (due to a bug in OpenGL)
Impact_Models_Names = {'WBM', 'MacPDM', 'PCR-GLOBWB', 'DBH', 'LPJmL'}; %% , 'H08'};
GCM_Models_Names = {'GFDL-ESM2M', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'MIROC-ESM-CHEM','NorESM1-M'};
Continent_Names = {'Global','N. A.', 'S. A.', 'Europe', 'Oceania', 'Africa', 'Asia', 'India'};

dir_data_in=[pwd '\Results - RCP8.5\']; % Directory to raed raw data from
dir_mat_out=[pwd '\Results - MultiModel\'];
excel_out_name_1=[pwd '\ISI-MIP GCM-GHM Extreme Discharge Change - RCP 8p5 20702099-19712000.xls']; xls_plcs=[5; 10; 15; 20; 25; 30]; %the row in which the results of each GCM will be written

min_NO_st_d=25; % Minimum number of available data for the grid to have a reliable calculation
Lat_n=360; Lon_n=720;

n_mods=25; %% n_mods=numel(GCM_Models_Names) * numel(Impact_Models_Names);

earth_R = 6378; lon_diff_miltiplier = ( Lon_bound(:,2) - Lon_bound(:,1) )' * (pi/180) ; % ( Lon2 - Lon1 ) in the formula
lat_sin_multiplier =  sind(Lat_bound(:,1)) - sind(Lat_bound(:,2)) ;   GridCell_Areas=NaN(size(Lat_bound,1), size(Lon_bound,1));
for ii=1:size(Lat_bound,1)
    for jj=1:size(Lon_bound,1)
        GridCell_Areas (ii,jj) = abs( (earth_R^2) * lon_diff_miltiplier (1,jj) * lat_sin_multiplier (ii,1) );
    end
end
GridCell_Areas=GridCell_Areas / 1e6; GridCell_Areas(isnan(WetCells))=NaN; Area_sum=nansum(nansum(GridCell_Areas));

%%% Relative Change in Multimodel-Average of Discharge Percentiles in RCP8.5 2070-2099 compared to historical 1971-2000 %%%
%%% Normalized relative change: change=(<Q2> - <Q1>)/(<Q2> + <Q1>) - Range is [-1 1]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Multi GHM/GCM Average Maps %%%
dir_out_fig=[pwd '\Figures - Multimodel - RCP8.5\' ]; % Directory to save Figures and Maps
dir_out_name_fig_tag='Disch_ISIMIP_AllModels_Suppl_rcp8p5_hist';
x_limit=[-180 180]; y_limit=[-60 90];
Lon_img=Lon'; Lat_img=Lat;
var_limit=[-0.6 0.6];

red_to_blue=[1.00	0.00	0.00;... % Red to White
    1.00	0.10	0.10;...
    1.00	0.20	0.20;...
    1.00	0.30	0.30;...
    1.00	0.40	0.40;...
    1.00	0.50	0.50;...
    1.00	0.60	0.60;...
    1.00	0.70	0.70;...
    1.00	0.80	0.80;...
    1.00	0.90	0.90;...
    0.90	0.90	1.00;...% White to Blue
    0.80	0.80	1.00;...
    0.70	0.70	1.00;...
    0.60	0.60	1.00;...
    0.50	0.50	1.00;...
    0.40	0.40	1.00;...
    0.30	0.30	1.00;...
    0.20	0.20	1.00;...
    0.10	0.10	1.00;...
    0.00	0.00	1.00];

dir_out_name_fig=[dir_out_fig 'Map_norm_Disch_Change_P95_' dir_out_name_fig_tag];
Variable=Multimodel_Q_Change_P95_hist_rcp8p5_norm;
figure % Columns=GCMs , Rows=GHMs
suptitle('Normalized change in 95th percentile of discharge - All Models - RCP8.5');
set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
for i=1:n_mods
    subplot( numel(Impact_Models_Names) , numel(GCM_Models_Names) , i )
    Var_ij=Variable(:,:,i);
    load coast
    geoshow(flipud(lat), flipud(long),'DisplayType', 'polygon', 'FaceColor', [0.94 0.97 1])
    xlim(x_limit); ylim(y_limit); hold on
    h1=imagesc(Lon_img, Lat_img, Var_ij , var_limit);
    set(h1,'alphadata',~isnan(Var_ij)) % Sets NaN values no color (colors them white)
    set(gcf, 'ColorMap', red_to_blue)
    set(gca,'YDir','normal') % Prevents fliping latitude directions because of the upper arrays being greater than lower arrays
    set(gca,'FontSize',16, 'FontName', 'MS Sens Serif') % Axis Numbers and rages Font
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
    set(gcf,'PaperPositionMode','auto') % The most important line for proper saving the figures ***
end   
%saveas(h1,dir_out_name_fig)
print(dir_out_name_fig,'-dpng','-r300')
close

dir_out_name_fig=[dir_out_fig 'Map_norm_Disch_Change_P05_' dir_out_name_fig_tag];
Variable=Multimodel_Q_Change_P05_hist_rcp8p5_norm;
figure % Columns=GCMs , Rows=GHMs
suptitle('Normalized change in 5th percentile of discharge - All Models - RCP8.5');
set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
for i=1:n_mods
    subplot( numel(Impact_Models_Names) , numel(GCM_Models_Names) , i )
    Var_ij=Variable(:,:,i);
    load coast
    geoshow(flipud(lat), flipud(long),'DisplayType', 'polygon', 'FaceColor', [0.94 0.97 1])
    xlim(x_limit); ylim(y_limit); hold on
    h1=imagesc(Lon_img, Lat_img, Var_ij , var_limit);
    set(h1,'alphadata',~isnan(Var_ij)) % Sets NaN values no color (colors them white)
    set(gcf, 'ColorMap', red_to_blue)
    set(gca,'YDir','normal') % Prevents fliping latitude directions because of the upper arrays being greater than lower arrays
    set(gca,'FontSize',16, 'FontName', 'MS Sens Serif') % Axis Numbers and rages Font
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
    set(gcf,'PaperPositionMode','auto') % The most important line for proper saving the figures ***
end   
%saveas(h1,dir_out_name_fig)
print(dir_out_name_fig,'-dpng','-r300')
close

dir_out_name_fig=[dir_out_fig 'Map_norm_Disch_Change_Median_' dir_out_name_fig_tag];
Variable=Multimodel_Q_Change_Med_hist_rcp8p5_norm;
figure % Columns=GCMs , Rows=GHMs
suptitle('Normalized change in median of discharge - All Models - RCP8.5');
set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
for i=1:n_mods
    subplot( numel(Impact_Models_Names) , numel(GCM_Models_Names) , i )
    Var_ij=Variable(:,:,i);
    load coast
    geoshow(flipud(lat), flipud(long),'DisplayType', 'polygon', 'FaceColor', [0.94 0.97 1])
    xlim(x_limit); ylim(y_limit); hold on
    h1=imagesc(Lon_img, Lat_img, Var_ij , var_limit);
    set(h1,'alphadata',~isnan(Var_ij)) % Sets NaN values no color (colors them white)
    set(gcf, 'ColorMap', red_to_blue)
    set(gca,'YDir','normal') % Prevents fliping latitude directions because of the upper arrays being greater than lower arrays
    set(gca,'FontSize',16, 'FontName', 'MS Sens Serif') % Axis Numbers and rages Font
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
    set(gcf,'PaperPositionMode','auto') % The most important line for proper saving the figures ***
end   
%saveas(h1,dir_out_name_fig)
print(dir_out_name_fig,'-dpng','-r300')
close


toc;
