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

load WetCells_WFD_GPCC_WBM
%feature('UseGenericOpengl', 1); % Sets OpenGL to GenericOpengl to prevent the texts in plots being upside-down (due to a bug in OpenGL)
Impact_Models_Names = {'WBM', 'MacPDM', 'PCR-GLOBWB', 'DBH', 'LPJmL'}; %% , 'H08'};
GCM_Models_Names = {'GFDL-ESM2M', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'MIROC-ESM-CHEM','NorESM1-M'};
Continent_Names = {'Global','N. A.', 'S. A.', 'Europe', 'Oceania', 'Africa', 'Asia', 'India'};

dir_data_in=[pwd '\Results - RCP8.5\']; % Directory to raed raw data from
dir_mat_out=[pwd '\Results - MultiModel\'];
excel_out_name_1=[pwd '\ISI-MIP GCM-GHM Extreme Discharge Change - RCP 8p5 20702099-19712000.xls']; xls_plcs=[5; 10; 15; 20; 25; 30]; %the row in which the results of each GCM will be written
% % excel_out_name_hist=[pwd '\ISI-MIP Extreme Discharge Trend 1971-2000.xls']; xls_plcs=[7; 19; 31; 43; 55]; %the row in which the results of each GCM will be written
% % excel_out_name_rcp8p5=[pwd '\ISI-MIP Extreme Discharge Trend 2070-2099.xls'];

min_NO_st_d=25; % Minimum number of available data for the grid to have a reliable calculation
Lat_n=360; Lon_n=720;
global_land_area=148.94; % Total global land area (km2)

n_mods=25; %% n_mods=numel(GCM_Models_Names) * numel(Impact_Models_Names);

%%% Normalized average Q - All Grids - All of the GCMs and GHMs %%%
Multimodel_Q_ave_Med_hist=NaN(Lat_n, Lon_n, numel(GCM_Models_Names) * numel(Impact_Models_Names)); %%% Stores results of all GCMs and GHMs %%% [ rows= Lat , columns= Lon ] * height= GCMs * GHMs
Multimodel_Q_ave_P05_hist=NaN(Lat_n, Lon_n, numel(GCM_Models_Names) * numel(Impact_Models_Names));
Multimodel_Q_ave_P95_hist=NaN(Lat_n, Lon_n, numel(GCM_Models_Names) * numel(Impact_Models_Names));
Multimodel_Q_ave_Med_rcp8p5=NaN(Lat_n, Lon_n, numel(GCM_Models_Names) * numel(Impact_Models_Names)); %%% Stores results of all GCMs and GHMs %%% [ rows= Lat , columns= Lon ] * height= GCMs * GHMs
Multimodel_Q_ave_P05_rcp8p5=NaN(Lat_n, Lon_n, numel(GCM_Models_Names) * numel(Impact_Models_Names));
Multimodel_Q_ave_P95_rcp8p5=NaN(Lat_n, Lon_n, numel(GCM_Models_Names) * numel(Impact_Models_Names));

%%% Two Sample t-test Significance - All Grids - All of the GCMs and GHMs %%%
Multimodel_T_test_Value_Q_Med_hist_rcp8p5=NaN(Lat_n, Lon_n, numel(GCM_Models_Names) * numel(Impact_Models_Names)); %%% Stores results of all GCMs and GHMs %%% [ rows= Lat , columns= Lon ] * height= GCMs * GHMs
Multimodel_T_test_Value_Q_P05_hist_rcp8p5=NaN(Lat_n, Lon_n, numel(GCM_Models_Names) * numel(Impact_Models_Names));
Multimodel_T_test_Value_Q_P95_hist_rcp8p5=NaN(Lat_n, Lon_n, numel(GCM_Models_Names) * numel(Impact_Models_Names));
Multimodel_T_test_Significance_Q_Med_hist_rcp8p5=NaN(Lat_n, Lon_n, numel(GCM_Models_Names) * numel(Impact_Models_Names)); %%% Stores results of all GCMs and GHMs %%% [ rows= Lat , columns= Lon ] * height= GCMs * GHMs
Multimodel_T_test_Significance_Q_P05_hist_rcp8p5=NaN(Lat_n, Lon_n, numel(GCM_Models_Names) * numel(Impact_Models_Names));
Multimodel_T_test_Significance_Q_P95_hist_rcp8p5=NaN(Lat_n, Lon_n, numel(GCM_Models_Names) * numel(Impact_Models_Names));

%%% Normalized change in Q - All Grids - All of the GCMs and GHMs %%%
Multimodel_Q_Change_Med_hist_rcp8p5_norm=NaN(Lat_n, Lon_n, numel(GCM_Models_Names) * numel(Impact_Models_Names)); %%% Stores results of all GCMs and GHMs %%% [ rows= Lat , columns= Lon ] * height= GCMs * GHMs
Multimodel_Q_Change_P05_hist_rcp8p5_norm=NaN(Lat_n, Lon_n, numel(GCM_Models_Names) * numel(Impact_Models_Names));
Multimodel_Q_Change_P95_hist_rcp8p5_norm=NaN(Lat_n, Lon_n, numel(GCM_Models_Names) * numel(Impact_Models_Names));

%%% Normalized change in Q - Continental Averages - All of the GCMs and GHMs %%%
Multimodel_C_Ave_Q_Change_Med_hist_rcp8p5_norm=NaN(numel(GCM_Models_Names), numel(Continent_Names), numel(Impact_Models_Names)); %%% [ rows= GCMs , columns= Continents ] * height= GHMs
Multimodel_C_Ave_Q_Change_P05_hist_rcp8p5_norm=NaN(numel(GCM_Models_Names), numel(Continent_Names), numel(Impact_Models_Names));
Multimodel_C_Ave_Disch_Change_P95_hist_rcp8p5_norm=NaN(numel(GCM_Models_Names), numel(Continent_Names), numel(Impact_Models_Names));

%%% Normalized change in Q - By Latitude Averages - All of the GCMs and GHMs %%%
Multimodel_Q_Change_Med_hist_rcp8p5_norm_Lat=NaN(numel(GCM_Models_Names), Lat_n, numel(Impact_Models_Names)); %%% [ rows= GCMs , columns= Latitudinal windows ] * height= GHMs
Multimodel_Q_Change_P05_hist_rcp8p5_norm_Lat=NaN(numel(GCM_Models_Names), Lat_n, numel(Impact_Models_Names));
Multimodel_Q_Change_P95_hist_rcp8p5_norm_Lat=NaN(numel(GCM_Models_Names), Lat_n, numel(Impact_Models_Names));

%%% Relative change in Q - Continental Averages - All of the GCMs and GHMs %%%
Multimodel_C_Ave_Q_Change_Med_hist_rcp8p5=NaN(numel(GCM_Models_Names), numel(Continent_Names), numel(Impact_Models_Names)); %%% [ rows= GCMs , columns= Continents ] * height= GHMs
Multimodel_C_Ave_Q_Change_P05_hist_rcp8p5=NaN(numel(GCM_Models_Names), numel(Continent_Names), numel(Impact_Models_Names));
Multimodel_C_Ave_Q_Change_P95_hist_rcp8p5=NaN(numel(GCM_Models_Names), numel(Continent_Names), numel(Impact_Models_Names));

%%% Relative change in Q - Quadrants Averages - All of the GCMs and GHMs %%%
% [Quadrant 1: P95 incrs,P05 incrs]  [Quadrant 2: P95 incrs,P05 decrs]  [Quadrant 3: P95 decrs, P05 incrs] [Quadrant 4: P95 decrs, P05 decrs]
% [columns 1-4] = ave. norm. change in P95 in quadrants 1-4   ,   [columns 5-8] = ave. norm. change in P05 in quadrants 1-4
Multimodel_1_Quadrant_Ave_P95_P05_hist_rcp8p5_norm=NaN(numel(GCM_Models_Names), 8, numel(Impact_Models_Names)); %%% [ rows= GCMs , columns= Quadrants ] * height= GHMs
Multimodel_1_IndDec_Ave_P95_P05_hist_rcp8p5_norm=NaN(numel(GCM_Models_Names), 4, numel(Impact_Models_Names)); %%% [ rows= GCMs , columns= Inc/Dec ] * height= GHMs
Multimodel_1_T_test_Significance_frac_P95_P05_hist_rcp8p5=NaN(numel(GCM_Models_Names), 3, numel(Impact_Models_Names)); %Column 1 = Total land area investigaed, Column 2-3 = Total investigated and total global land area fraction (after excluding the less significant change cells
Multimodel_3_T_test_Significance_frac_P95_P05_hist_rcp8p5=NaN(numel(GCM_Models_Names), 5, numel(Impact_Models_Names)); %Column 1 = Total land area investigaed, Column 2,3= total global land area with increase/decreas in flood, Column 3,4= total global land area with increase/decreas in drought

%%% Relative change in Q - By Latitude Averages - All of the GCMs and GHMs %%%
Multimodel_Q_Change_Med_hist_rcp8p5_Lat=NaN(numel(GCM_Models_Names), Lat_n, numel(Impact_Models_Names)); %%% [ rows= GCMs , columns= Latitudinal windows ] * height= GHMs
Multimodel_Q_Change_P05_hist_rcp8p5_Lat=NaN(numel(GCM_Models_Names), Lat_n, numel(Impact_Models_Names));
Multimodel_Q_Change_P95_hist_rcp8p5_Lat=NaN(numel(GCM_Models_Names), Lat_n, numel(Impact_Models_Names));

ss=0;
%for ghm_i=1:size(Impact_Models_Names,2) % For the total number of Imapc Models
for ghm_i=1:5
    
    Impact_Models_Name = Impact_Models_Names{1, ghm_i};
    
    file_name_dir_in=[dir_data_in 'Results_Disch_PercentileTrend_ISIMIP_' Impact_Models_Name '_hist_rcp8p5_1971-2000_2070-2099' '.mat']; % Directory and Name of .nc file to be loaded
    load (file_name_dir_in) % Loading the Data
    
    earth_R = 6378; lon_diff_miltiplier = ( Lon_bound(:,2) - Lon_bound(:,1) )' * (pi/180) ; % ( Lon2 - Lon1 ) in the formula
    lat_sin_multiplier =  sind(Lat_bound(:,1)) - sind(Lat_bound(:,2)) ;   GridCell_Areas=NaN(size(Lat_bound,1), size(Lon_bound,1));
    for ii=1:size(Lat_bound,1)
        for jj=1:size(Lon_bound,1)
            GridCell_Areas (ii,jj) = abs( (earth_R^2) * lon_diff_miltiplier (1,jj) * lat_sin_multiplier (ii,1) );
        end
    end
    GridCell_Areas=GridCell_Areas / 1e6; GridCell_Areas(isnan(GHM_Disch_ave_Median_hist(:,:,1)))=NaN; Area_sum=nansum(nansum(GridCell_Areas));
    
    Multimodel_Q_ave_Med_hist (:,:, ss+1:ss+numel(GCM_Models_Names)) = GHM_Disch_ave_Median_hist;
    Multimodel_Q_ave_P05_hist (:,:, ss+1:ss+numel(GCM_Models_Names)) = GHM_Disch_ave_P05_hist;
    Multimodel_Q_ave_P95_hist (:,:, ss+1:ss+numel(GCM_Models_Names)) = GHM_Disch_ave_P95_hist;
    Multimodel_Q_ave_Med_rcp8p5 (:,:, ss+1:ss+numel(GCM_Models_Names)) = GHM_Disch_ave_Median_rcp8p5;
    Multimodel_Q_ave_P05_rcp8p5 (:,:, ss+1:ss+numel(GCM_Models_Names)) = GHM_Disch_ave_P05_rcp8p5;
    Multimodel_Q_ave_P95_rcp8p5 (:,:, ss+1:ss+numel(GCM_Models_Names)) = GHM_Disch_ave_P95_rcp8p5;
    
    Multimodel_T_test_Value_Q_Med_hist_rcp8p5 (:,:, ss+1:ss+numel(GCM_Models_Names)) = GHM_T_test_Value_Disch_Med_rcp8p5_hist;
    Multimodel_T_test_Value_Q_P05_hist_rcp8p5 (:,:, ss+1:ss+numel(GCM_Models_Names)) = GHM_T_test_Value_Disch_P95_rcp8p5_hist;
    Multimodel_T_test_Value_Q_P95_hist_rcp8p5 (:,:, ss+1:ss+numel(GCM_Models_Names)) = GHM_T_test_Value_Disch_P95_rcp8p5_hist;
    Multimodel_T_test_Significance_Q_Med_hist_rcp8p5 (:,:, ss+1:ss+numel(GCM_Models_Names)) = GHM_T_test_Significance_Disch_Med_rcp8p5_hist;
    Multimodel_T_test_Significance_Q_P05_hist_rcp8p5 (:,:, ss+1:ss+numel(GCM_Models_Names)) = GHM_T_test_Significance_Disch_P05_rcp8p5_hist;
    Multimodel_T_test_Significance_Q_P95_hist_rcp8p5 (:,:, ss+1:ss+numel(GCM_Models_Names)) = GHM_T_test_Significance_Disch_P95_rcp8p5_hist;
    
    %%% Relative Change in Average of Discharge Percentiles in RCP8.5 2070-2099 compared to historical 1971-2000 %%%
    %%% Original relative change: change=(Q2-Q1)/Q1 * 100
    GHM_Disch_Change_Median_rcp8p5_hist=NaN(Lat_n, Lon_n, numel(GCM_Models_Names));
    GHM_Disch_Change_P05_rcp8p5_hist=NaN(Lat_n, Lon_n, numel(GCM_Models_Names));
    GHM_Disch_Change_P95_rcp8p5_hist=NaN(Lat_n, Lon_n, numel(GCM_Models_Names));
    for gcm_i=1:5
        
        for ii=1:Lat_n
            for jj=1:Lon_n
                if ~isnan(GHM_Disch_ave_Median_rcp8p5 (ii,jj,gcm_i)) && ~isnan(GHM_Disch_ave_Median_hist (ii,jj,gcm_i)) && GHM_Disch_ave_Median_hist (ii,jj,gcm_i)~=0
                    GHM_Disch_Change_Median_rcp8p5_hist(ii,jj,gcm_i) = (( GHM_Disch_ave_Median_rcp8p5 (ii,jj,gcm_i) - GHM_Disch_ave_Median_hist (ii,jj,gcm_i) ) / GHM_Disch_ave_Median_hist (ii,jj,gcm_i) ) *100;
                end
                if ~isnan(GHM_Disch_ave_P05_rcp8p5 (ii,jj,gcm_i)) && ~isnan(GHM_Disch_ave_P05_hist (ii,jj,gcm_i)) && GHM_Disch_ave_P05_hist (ii,jj,gcm_i)~=0
                    GHM_Disch_Change_P05_rcp8p5_hist(ii,jj,gcm_i) = (( GHM_Disch_ave_P05_rcp8p5 (ii,jj,gcm_i) - GHM_Disch_ave_P05_hist (ii,jj,gcm_i) ) / GHM_Disch_ave_P05_hist (ii,jj,gcm_i) ) *100;
                end
                if ~isnan(GHM_Disch_ave_P95_rcp8p5 (ii,jj,gcm_i)) && ~isnan(GHM_Disch_ave_P95_hist (ii,jj,gcm_i)) && GHM_Disch_ave_P95_hist (ii,jj,gcm_i)~=0
                    GHM_Disch_Change_P95_rcp8p5_hist(ii,jj,gcm_i) = (( GHM_Disch_ave_P95_rcp8p5 (ii,jj,gcm_i) - GHM_Disch_ave_P95_hist (ii,jj,gcm_i) ) / GHM_Disch_ave_P95_hist (ii,jj,gcm_i) ) *100;
                end
                
            end
        end
        
    end
    
    %%% Relative Change in Average of Discharge Percentiles in RCP8.5 2070-2099 compared to historical 1971-2000 %%%
    %%% Normalized relative change: change=(Q2-Q1)/(Q2+Q1) - Range is [-1 1]
    GHM_Disch_Change_Median_rcp8p5_hist_norm=NaN(Lat_n, Lon_n, numel(GCM_Models_Names));
    GHM_Disch_Change_P05_rcp8p5_hist_norm=NaN(Lat_n, Lon_n, numel(GCM_Models_Names));
    GHM_Disch_Change_P95_rcp8p5_hist_norm=NaN(Lat_n, Lon_n, numel(GCM_Models_Names));
    for gcm_i=1:5
        
        for ii=1:Lat_n
            for jj=1:Lon_n
                if ~isnan(GHM_Disch_ave_Median_rcp8p5 (ii,jj,gcm_i)) && ~isnan(GHM_Disch_ave_Median_hist (ii,jj,gcm_i)) && GHM_Disch_ave_Median_hist (ii,jj,gcm_i)~=0
                    GHM_Disch_Change_Median_rcp8p5_hist_norm(ii,jj,gcm_i) = (( GHM_Disch_ave_Median_rcp8p5(ii,jj,gcm_i) - GHM_Disch_ave_Median_hist(ii,jj,gcm_i) ) / ( GHM_Disch_ave_Median_rcp8p5(ii,jj,gcm_i) + GHM_Disch_ave_Median_hist(ii,jj,gcm_i) ) );
                end
                if ~isnan(GHM_Disch_ave_P05_rcp8p5(ii,jj,gcm_i)) && ~isnan(GHM_Disch_ave_P05_hist(ii,jj,gcm_i)) && GHM_Disch_ave_P05_hist(ii,jj,gcm_i)~=0
                    GHM_Disch_Change_P05_rcp8p5_hist_norm(ii,jj,gcm_i) = (( GHM_Disch_ave_P05_rcp8p5(ii,jj,gcm_i) - GHM_Disch_ave_P05_hist (ii,jj,gcm_i) ) / ( GHM_Disch_ave_P05_rcp8p5(ii,jj,gcm_i) + GHM_Disch_ave_P05_hist (ii,jj,gcm_i) ) );
                end
                if ~isnan(GHM_Disch_ave_P95_rcp8p5(ii,jj,gcm_i)) && ~isnan(GHM_Disch_ave_P95_hist(ii,jj,gcm_i)) && GHM_Disch_ave_P95_hist(ii,jj,gcm_i)~=0
                    GHM_Disch_Change_P95_rcp8p5_hist_norm(ii,jj,gcm_i) = (( GHM_Disch_ave_P95_rcp8p5(ii,jj,gcm_i) - GHM_Disch_ave_P95_hist(ii,jj,gcm_i) ) / ( GHM_Disch_ave_P95_rcp8p5(ii,jj,gcm_i) + GHM_Disch_ave_P95_hist(ii,jj,gcm_i) ) );
                end
                
            end
        end
        
    end
    
    Multimodel_Q_Change_Med_hist_rcp8p5_norm (:,:, ss+1:ss+numel(GCM_Models_Names)) = GHM_Disch_Change_Median_rcp8p5_hist_norm;
    Multimodel_Q_Change_P05_hist_rcp8p5_norm (:,:, ss+1:ss+numel(GCM_Models_Names)) = GHM_Disch_Change_P05_rcp8p5_hist_norm;
    Multimodel_Q_Change_P95_hist_rcp8p5_norm (:,:, ss+1:ss+numel(GCM_Models_Names)) = GHM_Disch_Change_P95_rcp8p5_hist_norm;
    ss=ss+numel(GCM_Models_Names);
    
    %%% Averaging the results Globally and by Continent %%%
    for gcm_i=1:numel(GCM_Models_Names)
        [C_Ave_Median, C_Ave_P05, C_Ave_P95]=func_TrendAve_CellAreaWeighted_3var_grid(GHM_Disch_Change_Median_rcp8p5_hist_norm(:,:,gcm_i), GHM_Disch_Change_P05_rcp8p5_hist_norm(:,:,gcm_i), GHM_Disch_Change_P95_rcp8p5_hist_norm(:,:,gcm_i), Lat, Lon, Lat_bound, Lon_bound);
        Multimodel_C_Ave_Q_Change_Med_hist_rcp8p5_norm(gcm_i, :, ghm_i)=C_Ave_Median';
        Multimodel_C_Ave_Q_Change_P05_hist_rcp8p5_norm(gcm_i, :, ghm_i)=C_Ave_P05';
        Multimodel_C_Ave_Disch_Change_P95_hist_rcp8p5_norm(gcm_i, :, ghm_i)=C_Ave_P95';
        
        [C_Ave_Median, C_Ave_P05, C_Ave_P95]=func_TrendAve_CellAreaWeighted_3var_grid(GHM_Disch_Change_Median_rcp8p5_hist(:,:,gcm_i), GHM_Disch_Change_P05_rcp8p5_hist(:,:,gcm_i), GHM_Disch_Change_P95_rcp8p5_hist(:,:,gcm_i), Lat, Lon, Lat_bound, Lon_bound);
        Multimodel_C_Ave_Q_Change_Med_hist_rcp8p5(gcm_i, :, ghm_i)=C_Ave_Median';
        Multimodel_C_Ave_Q_Change_P05_hist_rcp8p5(gcm_i, :, ghm_i)=C_Ave_P05';
        Multimodel_C_Ave_Q_Change_P95_hist_rcp8p5(gcm_i, :, ghm_i)=C_Ave_P95';
    end
    
    %%% Writing the results in excell file %%%
    xlswrite(excel_out_name_1, Multimodel_C_Ave_Q_Change_P05_hist_rcp8p5_norm(:, 1:7, ghm_i), 'Multimodel_Continents', ['D' num2str( (ghm_i) * 5 ) ':' 'J' num2str( ((ghm_i) * 5) + 4 )])
    xlswrite(excel_out_name_1, Multimodel_C_Ave_Disch_Change_P95_hist_rcp8p5_norm(:, 1:7, ghm_i), 'Multimodel_Continents', ['L' num2str( (ghm_i) * 5 ) ':' 'R' num2str( ((ghm_i) * 5) + 4 )])
    xlswrite(excel_out_name_1, Multimodel_C_Ave_Q_Change_Med_hist_rcp8p5_norm(:, 1:7, ghm_i), 'Multimodel_Continents', ['T' num2str( (ghm_i) * 5 ) ':' 'Z' num2str( ((ghm_i) * 5) + 4 )])
    
    %%% Averaging the results by latitude %%%
    Multimodel_Q_Change_Med_hist_rcp8p5_norm_Lat(:, :, ghm_i)=( reshape(nanmean(GHM_Disch_Change_Median_rcp8p5_hist_norm,2), [Lat_n, numel(GCM_Models_Names)]) )';
    Multimodel_Q_Change_P05_hist_rcp8p5_norm_Lat(:, :, ghm_i)=( reshape(nanmean(GHM_Disch_Change_P05_rcp8p5_hist_norm,2), [Lat_n, numel(GCM_Models_Names)]) )';
    Multimodel_Q_Change_P95_hist_rcp8p5_norm_Lat(:, :, ghm_i)=( reshape(nanmean(GHM_Disch_Change_P95_rcp8p5_hist_norm,2), [Lat_n, numel(GCM_Models_Names)]) )';
    
    Multimodel_Q_Change_Med_hist_rcp8p5_Lat(:, :, ghm_i)=( reshape(nanmean(GHM_Disch_Change_Median_rcp8p5_hist,2), [Lat_n, numel(GCM_Models_Names)]) )';
    Multimodel_Q_Change_P05_hist_rcp8p5_Lat(:, :, ghm_i)=( reshape(nanmean(GHM_Disch_Change_P05_rcp8p5_hist,2), [Lat_n, numel(GCM_Models_Names)]) )';
    Multimodel_Q_Change_P95_hist_rcp8p5_Lat(:, :, ghm_i)=( reshape(nanmean(GHM_Disch_Change_P95_rcp8p5_hist,2), [Lat_n, numel(GCM_Models_Names)]) )';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Normalized Change Averaged by Quadrants %%%
    for gcm_i=1:5
        %[Quadrant 1: P95 incrs,P05 incrs]  [Quadrant 2: P95 incrs,P05 decrs]  [Quadrant 3: P95 decrs, P05 incrs] [Quadrant 4: P95 decrs, P05 decrs]
        % Multimodel_Quadrant_Ave_P95_P05 : [columns 1-4] = ave. norm. change in P95 in quadrants 1-4   ,   [columns 5-8] = ave. norm. change in P05 in quadrants 1-4
        Var_P05=GHM_Disch_Change_P05_rcp8p5_hist_norm(:,:,gcm_i);
        Var_P95=GHM_Disch_Change_P95_rcp8p5_hist_norm(:,:,gcm_i);
        [Flood_Drought_Indicator_All]=func_FloodDroughtRisk_indicator(Var_P05, Var_P95);
        
        %%%% All cells %%%
        Var_P05_help=Var_P05; Var_P95_help=Var_P95; GridCell_Areas_help=GridCell_Areas; % The Quadrant experiencing increased Flood and increased Drought (indicator 1)
        Var_P05_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==2 | Flood_Drought_Indicator_All==3 | Flood_Drought_Indicator_All==4)=NaN; % excluding the cells that do not fit in this indicator
        Var_P95_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==2 | Flood_Drought_Indicator_All==3 | Flood_Drought_Indicator_All==4)=NaN;
        GridCell_Areas_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==2 | Flood_Drought_Indicator_All==3 | Flood_Drought_Indicator_All==4)=NaN;
        Multimodel_1_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(gcm_i, 1, ghm_i)=nansum(nansum( Var_P95_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help )); % Cell-Area-Weighted averaging
        Multimodel_1_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(gcm_i, 5, ghm_i)=nansum(nansum( Var_P05_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help ));
        
        Var_P05_help=Var_P05; Var_P95_help=Var_P95; GridCell_Areas_help=GridCell_Areas; % The Quadrant experiencing increased Flood and decresed Drought (indicator 2)
        Var_P05_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==1 | Flood_Drought_Indicator_All==3 | Flood_Drought_Indicator_All==4)=NaN; % excluding the cells that do not fit in this indicator
        Var_P95_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==1 | Flood_Drought_Indicator_All==3 | Flood_Drought_Indicator_All==4)=NaN;
        GridCell_Areas_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==1 | Flood_Drought_Indicator_All==3 | Flood_Drought_Indicator_All==4)=NaN;
        Multimodel_1_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(gcm_i, 2, ghm_i)=nansum(nansum( Var_P95_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help )); % Cell-Area-Weighted averaging
        Multimodel_1_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(gcm_i, 6, ghm_i)=nansum(nansum( Var_P05_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help ));
        
        Var_P05_help=Var_P05; Var_P95_help=Var_P95; GridCell_Areas_help=GridCell_Areas; % The Quadrant experiencing decresed Flood and increased Drought (indicator 3)
        Var_P05_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==1 | Flood_Drought_Indicator_All==2 | Flood_Drought_Indicator_All==4)=NaN; % excluding the cells that do not fit in this indicator
        Var_P95_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==1 | Flood_Drought_Indicator_All==2 | Flood_Drought_Indicator_All==4)=NaN;
        GridCell_Areas_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==1 | Flood_Drought_Indicator_All==2 | Flood_Drought_Indicator_All==4)=NaN;
        Multimodel_1_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(gcm_i, 3, ghm_i)=nansum(nansum( Var_P95_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help )); % Cell-Area-Weighted averaging
        Multimodel_1_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(gcm_i, 7, ghm_i)=nansum(nansum( Var_P05_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help ));
        
        Var_P05_help=Var_P05; Var_P95_help=Var_P95; GridCell_Areas_help=GridCell_Areas; % The Quadrant experiencing decresed Flood and decresed Drought (indicator 4)
        Var_P05_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==1 | Flood_Drought_Indicator_All==2 | Flood_Drought_Indicator_All==3)=NaN; % excluding the cells that do not fit in this indicator
        Var_P95_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==1 | Flood_Drought_Indicator_All==2 | Flood_Drought_Indicator_All==3)=NaN;
        GridCell_Areas_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==1 | Flood_Drought_Indicator_All==2 | Flood_Drought_Indicator_All==3)=NaN;
        Multimodel_1_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(gcm_i, 4, ghm_i)=nansum(nansum( Var_P95_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help )); % Cell-Area-Weighted averaging
        Multimodel_1_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(gcm_i, 8, ghm_i)=nansum(nansum( Var_P05_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help ));
        
        % All Quadrants (To calculate the total land area that is investigated, after excluding the cells with less change significance than the defined threshold
        Multimodel_1_T_test_Significance_frac_P95_P05_hist_rcp8p5(gcm_i, 1, ghm_i) = Area_sum ;
        Var_P95_help=GHM_T_test_Significance_Disch_P95_rcp8p5_hist(:,:,gcm_i); GridCell_Areas_help=GridCell_Areas;
        GridCell_Areas_help(isnan(Var_P95_help) | Var_P95_help==0 )=NaN;
        Multimodel_1_T_test_Significance_frac_P95_P05_hist_rcp8p5(gcm_i, 2, ghm_i)=nansum(nansum( GridCell_Areas_help )) / Area_sum ;
        Var_P05_help=GHM_T_test_Significance_Disch_P05_rcp8p5_hist(:,:,gcm_i); GridCell_Areas_help=GridCell_Areas;
        GridCell_Areas_help(isnan(Var_P05_help) | Var_P05_help==0 )=NaN;
        Multimodel_1_T_test_Significance_frac_P95_P05_hist_rcp8p5(gcm_i, 3, ghm_i)=nansum(nansum( GridCell_Areas_help )) / Area_sum ;
        
    end
    
    %%% Writing the results in excell file %%%
    xlswrite(excel_out_name_1, Multimodel_1_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(:,1,ghm_i), 'Multimodel_Quadrants', ['E' num2str( (ghm_i) * 5 ) ':' 'E' num2str( ((ghm_i) * 5) + 4 )] )
    xlswrite(excel_out_name_1, Multimodel_1_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(:,2,ghm_i), 'Multimodel_Quadrants', ['D' num2str( (ghm_i) * 5 ) ':' 'D' num2str( ((ghm_i) * 5) + 4 )] )
    xlswrite(excel_out_name_1, Multimodel_1_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(:,3,ghm_i), 'Multimodel_Quadrants', ['F' num2str( (ghm_i) * 5 ) ':' 'F' num2str( ((ghm_i) * 5) + 4 )] )
    xlswrite(excel_out_name_1, Multimodel_1_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(:,4,ghm_i), 'Multimodel_Quadrants', ['G' num2str( (ghm_i) * 5 ) ':' 'G' num2str( ((ghm_i) * 5) + 4 )] )
    xlswrite(excel_out_name_1, Multimodel_1_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(:,5,ghm_i) * -1, 'Multimodel_Quadrants', ['J' num2str( (ghm_i) * 5 ) ':' 'J' num2str( ((ghm_i) * 5) + 4 )] )
    xlswrite(excel_out_name_1, Multimodel_1_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(:,6,ghm_i) * -1, 'Multimodel_Quadrants', ['I' num2str( (ghm_i) * 5 ) ':' 'I' num2str( ((ghm_i) * 5) + 4 )] )
    xlswrite(excel_out_name_1, Multimodel_1_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(:,7,ghm_i) * -1, 'Multimodel_Quadrants', ['K' num2str( (ghm_i) * 5 ) ':' 'K' num2str( ((ghm_i) * 5) + 4 )] )
    xlswrite(excel_out_name_1, Multimodel_1_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(:,8,ghm_i) * -1, 'Multimodel_Quadrants', ['L' num2str( (ghm_i) * 5 ) ':' 'L' num2str( ((ghm_i) * 5) + 4 )] )
    xlswrite(excel_out_name_1, Multimodel_1_T_test_Significance_frac_P95_P05_hist_rcp8p5(:,1,ghm_i), 'Multimodel_Quadrants', ['N' num2str( (ghm_i) * 5 ) ':' 'N' num2str( ((ghm_i) * 5) + 4 )] )
    xlswrite(excel_out_name_1, Multimodel_1_T_test_Significance_frac_P95_P05_hist_rcp8p5(:,2,ghm_i), 'Multimodel_Quadrants', ['O' num2str( (ghm_i) * 5 ) ':' 'O' num2str( ((ghm_i) * 5) + 4 )] )
    xlswrite(excel_out_name_1, Multimodel_1_T_test_Significance_frac_P95_P05_hist_rcp8p5(:,3,ghm_i), 'Multimodel_Quadrants', ['P' num2str( (ghm_i) * 5 ) ':' 'P' num2str( ((ghm_i) * 5) + 4 )] )
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Normalized Change Averaged by Increase/Decrease %%%
    for gcm_i=1:5
        %[Column 1: P95 incrs]  [Column 2: P95 decrs]  [Column 3: P05 incrs]  [Column 4: P05 decrs]
        % Multimodel_IndDec_Ave_P95_P05 : column 1 = P95 change average in increasing cells, column 2 = P95 change average in decreasing cells, column 3 = P05 change average in increasing cells, column 4 = P05 change average in decreasing cells
        Var_P05=GHM_Disch_Change_P05_rcp8p5_hist_norm(:,:,gcm_i);
        Var_P95=GHM_Disch_Change_P95_rcp8p5_hist_norm(:,:,gcm_i);
        
        %%%% All cells %%%
        Var_P95_help=Var_P95; GridCell_Areas_help=GridCell_Areas; % Cells with Increased Flood (Quads 1 and 2)(Indicators 1 and 2)
        Var_P95_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==3 | Flood_Drought_Indicator_All==4)=NaN; % excluding the cells that do not fit in this indicators
        GridCell_Areas_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==3 | Flood_Drought_Indicator_All==4)=NaN;
        Multimodel_1_IndDec_Ave_P95_P05_hist_rcp8p5_norm(gcm_i, 1, ghm_i)=nansum(nansum( Var_P95_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help )); % Cell-Area-Weighted averaging
        
        Var_P95_help=Var_P95; GridCell_Areas_help=GridCell_Areas; % Cells with Decreased Flood (Quads 3 and 4)(Indicators 3 and 4)
        Var_P95_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==1 | Flood_Drought_Indicator_All==2)=NaN; % excluding the cells that do not fit in this indicators
        GridCell_Areas_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==1 | Flood_Drought_Indicator_All==2)=NaN;
        Multimodel_1_IndDec_Ave_P95_P05_hist_rcp8p5_norm(gcm_i, 2, ghm_i)=nansum(nansum( Var_P95_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help )); % Cell-Area-Weighted averaging
        
        Var_P05_help=Var_P05;GridCell_Areas_help=GridCell_Areas; % Cells with Increased Drought (Quads 1 and 3)(Indicators 1 and 3)
        Var_P05_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==2 | Flood_Drought_Indicator_All==4)=NaN; % excluding the cells that do not fit in this indicator
        GridCell_Areas_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==2 | Flood_Drought_Indicator_All==4)=NaN;
        Multimodel_1_IndDec_Ave_P95_P05_hist_rcp8p5_norm(gcm_i, 3, ghm_i)=nansum(nansum( Var_P05_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help )); % Cell-Area-Weighted averaging
        
        Var_P05_help=Var_P05; GridCell_Areas_help=GridCell_Areas; % Cells with Decreased Drought (Quads 2 and 4)(Indicators 2 and 4)
        Var_P05_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==1 | Flood_Drought_Indicator_All==3)=NaN; % excluding the cells that do not fit in this indicator
        GridCell_Areas_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==1 | Flood_Drought_Indicator_All==3)=NaN;
        Multimodel_1_IndDec_Ave_P95_P05_hist_rcp8p5_norm(gcm_i, 4, ghm_i)=nansum(nansum( Var_P05_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help )); % Cell-Area-Weighted averaging
        
        % All Quadrants (To calculate the total global land area, after excluding the cells with less change significance than the defined threshold
        Multimodel_3_T_test_Significance_frac_P95_P05_hist_rcp8p5(gcm_i, 1, ghm_i) = Area_sum ; % Column 1=Sum of invesitigated area
        Var_P95_help=GHM_T_test_Significance_Disch_P95_rcp8p5_hist(:,:,gcm_i); GridCell_Areas_help=GridCell_Areas;
        GridCell_Areas_help(isnan(Var_P95_help) | Var_P95_help==0 | Var_P95_help==-1)=NaN;
        Multimodel_3_T_test_Significance_frac_P95_P05_hist_rcp8p5(gcm_i, 2, ghm_i)=nansum(nansum( GridCell_Areas_help )) / global_land_area ; % Column 2= % of global land area with significant increased flood
        Var_P95_help=GHM_T_test_Significance_Disch_P95_rcp8p5_hist(:,:,gcm_i); GridCell_Areas_help=GridCell_Areas;
        GridCell_Areas_help(isnan(Var_P95_help) | Var_P95_help==0 | Var_P95_help==1)=NaN;
        Multimodel_3_T_test_Significance_frac_P95_P05_hist_rcp8p5(gcm_i, 3, ghm_i)=nansum(nansum( GridCell_Areas_help )) / global_land_area ; % Column 3= % of global land area with significant decreased flood
        Var_P05_help=GHM_T_test_Significance_Disch_P05_rcp8p5_hist(:,:,gcm_i); GridCell_Areas_help=GridCell_Areas;
        GridCell_Areas_help(isnan(Var_P05_help) | Var_P05_help==0 | Var_P05_help==-1)=NaN;
        Multimodel_3_T_test_Significance_frac_P95_P05_hist_rcp8p5(gcm_i, 4, ghm_i)=nansum(nansum( GridCell_Areas_help )) / global_land_area ; % Column 4= % of global land area with significant increased drought
        Var_P05_help=GHM_T_test_Significance_Disch_P05_rcp8p5_hist(:,:,gcm_i); GridCell_Areas_help=GridCell_Areas;
        GridCell_Areas_help(isnan(Var_P05_help) | Var_P05_help==0 | Var_P05_help==1)=NaN;
        Multimodel_3_T_test_Significance_frac_P95_P05_hist_rcp8p5(gcm_i, 5, ghm_i)=nansum(nansum( GridCell_Areas_help )) / global_land_area ; % Column 4= % of global land area with significant decreased drought
        
    end
    
    %%% Writing the results in excell file %%%
    xlswrite(excel_out_name_1, Multimodel_1_IndDec_Ave_P95_P05_hist_rcp8p5_norm(:,1,ghm_i), 'Multimodel_IncDec', ['D' num2str( (ghm_i) * 5 ) ':' 'D' num2str( ((ghm_i) * 5) + 4 )] )
    xlswrite(excel_out_name_1, Multimodel_1_IndDec_Ave_P95_P05_hist_rcp8p5_norm(:,2,ghm_i), 'Multimodel_IncDec', ['E' num2str( (ghm_i) * 5 ) ':' 'E' num2str( ((ghm_i) * 5) + 4 )] )
    xlswrite(excel_out_name_1, Multimodel_1_IndDec_Ave_P95_P05_hist_rcp8p5_norm(:,3,ghm_i) * -1, 'Multimodel_IncDec', ['F' num2str( (ghm_i) * 5 ) ':' 'F' num2str( ((ghm_i) * 5) + 4 )] )
    xlswrite(excel_out_name_1, Multimodel_1_IndDec_Ave_P95_P05_hist_rcp8p5_norm(:,4,ghm_i) * -1, 'Multimodel_IncDec', ['G' num2str( (ghm_i) * 5 ) ':' 'G' num2str( ((ghm_i) * 5) + 4 )] )
    xlswrite(excel_out_name_1, Multimodel_3_T_test_Significance_frac_P95_P05_hist_rcp8p5(:,1,ghm_i), 'Multimodel_IncDec', ['I' num2str( (ghm_i) * 5 ) ':' 'I' num2str( ((ghm_i) * 5) + 4 )] )
    xlswrite(excel_out_name_1, Multimodel_3_T_test_Significance_frac_P95_P05_hist_rcp8p5(:,2,ghm_i), 'Multimodel_IncDec', ['K' num2str( (ghm_i) * 5 ) ':' 'K' num2str( ((ghm_i) * 5) + 4 )] )
    xlswrite(excel_out_name_1, Multimodel_3_T_test_Significance_frac_P95_P05_hist_rcp8p5(:,3,ghm_i), 'Multimodel_IncDec', ['L' num2str( (ghm_i) * 5 ) ':' 'L' num2str( ((ghm_i) * 5) + 4 )] )
    xlswrite(excel_out_name_1, Multimodel_3_T_test_Significance_frac_P95_P05_hist_rcp8p5(:,4,ghm_i), 'Multimodel_IncDec', ['N' num2str( (ghm_i) * 5 ) ':' 'N' num2str( ((ghm_i) * 5) + 4 )] )
    xlswrite(excel_out_name_1, Multimodel_3_T_test_Significance_frac_P95_P05_hist_rcp8p5(:,5,ghm_i), 'Multimodel_IncDec', ['O' num2str( (ghm_i) * 5 ) ':' 'O' num2str( ((ghm_i) * 5) + 4 )] )
  
    % %     %%%%%%%%%%%%
% %     %%% Maps %%%
% %     %%%%%%%%%%%%
% %     dir_out_fig=[pwd '\Figures - Global MultiGCM - RCP8.5\' Impact_Models_Name '\']; % Directory to save Figures and Maps
% %     dir_out_name_fig_tag=['ISIMIP_' Impact_Models_Name '_MultiGCM_rcp8p5_hist'];
% %     x_limit=[-180 180]; y_limit=[-75 90];
% %     
% %     Lon_img=Lon'; Lat_img=Lat;
% %     
% %     %%% Average by Latitude %%%
% %     Variable=Multimodel_Q_Change_Med_hist_rcp8p5_norm_Lat(:, :, ghm_i);
% %     dir_out_name_fig=[dir_out_fig 'Lat_Disch_Change_Median_' dir_out_name_fig_tag];
% %     h=figure;
% %     hold on
% %     func_shadedErrorBar(Lat',Variable,{@mean,@std},{'g','LineWidth',4});
% %     title({'Latitudinal average of normalized change in median of discharge' ; ['ISI-MIP - GHM: ' Impact_Models_Name ' , GCM: multiple - [2070-2099 RCP8.5] compared to [1971-2000]']})
% %     legend('GCM Average +/- St.Dev.', 'location', 'NorthWest');
% %     xlabel('Latitude ºN'); ylabel('Average of normalized change in discharge');
% %     set(gca,'FontSize',36, 'FontName', 'MS Sens Serif') % Axis Numbers and ranges Font
% %     set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
% %     set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
% %     xlim([-60,85])
% %     ylim([-1,1])
% %     plot([-90,90], [0,0],'LineWidth',0.5,'Color','black');
% %     set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
% %     set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
% %     set(gcf,'PaperPositionMode','auto') % The most important line for proper saving the figures ***
% %     saveas(h,dir_out_name_fig)
% %     print(dir_out_name_fig,'-dpng','-r300')
% %     close
% %     
% %     Variable=Multimodel_Q_Change_P05_hist_rcp8p5_norm_Lat(:, :, ghm_i);
% %     dir_out_name_fig=[dir_out_fig 'Lat_Disch_Change_P05_' dir_out_name_fig_tag];
% %     h=figure;
% %     hold on
% %     func_shadedErrorBar(Lat',Variable,{@mean,@std},{'b','LineWidth',4});
% %     title({'Latitudinal average of normalized change in 5th percentile of Discharge' ; ['ISI-MIP - GHM: ' Impact_Models_Name ' , GCM: multiple - [2070-2099 RCP8.5] compared to [1971-2000]']})
% %     legend('GCM Average +/- St.Dev.', 'location', 'NorthWest');
% %     xlabel('Latitude ºN'); ylabel('Average of normalized change in discharge');
% %     set(gca,'FontSize',36, 'FontName', 'MS Sens Serif') % Axis Numbers and ranges Font
% %     set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
% %     set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
% %     xlim([-60,85])
% %     ylim([-1,1])
% %     plot([-90,90], [0,0],'LineWidth',0.5,'Color','black');
% %     set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
% %     set(gcf,'PaperPositionMode','auto') % The most important line for proper saving the figures ***
% %     saveas(h,dir_out_name_fig)
% %     print(dir_out_name_fig,'-dpng','-r300')
% %     close
% %     
% %     Variable=Multimodel_Q_Change_P95_hist_rcp8p5_norm_Lat(:, :, ghm_i);
% %     dir_out_name_fig=[dir_out_fig 'Lat_Disch_Change_P95_' dir_out_name_fig_tag];
% %     h=figure;
% %     hold on
% %     func_shadedErrorBar(Lat',Variable,{@mean,@std},{'r','LineWidth',4});
% %     title({'Latitudinal average of normalized change in 95th percentile of discharge' ; ['ISI-MIP - GHM: ' Impact_Models_Name ' , GCM: multiple - [2070-2099 RCP8.5] compared to [1971-2000]']})
% %     legend('GCM Average +/- St.Dev.', 'location', 'NorthWest');
% %     xlabel('Latitude ºN'); ylabel('Average of normalized change in discharge');
% %     set(gca,'FontSize',36, 'FontName', 'MS Sens Serif') % Axis Numbers and ranges Font
% %     set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
% %     set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
% %     xlim([-60,85])
% %     ylim([-1,1])
% %     plot([-90,90], [0,0],'LineWidth',0.5,'Color','black');
% %     set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
% %     set(gcf,'PaperPositionMode','auto') % The most important line for proper saving the figures ***
% %     saveas(h,dir_out_name_fig)
% %     print(dir_out_name_fig,'-dpng','-r300')
% %     close
% %     
% %     Variable1=Multimodel_Q_Change_P05_hist_rcp8p5_norm_Lat(:, :, ghm_i);
% %     Variable2=Multimodel_Q_Change_P95_hist_rcp8p5_norm_Lat(:, :, ghm_i);
% %     dir_out_name_fig=[dir_out_fig 'Lat_Disch_Change_P05_P95_' dir_out_name_fig_tag];
% %     h=figure;
% %     hold on
% %     func_shadedErrorBar(Lat',Variable1,{@mean,@std},{'b','LineWidth',4},1);
% %     func_shadedErrorBar(Lat',Variable2,{@mean,@std},{'r','LineWidth',4},1);
% %     title({'Latitudinal average of normalized change in 5th and 95th percentile of discharge' ; ['ISI-MIP - GHM: ' Impact_Models_Name ' , GCM: multiple - [2070-2099 RCP8.5] compared to [1971-2000]']})
% %     legend('GCM Average P05 +/- St.Dev.', 'GCM Average P95 +/- St.Dev.', 'location', 'NorthWest');
% %     xlabel('Latitude ºN'); ylabel('Average of normalized change in discharge');
% %     set(gca,'FontSize',36, 'FontName', 'MS Sens Serif') % Axis Numbers and ranges Font
% %     set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
% %     set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
% %     xlim([-60,85])
% %     ylim([-1,1])
% %     plot([-90,90], [0,0],'LineWidth',0.5,'Color','black');
% %     set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
% %     set(gcf,'PaperPositionMode','auto') % The most important line for proper saving the figures ***
% %     saveas(h,dir_out_name_fig)
% %     print(dir_out_name_fig,'-dpng','-r300')
% %     close
% %     
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     %%% Multi P05 vs P95 Scatter %%%
% %     help_P95_norm=reshape(nanmean(GHM_Disch_Change_P95_rcp8p5_hist_norm,3), [Lat_n*Lon_n,1]);
% %     help_P05_norm=reshape(nanmean(GHM_Disch_Change_P05_rcp8p5_hist_norm,3), [Lat_n*Lon_n,1]);
% %     GCM_P95_vs_P05_norm=help_P05_norm;
% %     GCM_P95_vs_P05_norm(:,2)=help_P95_norm;
% %     
% %     GCM_P95_vs_P05_norm=GCM_P95_vs_P05_norm(~isnan(GCM_P95_vs_P05_norm(:,2)),:); % Eliminating the rows with NaN values
% %     GCM_P95_vs_P05_norm=GCM_P95_vs_P05_norm(~isnan(GCM_P95_vs_P05_norm(:,1)),:);
% %     
% %     dir_out_name_fig=[dir_out_fig 'Scatter_Disch_Change_norm_P05vsP95_' dir_out_name_fig_tag];
% %     h1=scatter(GCM_P95_vs_P05_norm(:,1), GCM_P95_vs_P05_norm(:,2) * -1, 20 , 'fill');
% %     title({'Normalized change in 5th versus 95th percentile of discharge - Average of GCMs' ; ['ISI-MIP - GHM: ' Impact_Models_Name ' , GCM: multiple - [2070-2099 RCP8.5] compared to [1971-2000]']})
% %     xlabel({'Normalized change in 5th percentile of discharge '; '[2070-2099 RCP8.5] compared to [1971-2000]'});
% %     ylabel({'Normalized change in 95th percentile of discharge'; '[2070-2099 RCP8.5] compared to [1971-2000]'});
% %     xlim([-1,1])
% %     ylim([-1,1])
% %     set(gca,'XTickLabel',(1:-0.5:-1))
% %     set(gca,'FontSize',34, 'FontName', 'MS Sens Serif') % Axis Numbers and rages Font
% %     set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
% %     set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
% %     hold on
% %     plot([0,0], [-1,1],'LineWidth',1,'Color','black');
% %     hold on
% %     plot([-1,1], [0,0],'LineWidth',1,'Color','black');
% %     axis square
% %     set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
% %     set(gcf,'PaperPositionMode','auto') % The most important line for proper saving the figures ***
% %     saveas(h1,dir_out_name_fig)
% %     print(dir_out_name_fig,'-dpng','-r300')
% %     close
% %     
% %     
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     dir_out_name_fig=[dir_out_fig 'Map_Disch_Change_Median_' dir_out_name_fig_tag];
% %     Variable=nanmean(GHM_Disch_Change_Median_rcp8p5_hist, 3);
% %     %limits = max(abs(quantile(Variable(:), [0.05 0.95])));
% %     %var_limit=[-limits limits];
% %     var_limit=[-120 120];
% %     title_text={'Relative change in median of discharge (%) - Average of GCMs' ; ['ISI-MIP - GHM: ' Impact_Models_Name ' , GCM: multiple - [2070-2099 RCP8.5] compared to [1971-2000]']};
% %     [~]=func_GeoshowMap_Save(Variable, Lat_img, Lon_img, var_limit, x_limit, y_limit, dir_out_name_fig, title_text);
% %     
% %     dir_out_name_fig=[dir_out_fig 'Map_Disch_Change_P05_' dir_out_name_fig_tag];
% %     Variable=nanmean(GHM_Disch_Change_P05_rcp8p5_hist, 3);
% %     %limits = max(abs(quantile(Variable(:), [0.05 0.95])));
% %     %var_limit=[-limits limits];
% %     var_limit=[-120 120];
% %     title_text={'Relative change in 5th percentile of discharge (%) - Average of GCMs' ; ['ISI-MIP - GHM: ' Impact_Models_Name ' , GCM: multiple - [2070-2099 RCP8.5] compared to [1971-2000]']};
% %     [~]=func_GeoshowMap_Save(Variable, Lat_img, Lon_img, var_limit, x_limit, y_limit, dir_out_name_fig, title_text);
% %     
% %     dir_out_name_fig=[dir_out_fig 'Map_Disch_Change_P95_' dir_out_name_fig_tag];
% %     Variable=nanmean(GHM_Disch_Change_P95_rcp8p5_hist, 3);
% %     %limits = max(abs(quantile(Variable(:), [0.05 0.95])));
% %     %var_limit=[-limits limits];
% %     var_limit=[-120 120];
% %     title_text={'Relative change in 95th percentile of discharge (%) - Average of GCMs' ; ['ISI-MIP - GHM: ' Impact_Models_Name ' , GCM: multiple - [2070-2099 RCP8.5] compared to [1971-2000]']};
% %     [~]=func_GeoshowMap_Save(Variable, Lat_img, Lon_img, var_limit, x_limit, y_limit, dir_out_name_fig, title_text);
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     dir_out_name_fig=[dir_out_fig 'Map_norm_Disch_Change_Median_' dir_out_name_fig_tag];
% %     Variable=nanmean(GHM_Disch_Change_Median_rcp8p5_hist_norm, 3);
% %     var_limit=[-1 1];
% %     title_text={'Normalized change in median of discharge - Average of GCMs' ; ['ISI-MIP - GHM: ' Impact_Models_Name ' , GCM: multiple - [2070-2099 RCP8.5] compared to [1971-2000]']};
% %     [~]=func_GeoshowMap_Save(Variable, Lat_img, Lon_img, var_limit, x_limit, y_limit, dir_out_name_fig, title_text);
% %     
% %     dir_out_name_fig=[dir_out_fig 'Map_norm_Disch_Change_P05_' dir_out_name_fig_tag];
% %     Variable=nanmean(GHM_Disch_Change_P05_rcp8p5_hist_norm, 3);
% %     var_limit=[-1 1];
% %     title_text={'Normalized change in 5th percentile of discharge - Average of GCMs' ; ['ISI-MIP - GHM: ' Impact_Models_Name ' , GCM: multiple - [2070-2099 RCP8.5] compared to [1971-2000]']};
% %     [~]=func_GeoshowMap_Save(Variable, Lat_img, Lon_img, var_limit, x_limit, y_limit, dir_out_name_fig, title_text);
% %     
% %     dir_out_name_fig=[dir_out_fig 'Map_norm_Disch_Change_P95_' dir_out_name_fig_tag];
% %     Variable=nanmean(GHM_Disch_Change_P95_rcp8p5_hist_norm, 3);
% %     var_limit=[-1 1];
% %     title_text={'Normalized change in 95th percentile of discharge - Average of GCMs' ; ['ISI-MIP - GHM: ' Impact_Models_Name ' , GCM: multiple - [2070-2099 RCP8.5] compared to [1971-2000]']};
% %     [~]=func_GeoshowMap_Save(Variable, Lat_img, Lon_img, var_limit, x_limit, y_limit, dir_out_name_fig, title_text);
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
end

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
help_Pmed_hist = nanmean(Multimodel_Q_ave_Med_hist,3);
help_P05_hist  = nanmean(Multimodel_Q_ave_P05_hist,3);
help_P95_hist = nanmean(Multimodel_Q_ave_P95_hist,3);
help_Pmed_rcp = nanmean(Multimodel_Q_ave_Med_rcp8p5,3);
help_P05_rcp = nanmean(Multimodel_Q_ave_P05_rcp8p5,3);
help_P95_rcp = nanmean(Multimodel_Q_ave_P95_rcp8p5,3);

Multimodel_2_Q_Change_Med_norm=NaN(Lat_n, Lon_n);
Multimodel_2_Q_Change_P05_norm=NaN(Lat_n, Lon_n);
Multimodel_2_Q_Change_P95_norm=NaN(Lat_n, Lon_n);
for gcm_i=1:5
    
    for ii=1:Lat_n
        for jj=1:Lon_n
            if ~isnan(help_Pmed_rcp (ii,jj)) && ~isnan(help_Pmed_hist (ii,jj)) && help_Pmed_hist (ii,jj)~=0
                Multimodel_2_Q_Change_Med_norm(ii,jj) = (( help_Pmed_rcp(ii,jj) - help_Pmed_hist(ii,jj) ) / ( help_Pmed_rcp(ii,jj) + help_Pmed_hist(ii,jj) ) );
            end
            if ~isnan(help_P05_rcp(ii,jj)) && ~isnan(help_P05_hist(ii,jj)) && help_P05_hist(ii,jj)~=0
                Multimodel_2_Q_Change_P05_norm(ii,jj) = (( help_P05_rcp(ii,jj) - help_P05_hist (ii,jj) ) / ( help_P05_rcp(ii,jj) + help_P05_hist (ii,jj) ) );
            end
            if ~isnan(help_P95_rcp(ii,jj)) && ~isnan(help_P95_hist(ii,jj)) && help_P95_hist(ii,jj)~=0
                Multimodel_2_Q_Change_P95_norm(ii,jj) = (( help_P95_rcp(ii,jj) - help_P95_hist(ii,jj) ) / ( help_P95_rcp(ii,jj) + help_P95_hist(ii,jj) ) );
            end
            
        end
    end
    
end

[Multimodel_2_C_Ave_Q_Change_Med_norm, Multimodel_2_C_Ave_Q_Change_P05_norm, Multimodel_2_C_Ave_Q_Change_P95_norm]= ...
    func_TrendAve_CellAreaWeighted_3var_grid(Multimodel_2_Q_Change_Med_norm, Multimodel_2_Q_Change_P05_norm, Multimodel_2_Q_Change_P95_norm, Lat, Lon, Lat_bound, Lon_bound);

xlswrite(excel_out_name_1, Multimodel_2_C_Ave_Q_Change_Med_norm(1:7, 1)', 'Multimodel_Continents', 'D37:J37')
xlswrite(excel_out_name_1, Multimodel_2_C_Ave_Q_Change_P05_norm(1:7, 1)', 'Multimodel_Continents', 'L37:R37')
xlswrite(excel_out_name_1, Multimodel_2_C_Ave_Q_Change_P95_norm(1:7, 1)', 'Multimodel_Continents', 'T37:Z37')

%%% Normalized Change Averaged by Quadrants %%%
Multimodel_2_Quadrant_Ave_P95_P05_hist_rcp8p5_norm=NaN(1, 8); 
Var_P05=Multimodel_2_Q_Change_P05_norm;
Var_P95=Multimodel_2_Q_Change_P95_norm;
[Flood_Drought_Indicator_All]=func_FloodDroughtRisk_indicator(Var_P05, Var_P95);

Var_P05_help=Var_P05; Var_P95_help=Var_P95; GridCell_Areas_help=GridCell_Areas; % The Quadrant experiencing increased Flood and increased Drought (indicator 1)
Var_P05_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==2 | Flood_Drought_Indicator_All==3 | Flood_Drought_Indicator_All==4)=NaN; % excluding the cells that do not fit in this indicator
Var_P95_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==2 | Flood_Drought_Indicator_All==3 | Flood_Drought_Indicator_All==4)=NaN;
GridCell_Areas_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==2 | Flood_Drought_Indicator_All==3 | Flood_Drought_Indicator_All==4)=NaN;
Multimodel_2_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1, 1)=nansum(nansum( Var_P95_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help )); % Cell-Area-Weighted averaging
Multimodel_2_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1, 5)=nansum(nansum( Var_P05_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help ));
Var_P05_help=Var_P05; Var_P95_help=Var_P95; GridCell_Areas_help=GridCell_Areas; % The Quadrant experiencing increased Flood and decresed Drought (indicator 2)
Var_P05_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==1 | Flood_Drought_Indicator_All==3 | Flood_Drought_Indicator_All==4)=NaN; % excluding the cells that do not fit in this indicator
Var_P95_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==1 | Flood_Drought_Indicator_All==3 | Flood_Drought_Indicator_All==4)=NaN;
GridCell_Areas_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==1 | Flood_Drought_Indicator_All==3 | Flood_Drought_Indicator_All==4)=NaN;
Multimodel_2_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1, 2)=nansum(nansum( Var_P95_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help )); % Cell-Area-Weighted averaging
Multimodel_2_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1, 6)=nansum(nansum( Var_P05_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help ));
Var_P05_help=Var_P05; Var_P95_help=Var_P95; GridCell_Areas_help=GridCell_Areas; % The Quadrant experiencing decresed Flood and increased Drought (indicator 3)
Var_P05_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==1 | Flood_Drought_Indicator_All==2 | Flood_Drought_Indicator_All==4)=NaN; % excluding the cells that do not fit in this indicator
Var_P95_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==1 | Flood_Drought_Indicator_All==2 | Flood_Drought_Indicator_All==4)=NaN;
GridCell_Areas_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==1 | Flood_Drought_Indicator_All==2 | Flood_Drought_Indicator_All==4)=NaN;
Multimodel_2_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1, 3)=nansum(nansum( Var_P95_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help )); % Cell-Area-Weighted averaging
Multimodel_2_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1, 7)=nansum(nansum( Var_P05_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help ));
Var_P05_help=Var_P05; Var_P95_help=Var_P95; GridCell_Areas_help=GridCell_Areas; % The Quadrant experiencing decresed Flood and decresed Drought (indicator 4)
Var_P05_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==1 | Flood_Drought_Indicator_All==2 | Flood_Drought_Indicator_All==3)=NaN; % excluding the cells that do not fit in this indicator
Var_P95_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==1 | Flood_Drought_Indicator_All==2 | Flood_Drought_Indicator_All==3)=NaN;
GridCell_Areas_help(isnan(Flood_Drought_Indicator_All) | Flood_Drought_Indicator_All==1 | Flood_Drought_Indicator_All==2 | Flood_Drought_Indicator_All==3)=NaN;
Multimodel_2_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1, 4)=nansum(nansum( Var_P95_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help )); % Cell-Area-Weighted averaging
Multimodel_2_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1, 8)=nansum(nansum( Var_P05_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help ));

%%% Writing the results in excell file %%%
xlswrite(excel_out_name_1, Multimodel_2_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1,1), 'Multimodel_Quadrants', 'E37')
xlswrite(excel_out_name_1, Multimodel_2_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1,2), 'Multimodel_Quadrants', 'D37')
xlswrite(excel_out_name_1, Multimodel_2_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1,3), 'Multimodel_Quadrants', 'F37')
xlswrite(excel_out_name_1, Multimodel_2_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1,4), 'Multimodel_Quadrants', 'G37')
xlswrite(excel_out_name_1, Multimodel_2_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1,5) * -1, 'Multimodel_Quadrants', 'J37')
xlswrite(excel_out_name_1, Multimodel_2_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1,6) * -1, 'Multimodel_Quadrants', 'I37')
xlswrite(excel_out_name_1, Multimodel_2_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1,7) * -1, 'Multimodel_Quadrants', 'K37')
xlswrite(excel_out_name_1, Multimodel_2_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1,8) * -1, 'Multimodel_Quadrants', 'L37')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Multi GHM/GCM Average Maps %%%
dir_out_fig=[pwd '\Figures - Multimodel - RCP8.5\' ]; % Directory to save Figures and Maps
dir_out_name_fig_tag='Disch_ISIMIP_GHM_GCM_averaged_rcp8p5_hist';
x_limit=[-180 180]; y_limit=[-60 90];
Lon_img=Lon'; Lat_img=Lat;

dir_out_name_fig=[dir_out_fig 'Map_norm_Disch_Change_Median_' dir_out_name_fig_tag];
Variable=nanmean(Multimodel_Q_Change_Med_hist_rcp8p5_norm, 3);
var_limit=[-0.6 0.6];
title_text={'Normalized change in median of discharge - Multimodel average of GHMs and GCMs' ; ('ISI-MIP - GHM: multiple , GCM: multiple - [2070-2099 RCP8.5] compared to [1971-2000]')};
[~]=func_GeoshowMap_Save_RedBlue(Variable, Lat_img, Lon_img, var_limit, x_limit, y_limit, dir_out_name_fig, title_text);

dir_out_name_fig=[dir_out_fig 'Map_norm_Disch_Change_P05_' dir_out_name_fig_tag];
Variable=nanmean(Multimodel_Q_Change_P05_hist_rcp8p5_norm, 3);
var_limit=[-0.6 0.6];
title_text={'Normalized change in 5th percentile of discharge - Multimodel average of GHMs and GCMs' ; ('ISI-MIP - GHM: multiple , GCM: multiple - [2070-2099 RCP8.5] compared to [1971-2000]')};
[~]=func_GeoshowMap_Save_RedBlue(Variable, Lat_img, Lon_img, var_limit, x_limit, y_limit, dir_out_name_fig, title_text);

dir_out_name_fig=[dir_out_fig 'Map_norm_Disch_Change_P95_' dir_out_name_fig_tag];
Variable=nanmean(Multimodel_Q_Change_P95_hist_rcp8p5_norm, 3);
var_limit=[-0.6 0.6];
title_text={'Normalized change in 95th percentile of discharge - Multimodel average of GHMs and GCMs' ; ('ISI-MIP - GHM: multiple , GCM: multiple - [2070-2099 RCP8.5] compared to [1971-2000]')};
[~]=func_GeoshowMap_Save_RedBlue(Variable, Lat_img, Lon_img, var_limit, x_limit, y_limit, dir_out_name_fig, title_text);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Flood/Drought Risk Spots %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Model Agreement on Trend Significance %%%

Multimodel_No_Data_P05 = nansum( ~isnan(Multimodel_Q_Change_P05_hist_rcp8p5_norm) , 3);
Multimodel_No_Data_P05(Multimodel_No_Data_P05 < 20 ) =NaN;
Multimodel_No_Data_P95 = nansum( ~isnan(Multimodel_Q_Change_P95_hist_rcp8p5_norm) , 3);
Multimodel_No_Data_P95(Multimodel_No_Data_P95 < 20 ) =NaN;

dir_out_fig=[pwd '\Figures - Multimodel - RCP8.5\' ]; % Directory to save Figures and Maps
dir_out_name_fig_tag='Disch_ISIMIP_GHM_GCM_averaged_rcp8p5_hist';
x_limit=[-180 180]; y_limit=[-60 90];
Lon_img=Lon'; Lat_img=Lat;

dir_out_name_fig=[dir_out_fig 'Map_Flood_Drought_Risk_' dir_out_name_fig_tag];
Var_P05=nanmean(Multimodel_Q_Change_P05_hist_rcp8p5_norm,3);
Var_P95=nanmean(Multimodel_Q_Change_P95_hist_rcp8p5_norm,3);
title_text={'Change in flood and drought risk - Multimodel average of GHMs and GCMs' ; ('ISI-MIP - [2070-2099 RCP8.5] compared to [1971-2000]')};
[~, Multimodel_Flood_Drought_Indicator_All]=func_FloodDroughtRisk_Map_40c(Var_P05, Var_P95, Lat_img, Lon_img, x_limit, y_limit, dir_out_name_fig, title_text);

help_sig_P05=sum(Multimodel_T_test_Significance_Q_P05_hist_rcp8p5,3); % This sums all the postive and negative significances - But a negative sig WILL decrease the sum pf postives or vise versa
help_sig_P95=sum(Multimodel_T_test_Significance_Q_P95_hist_rcp8p5,3);
Multimodel_T_test_Significance_Q_P05_Sum=zeros(Lat_n, Lon_n);  % This sums all the postive and negative significances - But a negative sig WILL NOT decrease the sum pf postives or vise versa
Multimodel_T_test_Significance_Q_P05_Sum(isnan(help_sig_P05))=NaN;
Multimodel_T_test_Significance_Q_P95_Sum=zeros(Lat_n, Lon_n);
Multimodel_T_test_Significance_Q_P95_Sum(isnan(help_sig_P95))=NaN;

for t=1:size(Multimodel_T_test_Significance_Q_P05_hist_rcp8p5,3)
    
    for ii=1:Lat_n
        for jj=1:Lon_n
            
            if  (help_sig_P05(ii,jj) > 0) && (Multimodel_T_test_Significance_Q_P05_hist_rcp8p5(ii,jj,t)==1) % the help_sig_P05 is used to locate the postitive-sig cells and only increase the number of agreement between models in postive sig in only those cells
                Multimodel_T_test_Significance_Q_P05_Sum(ii,jj)=Multimodel_T_test_Significance_Q_P05_Sum(ii,jj)+1;
            elseif (help_sig_P05(ii,jj) < 0) && (Multimodel_T_test_Significance_Q_P05_hist_rcp8p5(ii,jj,t)==-1) % the help_sig_P05 is used to locate the negative-sig cells and only increase the number of agreement between models in negative sig in only those cells
                Multimodel_T_test_Significance_Q_P05_Sum(ii,jj)=Multimodel_T_test_Significance_Q_P05_Sum(ii,jj)-1;
            end
            
            if  (help_sig_P95(ii,jj) > 0) && (Multimodel_T_test_Significance_Q_P95_hist_rcp8p5(ii,jj,t)==1)
                Multimodel_T_test_Significance_Q_P95_Sum(ii,jj)=Multimodel_T_test_Significance_Q_P95_Sum(ii,jj)+1;
            elseif (help_sig_P95(ii,jj) < 0) && (Multimodel_T_test_Significance_Q_P95_hist_rcp8p5(ii,jj,t)==-1)
                Multimodel_T_test_Significance_Q_P95_Sum(ii,jj)=Multimodel_T_test_Significance_Q_P95_Sum(ii,jj)-1;
            end
            
        end
    end
    
end
Multimodel_T_test_Significance_Q_P05_Sum_10=Multimodel_T_test_Significance_Q_P05_Sum;
Multimodel_T_test_Significance_Q_P95_Sum_10=Multimodel_T_test_Significance_Q_P95_Sum;
Multimodel_T_test_Significance_Q_P05_Sum_10(Multimodel_T_test_Significance_Q_P05_Sum_10 < 10 & Multimodel_T_test_Significance_Q_P05_Sum_10 > -10)=NaN;% Cells with less than 10 out of 20 models agreeing on the significance are excluded
Multimodel_T_test_Significance_Q_P95_Sum_10(Multimodel_T_test_Significance_Q_P95_Sum_10 < 10 & Multimodel_T_test_Significance_Q_P95_Sum_10 > -10)=NaN;% Cells with less than 10 out of 20 models agreeing on the significance are excluded

dir_out_name_fig=[dir_out_fig 'Map_Flood_Drought_Risk_Sig10_' dir_out_name_fig_tag];
Var_P05=nanmean(Multimodel_Q_Change_P05_hist_rcp8p5_norm,3);
Var_P05(isnan(Multimodel_T_test_Significance_Q_P05_Sum_10))=NaN;
Var_P95=nanmean(Multimodel_Q_Change_P95_hist_rcp8p5_norm,3);
Var_P95(isnan(Multimodel_T_test_Significance_Q_P95_Sum_10))=NaN;
title_text={'Change in flood and drought risk - Multimodel average of GHMs and GCMs' ; ['Gird cells with changes significant at 95% confidence level in at least 10 out of ' num2str(n_mods) ' models']; 'ISI-MIP - [2070-2099 RCP8.5] compared to [1971-2000]'};
[~, ~]=func_FloodDroughtRisk_Map_4c(Var_P05, Var_P95, Lat_img, Lon_img, x_limit, y_limit, dir_out_name_fig, title_text);

Multimodel_T_test_Significance_Q_P05_Sum_half=Multimodel_T_test_Significance_Q_P05_Sum;
Multimodel_T_test_Significance_Q_P95_Sum_half=Multimodel_T_test_Significance_Q_P95_Sum;
Multimodel_T_test_Significance_Q_P05_Sum_half(Multimodel_T_test_Significance_Q_P05_Sum_half < (ceil(Multimodel_No_Data_P05/2)) & Multimodel_T_test_Significance_Q_P05_Sum_half > (-1 * ceil(Multimodel_No_Data_P05/2)))=NaN;% Cells with less than half of the models agreeing on the significance are excluded
Multimodel_T_test_Significance_Q_P95_Sum_half(Multimodel_T_test_Significance_Q_P95_Sum_half < (ceil(Multimodel_No_Data_P95/2)) & Multimodel_T_test_Significance_Q_P95_Sum_half > (-1 * ceil(Multimodel_No_Data_P95/2)))=NaN;% Cells with less than half of the models agreeing on the significance are excluded

dir_out_name_fig=[dir_out_fig 'Map_Flood_Drought_Risk_Sig_half_' dir_out_name_fig_tag];
Var_P05=nanmean(Multimodel_Q_Change_P05_hist_rcp8p5_norm,3);
Var_P05(isnan(Multimodel_T_test_Significance_Q_P05_Sum_half))=NaN;
Var_P95=nanmean(Multimodel_Q_Change_P95_hist_rcp8p5_norm,3);
Var_P95(isnan(Multimodel_T_test_Significance_Q_P95_Sum_half))=NaN;
title_text={'Change in flood and drought risk - Multimodel average of GHMs and GCMs' ; 'Gird cells with changes significant at 95% confidence level in at least half of the models'; 'ISI-MIP - [2070-2099 RCP8.5] compared to [1971-2000]'};
[~, Multimodel_Flood_Drought_Indicator_Sig_half]=func_FloodDroughtRisk_Map_4c(Var_P05, Var_P95, Lat_img, Lon_img, x_limit, y_limit, dir_out_name_fig, title_text);

Variable1 = ( reshape(nanmean( abs(Multimodel_T_test_Significance_Q_P05_Sum) ,2), [Lat_n, 1]) )';
Variable2 = ( reshape(nanmean( abs(Multimodel_T_test_Significance_Q_P95_Sum) ,2), [Lat_n, 1]) )';
dir_out_name_fig=[dir_out_fig 'Lat_Sig_Model_No_P05_P95_' dir_out_name_fig_tag];
h=figure;
hold on
plot(Lat', Variable1,'LineWidth',4,'Color','red');
plot(Lat', Variable2,'LineWidth',4,'Color','blue');
title({['Latitudinal average number of models (out of ' num2str(n_mods) ') agreeing on the significance of change'] ; ('ISI-MIP - Change in 5th and 95th percentile of discharge - [2070-2099 RCP8.5] vs. [1971-2000]')})
legend('P05 - ave. no. of models agreeing', 'P95 - ave. no. of models agreeing', 'location', 'SouthEast');
xlabel('Latitude ºN'); ylabel(['Number of models (out of ' num2str(n_mods) ')']);
set(gca,'FontSize',36, 'FontName', 'MS Sens Serif') % Axis Numbers and ranges Font
set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
xlim([-60,85])
ylim([0,n_mods])
plot([-90,90], [5,5],':','LineWidth',0.5,'Color','black');
plot([-90,90], [10,10],':','LineWidth',0.5,'Color','black');
plot([-90,90], [15,15],':','LineWidth',0.5,'Color','black');
plot([-90,90], [20,20],':','LineWidth',0.5,'Color','black');
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
set(gcf,'PaperPositionMode','auto') % The most important line for proper saving the figures ***
saveas(h,dir_out_name_fig)
print(dir_out_name_fig,'-dpng','-r300')
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Model Agreement on Trend Direction %%%
Multimodel_Q_Trend_Dir_P05=sign(Multimodel_Q_Change_P05_hist_rcp8p5_norm); % Direction of Trend in the gridcell
Multimodel_Q_Trend_Dir_P95=sign(Multimodel_Q_Change_P95_hist_rcp8p5_norm);

help_dir_P05=sum(Multimodel_Q_Trend_Dir_P05,3); % This sums all the postive and negative significances - But a negative sig WILL decrease the sum pf postives or vise versa
help_dir_P95=sum(Multimodel_Q_Trend_Dir_P95,3);
Multimodel_Q_Trend_Dir_P05_Sum=zeros(Lat_n, Lon_n);  % This sums all the postive and negative significances - But a negative sig WILL NOT decrease the sum pf postives or vise versa
Multimodel_Q_Trend_Dir_P05_Sum(isnan(help_dir_P05))=NaN;
Multimodel_Q_Trend_Dir_P95_Sum=zeros(Lat_n, Lon_n);
Multimodel_Q_Trend_Dir_P95_Sum(isnan(help_dir_P95))=NaN;

for t=1:size(Multimodel_Q_Trend_Dir_P05,3)
    
    for ii=1:Lat_n
        for jj=1:Lon_n
            
            if  (help_dir_P05(ii,jj) > 0) && (Multimodel_Q_Trend_Dir_P05(ii,jj,t)==1) % the help_dir_P05 is used to locate the postitive-sig cells and only increase the number of agreement between models in postive sig in only those cells
                Multimodel_Q_Trend_Dir_P05_Sum(ii,jj)=Multimodel_Q_Trend_Dir_P05_Sum(ii,jj)+1;
            elseif (help_dir_P05(ii,jj) < 0) && (Multimodel_Q_Trend_Dir_P05(ii,jj,t)==-1) % the help_dir_P05 is used to locate the negative-sig cells and only increase the number of agreement between models in negative sig in only those cells
                Multimodel_Q_Trend_Dir_P05_Sum(ii,jj)=Multimodel_Q_Trend_Dir_P05_Sum(ii,jj)-1;
            end
            
            if  (help_dir_P95(ii,jj) > 0) && (Multimodel_Q_Trend_Dir_P95(ii,jj,t)==1)
                Multimodel_Q_Trend_Dir_P95_Sum(ii,jj)=Multimodel_Q_Trend_Dir_P95_Sum(ii,jj)+1;
            elseif (help_dir_P95(ii,jj) < 0) && (Multimodel_Q_Trend_Dir_P95(ii,jj,t)==-1)
                Multimodel_Q_Trend_Dir_P95_Sum(ii,jj)=Multimodel_Q_Trend_Dir_P95_Sum(ii,jj)-1;
            end
            
        end
    end
    
end
Multimodel_Q_Trend_Dir_P05_Sum_half=Multimodel_Q_Trend_Dir_P05_Sum;
Multimodel_Q_Trend_Dir_P95_Sum_half=Multimodel_Q_Trend_Dir_P95_Sum;
Multimodel_Q_Trend_Dir_P05_Sum_half(Multimodel_Q_Trend_Dir_P05_Sum_half < (ceil(Multimodel_No_Data_P05/2)) & Multimodel_Q_Trend_Dir_P05_Sum_half > (-1 * ceil(Multimodel_No_Data_P05/2)))=NaN; % Cells with less than half of the models agreeing on the trend direction are excluded
Multimodel_Q_Trend_Dir_P95_Sum_half(Multimodel_Q_Trend_Dir_P95_Sum_half < (ceil(Multimodel_No_Data_P05/2)) & Multimodel_Q_Trend_Dir_P95_Sum_half > (-1 * ceil(Multimodel_No_Data_P05/2)))=NaN; % Cells with less than half of the models agreeing on the trend direction are excluded

dir_out_name_fig=[dir_out_fig 'Map_Flood_Drought_Risk_Dir_half_' dir_out_name_fig_tag];
Var_P05=nanmean(Multimodel_Q_Change_P05_hist_rcp8p5_norm,3);
Var_P05(isnan(Multimodel_Q_Trend_Dir_P05_Sum_half))=NaN;
Var_P95=nanmean(Multimodel_Q_Change_P95_hist_rcp8p5_norm,3);
Var_P95(isnan(Multimodel_Q_Trend_Dir_P95_Sum_half))=NaN;
title_text={'Change in flood and drought risk - Multimodel average of GHMs and GCMs' ; 'Gird cells with agreement on direction of trend in at least half of the models'; 'ISI-MIP - [2070-2099 RCP8.5] compared to [1971-2000]'};
[~, Multimodel_Flood_Drought_Indicator_Dir_half]=func_FloodDroughtRisk_Map_4c(Var_P05, Var_P95, Lat_img, Lon_img, x_limit, y_limit, dir_out_name_fig, title_text);

Multimodel_Q_Trend_Dir_P05_Sum_twothirds=Multimodel_Q_Trend_Dir_P05_Sum;
Multimodel_Q_Trend_Dir_P95_Sum_twothirds=Multimodel_Q_Trend_Dir_P95_Sum;
Multimodel_Q_Trend_Dir_P05_Sum_twothirds(Multimodel_Q_Trend_Dir_P05_Sum_twothirds < (ceil(Multimodel_No_Data_P05*(2/3))) & Multimodel_Q_Trend_Dir_P05_Sum_twothirds > -1 * (ceil(Multimodel_No_Data_P05*(2/3))))=NaN; % Cells with less than two-thirds of the models agreeing on the trend direction are excluded
Multimodel_Q_Trend_Dir_P95_Sum_twothirds(Multimodel_Q_Trend_Dir_P95_Sum_twothirds < (ceil(Multimodel_No_Data_P05*(2/3))) & Multimodel_Q_Trend_Dir_P95_Sum_twothirds > -1 * (ceil(Multimodel_No_Data_P05*(2/3))))=NaN; % Cells with less than two-thirds of the models agreeing on the trend direction are excluded

dir_out_name_fig=[dir_out_fig 'Map_Flood_Drought_Risk_Dir_twothirds_' dir_out_name_fig_tag];
Var_P05=nanmean(Multimodel_Q_Change_P05_hist_rcp8p5_norm,3);
Var_P05(isnan(Multimodel_Q_Trend_Dir_P05_Sum_twothirds))=NaN;
Var_P95=nanmean(Multimodel_Q_Change_P95_hist_rcp8p5_norm,3);
Var_P95(isnan(Multimodel_Q_Trend_Dir_P95_Sum_twothirds))=NaN;
title_text={'Change in flood and drought risk - Multimodel average of GHMs and GCMs' ; 'Gird cells with agreement on direction of trend in at least two-thirds of the models'; 'ISI-MIP - [2070-2099 RCP8.5] compared to [1971-2000]'};
[~, Multimodel_Flood_Drought_Indicator_Dir_twothirds]=func_FloodDroughtRisk_Map_4c(Var_P05, Var_P95, Lat_img, Lon_img, x_limit, y_limit, dir_out_name_fig, title_text);

Variable1 = ( reshape(nanmean( abs(Multimodel_Q_Trend_Dir_P05_Sum) ,2), [Lat_n, 1]) )';
Variable2 = ( reshape(nanmean( abs(Multimodel_Q_Trend_Dir_P95_Sum) ,2), [Lat_n, 1]) )';
dir_out_name_fig=[dir_out_fig 'Lat_Dir_Model_No_P05_P95_' dir_out_name_fig_tag];
h=figure;
hold on
plot(Lat', Variable1,'LineWidth',4,'Color','red');
plot(Lat', Variable2,'LineWidth',4,'Color','blue');
title({['Latitudinal average number of models (out of ' num2str(n_mods) ') agreeing on the direction of change'] ; ('ISI-MIP - Change in 5th and 95th percentile of discharge - [2070-2099 RCP8.5] vs. [1971-2000]')})
legend('P05 - ave. no. of models agreeing', 'P95 - ave. no. of models agreeing', 'location', 'SouthEast');
xlabel('Latitude ºN'); ylabel(['Number of models (out of ' num2str(n_mods) ')']);
set(gca,'FontSize',36, 'FontName', 'MS Sens Serif') % Axis Numbers and ranges Font
set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
xlim([-60,85])
ylim([0,n_mods])
plot([-90,90], [5,5],':','LineWidth',0.5,'Color','black');
plot([-90,90], [10,10],':','LineWidth',0.5,'Color','black');
plot([-90,90], [15,15],':','LineWidth',0.5,'Color','black');
plot([-90,90], [20,20],':','LineWidth',0.5,'Color','black');
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
set(gcf,'PaperPositionMode','auto') % The most important line for proper saving the figures ***
saveas(h,dir_out_name_fig)
print(dir_out_name_fig,'-dpng','-r300')
close

%%% Multi P05 vs P95 Scatter %%%
help_multi_P95_norm=reshape(nanmean(Multimodel_Q_Change_P05_hist_rcp8p5_norm,3), [Lat_n*Lon_n,1]);
help_multi_P05_norm=reshape(nanmean(Multimodel_Q_Change_P95_hist_rcp8p5_norm,3), [Lat_n*Lon_n,1]);
GCM_GHM_P95_vs_P05_norm=help_multi_P05_norm;
GCM_GHM_P95_vs_P05_norm(:,2)=help_multi_P95_norm;

GCM_GHM_P95_vs_P05_norm=GCM_GHM_P95_vs_P05_norm(~isnan(GCM_GHM_P95_vs_P05_norm(:,2)),:); % Eliminating the rows with NaN values
GCM_GHM_P95_vs_P05_norm=GCM_GHM_P95_vs_P05_norm(~isnan(GCM_GHM_P95_vs_P05_norm(:,1)),:);

dir_out_name_fig=[dir_out_fig 'Scatter_Disch_Change_norm_P05vsP95_' dir_out_name_fig_tag];
h1=scatter(GCM_GHM_P95_vs_P05_norm(:,1), GCM_GHM_P95_vs_P05_norm(:,2) * -1, 40 , 'x'); %, 'MarkerEdgeColor', [0.9 0.55 0.05], 'MarkerFaceColor',[0.9 0.55 0.05]);
title({'Normalized change in 5th versus 95th percentile of discharge - Multimodel average of GHMs and GCMs' ; ('ISI-MIP - GHM: multiple , GCM: multiple - [2070-2099 RCP8.5] compared to [1971-2000]')})
xlabel({'Normalized change in drought risk (P05 *-1)'});
ylabel({'Normalized change in flood risk (P95)'});
xlim([-1,1])
ylim([-1,1])
set(gca,'XTickLabel',(-1:0.5:1))
set(gca,'FontSize',34, 'FontName', 'MS Sens Serif') % Axis Numbers and rages Font
set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
hold on
plot([0,0], [-1,1],'LineWidth',1,'Color','black');
hold on
plot([-1,1], [0,0],'LineWidth',1,'Color','black');
axis square
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
set(gcf,'PaperPositionMode','auto') % The most important line for proper saving the figures ***
saveas(h1,dir_out_name_fig)
print(dir_out_name_fig,'-dpng','-r300')
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % for gcm_i=1:size(GCM_Models_Names,2)
% %     
% %     Model_Name=GCM_Models_Names{1, gcm_i};
% %     
% %     dir_out_fig=[pwd '\Figures - Global MultiGHM - RCP8.5\' Model_Name '\']; % Directory to save Figures and Maps
% %     dir_out_name_fig_tag=['Disch_ISIMIP_' Model_Name '_MultiGHM_rcp8p5_hist'];
% %     
% %     %%% Average by Latitude %%%
% %     Variable= (reshape (Multimodel_Q_Change_Med_hist_rcp8p5_norm_Lat(gcm_i, :, :), [Lat_n, numel(Impact_Models_Names)]))' ;
% %     dir_out_name_fig=[dir_out_fig 'Lat_Disch_Change_Median_' dir_out_name_fig_tag];
% %     h=figure;
% %     hold on
% %     func_shadedErrorBar(Lat',Variable,{@mean,@std},{'g','LineWidth',4});
% %     title({'Latitudinal average of normalized change in median of discharge' ; ['ISI-MIP - GHM: multiple (Average +/- St.Dev.) , GCM: ' Model_Name ' - [2070-2099 RCP8.5] vs. [1971-2000]']})
% %     legend('GHM Average +/- St.Dev.','GHM Average', 'location', 'NorthWest');
% %     xlabel('Latitude ºN'); ylabel('Average of normalized change in discharge');
% %     set(gca,'FontSize',36, 'FontName', 'MS Sens Serif') % Axis Numbers and ranges Font
% %     set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
% %     set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
% %     xlim([-60,85])
% %     ylim([-1,1])
% %     plot([-90,90], [0,0],'LineWidth',0.5,'Color','black');
% %     set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
% %     set(gcf,'PaperPositionMode','auto') % The most important line for proper saving the figures ***
% %     saveas(h,dir_out_name_fig)
% %     print(dir_out_name_fig,'-dpng','-r300')
% %     close
% %     
% % end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Multi GHM/GCM Average by Latitude %%%
%%% GHM: Single , GCM: multiple (Average +/- St.Dev.) %%%
dir_out_fig=[pwd '\Figures - Multimodel - RCP8.5\' ]; % Directory to save Figures and Maps
dir_out_name_fig_tag='ISIMIP_GHMsingle_GCMaveraged_rcp8p5_hist';

dir_out_name_fig=[dir_out_fig 'Lat_Disch_Change_Median_' dir_out_name_fig_tag];
h=figure;
hold on
colors={'b','r','g','m','c','y','k'};
for ghm_i=1:size(Impact_Models_Names,2)
    Variable= Multimodel_Q_Change_Med_hist_rcp8p5_norm_Lat(:, :, ghm_i);
    func_shadedErrorBar(Lat',Variable,{@mean,@std},{colors{ghm_i},'LineWidth',4},1);
end
title({'Latitudinal average of normalized change in median of discharge' ; ('ISI-MIP - GHM: single , GCM: multiple (Average +/- St.Dev.) - [2070-2099 RCP8.5] vs. [1971-2000]')})
legend([Impact_Models_Names{1}], [Impact_Models_Names{2}], [Impact_Models_Names{3}], [Impact_Models_Names{4}], [Impact_Models_Names{5}], 'location', 'NorthWest');
xlabel('Latitude ºN'); ylabel('Average of normalized change in discharge');
set(gca,'FontSize',30, 'FontName', 'MS Sens Serif') % Axis Numbers and ranges Font
set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
xlim([-60,85])
ylim([-1,1])
plot([-90,90], [0,0],'LineWidth',0.5,'Color','black');
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
set(gcf,'PaperPositionMode','auto') % The most important line for proper saving the figures ***
saveas(h,dir_out_name_fig)
print(dir_out_name_fig,'-dpng','-r300')
close

dir_out_name_fig=[dir_out_fig 'Lat_Disch_Change_P05_' dir_out_name_fig_tag];
h=figure;
hold on
colors={'b','r','g','m','c','y','k'};
for ghm_i=1:size(Impact_Models_Names,2)
    Variable= Multimodel_Q_Change_P05_hist_rcp8p5_norm_Lat(:, :, ghm_i) * -1;
    func_shadedErrorBar(Lat',Variable,{@mean,@std},{colors{ghm_i},'LineWidth',4},1);
end
title({'Latitudinal average of normalized change in drought risk (-1 * 5th percentile of discharge)' ; ('ISI-MIP - GHM: single , GCM: multiple (Average +/- St.Dev.) - [2070-2099 RCP8.5] vs. [1971-2000]')})
legend([Impact_Models_Names{1}], [Impact_Models_Names{2}], [Impact_Models_Names{3}], [Impact_Models_Names{4}], [Impact_Models_Names{5}], 'location', 'SouthWest');
xlabel('Latitude ºN'); ylabel('Average of normalized change in drought');
set(gca,'FontSize',30, 'FontName', 'MS Sens Serif') % Axis Numbers and ranges Font
set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
xlim([-60,85])
ylim([-1,1])
plot([-90,90], [0,0],'LineWidth',0.5,'Color','black');
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
set(gcf,'PaperPositionMode','auto') % The most important line for proper saving the figures ***
saveas(h,dir_out_name_fig)
print(dir_out_name_fig,'-dpng','-r300')
close

dir_out_name_fig=[dir_out_fig 'Lat_Disch_Change_P95_' dir_out_name_fig_tag];
h=figure;
hold on
colors={'b','r','g','m','c','y','k'};
for ghm_i=1:size(Impact_Models_Names,2)
    Variable= Multimodel_Q_Change_P95_hist_rcp8p5_norm_Lat(:, :, ghm_i);
    func_shadedErrorBar(Lat',Variable,{@mean,@std},{colors{ghm_i},'LineWidth',4},1);
end
title({'Latitudinal average of normalized change in flood risk (95th percentile of discharge)' ; ('ISI-MIP - GHM: single , GCM: multiple (Average +/- St.Dev.) - [2070-2099 RCP8.5] vs. [1971-2000]')})
legend([Impact_Models_Names{1}], [Impact_Models_Names{2}], [Impact_Models_Names{3}], [Impact_Models_Names{4}], [Impact_Models_Names{5}], 'location', 'NorthWest');
xlabel('Latitude ºN'); ylabel('Average of normalized change in flood');
set(gca,'FontSize',30, 'FontName', 'MS Sens Serif') % Axis Numbers and ranges Font
set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
xlim([-60,85])
ylim([-1,1])
plot([-90,90], [0,0],'LineWidth',0.5,'Color','black');
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
set(gcf,'PaperPositionMode','auto') % The most important line for proper saving the figures ***
saveas(h,dir_out_name_fig)
print(dir_out_name_fig,'-dpng','-r300')
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Multi GHM/GCM Average by Latitude %%%
%%% GHM: multiple (Average +/- St.Dev.) , GCM: Single %%%
dir_out_fig=[pwd '\Figures - Multimodel - RCP8.5\' ]; % Directory to save Figures and Maps
dir_out_name_fig_tag='ISIMIP_GHMaveraged_GCMsingle_rcp8p5_hist';

dir_out_name_fig=[dir_out_fig 'Lat_Disch_Change_Median_' dir_out_name_fig_tag];
h=figure;
hold on
colors={'b','r','g','m','c','y','k'};
for gcm_i=1:size(GCM_Models_Names,2)
    Variable= (reshape (Multimodel_Q_Change_Med_hist_rcp8p5_norm_Lat(gcm_i, :, :), [Lat_n, numel(Impact_Models_Names)]))' ;
    func_shadedErrorBar(Lat',Variable,{@mean,@std},{colors{gcm_i},'LineWidth',4},1);
end
title({'Latitudinal average of normalized change in median of discharge' ; ('ISI-MIP - GHM: multiple (Average +/- St.Dev.) , GCM: single - [2070-2099 RCP8.5] vs. [1971-2000]')})
legend([GCM_Models_Names{1}], [GCM_Models_Names{2}], [GCM_Models_Names{3}], [GCM_Models_Names{4}], [GCM_Models_Names{5}], 'location', 'NorthWest');
xlabel('Latitude ºN'); ylabel('Average of normalized change in discharge');
set(gca,'FontSize',30, 'FontName', 'MS Sens Serif') % Axis Numbers and ranges Font
set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
xlim([-60,85])
ylim([-1,1])
plot([-90,90], [0,0],'LineWidth',0.5,'Color','black');
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
set(gcf,'PaperPositionMode','auto') % The most important line for proper saving the figures ***
saveas(h,dir_out_name_fig)
print(dir_out_name_fig,'-dpng','-r300')
close

dir_out_name_fig=[dir_out_fig 'Lat_Disch_Change_P05_' dir_out_name_fig_tag];
h=figure;
hold on
colors={'b','r','g','m','c','y','k'};
for gcm_i=1:size(GCM_Models_Names,2)
    Variable= (reshape (Multimodel_Q_Change_P05_hist_rcp8p5_norm_Lat(gcm_i, :, :), [Lat_n, numel(Impact_Models_Names)]))' * -1 ;
    func_shadedErrorBar(Lat',Variable,{@mean,@std},{colors{gcm_i},'LineWidth',4},1);
end
title({'Latitudinal average of normalized change in drought risk (-1 * 5th percentile of discharge)' ; ('ISI-MIP - GHM: multiple (Average +/- St.Dev.) , GCM: single - [2070-2099 RCP8.5] vs. [1971-2000]')})
legend([GCM_Models_Names{1}], [GCM_Models_Names{2}], [GCM_Models_Names{3}], [GCM_Models_Names{4}], [GCM_Models_Names{5}], 'location', 'SouthWest');
xlabel('Latitude ºN'); ylabel('Average of normalized change in drought');
set(gca,'FontSize',30, 'FontName', 'MS Sens Serif') % Axis Numbers and ranges Font
set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
xlim([-60,85])
ylim([-1,1])
plot([-90,90], [0,0],'LineWidth',0.5,'Color','black');
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
set(gcf,'PaperPositionMode','auto') % The most important line for proper saving the figures ***
saveas(h,dir_out_name_fig)
print(dir_out_name_fig,'-dpng','-r300')
close

dir_out_name_fig=[dir_out_fig 'Lat_Disch_Change_P95_' dir_out_name_fig_tag];
h=figure;
hold on
colors={'b','r','g','m','c','y','k'};
for gcm_i=1:size(GCM_Models_Names,2)
    Variable= (reshape (Multimodel_Q_Change_P95_hist_rcp8p5_norm_Lat(gcm_i, :, :), [Lat_n, numel(Impact_Models_Names)]))' ;
    func_shadedErrorBar(Lat',Variable,{@mean,@std},{colors{gcm_i},'LineWidth',4},1);
end
title({'Latitudinal average of normalized change in flood risk (95th percentile of discharge)' ; ('ISI-MIP - GHM: multiple (Average +/- St.Dev.) , GCM: single - [2070-2099 RCP8.5] vs. [1971-2000]')})
legend([GCM_Models_Names{1}], [GCM_Models_Names{2}], [GCM_Models_Names{3}], [GCM_Models_Names{4}], [GCM_Models_Names{5}], 'location', 'NorthWest');
xlabel('Latitude ºN'); ylabel('Average of normalized change in flood');
set(gca,'FontSize',30, 'FontName', 'MS Sens Serif') % Axis Numbers and ranges Font
set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
xlim([-60,85])
ylim([-1,1])
plot([-90,90], [0,0],'LineWidth',0.5,'Color','black');
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
set(gcf,'PaperPositionMode','auto') % The most important line for proper saving the figures ***
saveas(h,dir_out_name_fig)
print(dir_out_name_fig,'-dpng','-r300')
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Multi GHM/GCM Average by Latitude %%%
%%% GHM and GCM: multiple (Average +/- St.Dev.) %%%
dir_out_fig=[pwd '\Figures - Multimodel - RCP8.5\' ]; % Directory to save Figures and Maps
dir_out_name_fig_tag='ISIMIP_GHM_GCM_averaged_rcp8p5_hist';

Variable1 = permute(Multimodel_Q_Change_P05_hist_rcp8p5_norm_Lat,[1 3 2]);
Variable1 = reshape(Variable1,[],size(Multimodel_Q_Change_P05_hist_rcp8p5_norm_Lat,2),1) * -1;
Variable2 = permute(Multimodel_Q_Change_P95_hist_rcp8p5_norm_Lat,[1 3 2]);
Variable2 = reshape(Variable2,[],size(Multimodel_Q_Change_P95_hist_rcp8p5_norm_Lat,2),1);
dir_out_name_fig=[dir_out_fig 'Lat_Disch_Change_P05_P95_' dir_out_name_fig_tag];
h=figure;
hold on
func_shadedErrorBar(Lat',Variable1,{@mean,@std},{'r','LineWidth',4},1);
func_shadedErrorBar(Lat',Variable2,{@mean,@std},{'b','LineWidth',4},1);
title({'Latitudinal average of normalized change in drought and flood risk (5th (* -1) and 95th percentile of discharge)' ; ('ISI-MIP - multimodel average of all GHMs and GCMs - [2070-2099 RCP8.5] compared to [1971-2000]')})
legend('Drought risk (-1* P05) multimodel (average +/- st.dev.)', 'Flood risk (P95)           multimodel (average +/- st.dev.)', 'location', 'NorthWest');
xlabel('Latitude ºN'); ylabel('Average of normalized change in drought and flood');
set(gca,'FontSize',36, 'FontName', 'MS Sens Serif') % Axis Numbers and ranges Font
set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
xlim([-60,85])
ylim([-1,1])
plot([-90,90], [0,0],'LineWidth',0.5,'Color','black');
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
set(gcf,'PaperPositionMode','auto') % The most important line for proper saving the figures ***
saveas(h,dir_out_name_fig)
print(dir_out_name_fig,'-dpng','-r300')
close

Variable1 = permute(Multimodel_Q_Change_P05_hist_rcp8p5_norm_Lat,[1 3 2]);
Variable1 = reshape(Variable1,[],size(Multimodel_Q_Change_P05_hist_rcp8p5_norm_Lat,2),1) * -1;
dir_out_name_fig=[dir_out_fig 'Lat_Disch_Change_P05_' dir_out_name_fig_tag];
h=figure;
hold on
func_shadedErrorBar(Lat',Variable1,{@mean,@std},{'r','LineWidth',4},1);
title({'Latitudinal average of normalized change in drought risk (-1 * 5th percentile of discharge)' ; ('ISI-MIP - multimodel average of all GHMs and GCMs - [2070-2099 RCP8.5] compared to [1971-2000]')})
legend('Drought risk (-1* P05) multimodel (average +/- st.dev.)', 'location', 'NorthWest');
xlabel('Latitude ºN'); ylabel('Average of normalized change in drought');
set(gca,'FontSize',36, 'FontName', 'MS Sens Serif') % Axis Numbers and ranges Font
set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
xlim([-60,85])
ylim([-1,1])
plot([-90,90], [0,0],'LineWidth',0.5,'Color','black');
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
set(gcf,'PaperPositionMode','auto') % The most important line for proper saving the figures ***
saveas(h,dir_out_name_fig)
print(dir_out_name_fig,'-dpng','-r300')
close

Variable1 = permute(Multimodel_Q_Change_P95_hist_rcp8p5_norm_Lat,[1 3 2]);
Variable1 = reshape(Variable1,[],size(Multimodel_Q_Change_P95_hist_rcp8p5_norm_Lat,2),1);
dir_out_name_fig=[dir_out_fig 'Lat_Disch_Change_P95_' dir_out_name_fig_tag];
h=figure;
hold on
func_shadedErrorBar(Lat',Variable1,{@mean,@std},{'b','LineWidth',4},1);
title({'Latitudinal average of normalized change in flood risk (95th percentile of discharge)' ; ('ISI-MIP - multimodel average of all GHMs and GCMs - [2070-2099 RCP8.5] compared to [1971-2000]')})
legend('Flood risk (P95) multimodel (average +/- st.dev.)', 'location', 'NorthWest');
xlabel('Latitude ºN'); ylabel('Average of normalized change in flood');
set(gca,'FontSize',36, 'FontName', 'MS Sens Serif') % Axis Numbers and ranges Font
set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
xlim([-60,85])
ylim([-1,1])
plot([-90,90], [0,0],'LineWidth',0.5,'Color','black');
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
set(gcf,'PaperPositionMode','auto') % The most important line for proper saving the figures ***
saveas(h,dir_out_name_fig)
print(dir_out_name_fig,'-dpng','-r300')
close

Variable1 = permute(Multimodel_Q_Change_Med_hist_rcp8p5_norm_Lat,[1 3 2]);
Variable1 = reshape(Variable1,[],size(Multimodel_Q_Change_Med_hist_rcp8p5_norm_Lat,2),1);
dir_out_name_fig=[dir_out_fig 'Lat_Disch_Change_Median_' dir_out_name_fig_tag];
h=figure;
hold on
func_shadedErrorBar(Lat',Variable1,{@mean,@std},{'g','LineWidth',4},1);
title({'Latitudinal average of normalized change in median percentile of discharge' ; ('ISI-MIP - multimodel average of all GHMs and GCMs - [2070-2099 RCP8.5] compared to [1971-2000]')})
legend('P_m_e_d_i_a_n multimodel (average +/- st.dev.)', 'location', 'NorthWest');
xlabel('Latitude ºN'); ylabel('Average of normalized change in discharge');
set(gca,'FontSize',36, 'FontName', 'MS Sens Serif') % Axis Numbers and ranges Font
set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
xlim([-60,85])
ylim([-1,1])
plot([-90,90], [0,0],'LineWidth',0.5,'Color','black');
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
set(gcf,'PaperPositionMode','auto') % The most important line for proper saving the figures ***
saveas(h,dir_out_name_fig)
print(dir_out_name_fig,'-dpng','-r300')
close

Variable1 = ( reshape(nanmean(Multimodel_2_Q_Change_P05_norm,2), [Lat_n, 1]) )' * -1;
Variable2 = ( reshape(nanmean(Multimodel_2_Q_Change_P95_norm,2), [Lat_n, 1]) )';
dir_out_name_fig=[dir_out_fig 'Lat_Disch_Change_P05_P95_' dir_out_name_fig_tag '_2'];
h=figure;
hold on
plot(Lat', Variable1,'LineWidth',4,'Color','red');
plot(Lat', Variable2,'LineWidth',4,'Color','blue');
title({'Latitudinal average of normalized change in drought and flood risk (5th (* -1) and 95th percentile of discharge)' ; ('ISI-MIP - Change in multimodel-averaged discharge (all GHMs and GCMs) - [2070-2099 RCP8.5] vs. [1971-2000]')})
legend('Drought risk (-1* P05) multimodel (average +/- st.dev.)', 'Flood risk (P95)           multimodel (average +/- st.dev.)', 'location', 'NorthWest');
xlabel('Latitude ºN'); ylabel('Average of normalized change in drought and flood');
set(gca,'FontSize',36, 'FontName', 'MS Sens Serif') % Axis Numbers and ranges Font
set(findall(gcf,'type','text'),'FontSize',24, 'FontName', 'MS Sens Serif') % Text (title and axis-lable) font
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
xlim([-60,85])
ylim([-1,1])
plot([-90,90], [0,0],'LineWidth',0.5,'Color','black');
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % Maximize the figure window
set(gcf,'PaperPositionMode','auto') % The most important line for proper saving the figures ***
saveas(h,dir_out_name_fig)
print(dir_out_name_fig,'-dpng','-r300')
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Population and Land Area Affected %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlswrite(excel_out_name_1, Area_sum, 'Population_Area_Affected', 'C3')

load Population_HYDE_CIESIN_2015
Population_HYDE_CIESIN_2015=Population_HYDE_CIESIN_2015 / 1e6; % to convert to million people
Population_HYDE_CIESIN_2015(isnan(WetCells))=NaN; % excluding the cells that are not investigated in this study
Pop_sum=nansum(nansum(Population_HYDE_CIESIN_2015)); % Total population of the earth - to be used for percentage of population affected calculations
xlswrite(excel_out_name_1, Pop_sum, 'Population_Area_Affected', 'C2')

%%%% Population and Land Area - All cells %%%
Pop_Flood_inc_Drought_inc_All=Population_HYDE_CIESIN_2015; % Population experiencing increased Flood and increased Drought (indicator 1)
Pop_Flood_inc_Drought_inc_All(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==2 | Multimodel_Flood_Drought_Indicator_All==3 | Multimodel_Flood_Drought_Indicator_All==4)=NaN;
Pop_Flood_inc_Drought_dec_All=Population_HYDE_CIESIN_2015; % Population experiencing increased Flood and decresed Drought (indicator 2)
Pop_Flood_inc_Drought_dec_All(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==1 | Multimodel_Flood_Drought_Indicator_All==3 | Multimodel_Flood_Drought_Indicator_All==4)=NaN;
Pop_Flood_dec_Drought_inc_All=Population_HYDE_CIESIN_2015; % Population experiencing decresed Flood and increased Drought (indicator 3)
Pop_Flood_dec_Drought_inc_All(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==1 | Multimodel_Flood_Drought_Indicator_All==2 | Multimodel_Flood_Drought_Indicator_All==4)=NaN;
Pop_Flood_dec_Drought_dec_All=Population_HYDE_CIESIN_2015; % Population experiencing decresed Flood and decresed Drought (indicator 4)
Pop_Flood_dec_Drought_dec_All(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==1 | Multimodel_Flood_Drought_Indicator_All==2 | Multimodel_Flood_Drought_Indicator_All==3)=NaN;

xlswrite(excel_out_name_1, nansum(nansum(Pop_Flood_inc_Drought_inc_All)), 'Population_Area_Affected', 'E7')
xlswrite(excel_out_name_1, nansum(nansum(Pop_Flood_inc_Drought_dec_All)), 'Population_Area_Affected', 'D7')
xlswrite(excel_out_name_1, nansum(nansum(Pop_Flood_dec_Drought_inc_All)), 'Population_Area_Affected', 'F7')
xlswrite(excel_out_name_1, nansum(nansum(Pop_Flood_dec_Drought_dec_All)), 'Population_Area_Affected', 'G7')

Area_Flood_inc_Drought_inc_All=GridCell_Areas; % Land Area experiencing increased Flood and increased Drought (indicator 1)
Area_Flood_inc_Drought_inc_All(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==2 | Multimodel_Flood_Drought_Indicator_All==3 | Multimodel_Flood_Drought_Indicator_All==4)=NaN;
Area_Flood_inc_Drought_dec_All=GridCell_Areas; % Land Area experiencing increased Flood and decresed Drought (indicator 2)
Area_Flood_inc_Drought_dec_All(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==1 | Multimodel_Flood_Drought_Indicator_All==3 | Multimodel_Flood_Drought_Indicator_All==4)=NaN;
Area_Flood_dec_Drought_inc_All=GridCell_Areas; % Land Area experiencing decresed Flood and increased Drought (indicator 3)
Area_Flood_dec_Drought_inc_All(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==1 | Multimodel_Flood_Drought_Indicator_All==2 | Multimodel_Flood_Drought_Indicator_All==4)=NaN;
Area_Flood_dec_Drought_dec_All=GridCell_Areas; % Land Area experiencing decresed Flood and decresed Drought (indicator 4)
Area_Flood_dec_Drought_dec_All(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==1 | Multimodel_Flood_Drought_Indicator_All==2 | Multimodel_Flood_Drought_Indicator_All==3)=NaN;

xlswrite(excel_out_name_1, nansum(nansum(Area_Flood_inc_Drought_inc_All)), 'Population_Area_Affected', 'E9')
xlswrite(excel_out_name_1, nansum(nansum(Area_Flood_inc_Drought_dec_All)), 'Population_Area_Affected', 'D9')
xlswrite(excel_out_name_1, nansum(nansum(Area_Flood_dec_Drought_inc_All)), 'Population_Area_Affected', 'F9')
xlswrite(excel_out_name_1, nansum(nansum(Area_Flood_dec_Drought_dec_All)), 'Population_Area_Affected', 'G9')

%%%% Population and Land Area - Cells with agreement in half of the models on trend significance  %%%
Pop_Flood_inc_Drought_inc_Sig_half=Population_HYDE_CIESIN_2015; % Population experiencing increased Flood and increased Drought (indicator 1)
Pop_Flood_inc_Drought_inc_Sig_half(isnan(Multimodel_Flood_Drought_Indicator_Sig_half) | Multimodel_Flood_Drought_Indicator_Sig_half==2 | Multimodel_Flood_Drought_Indicator_Sig_half==3 | Multimodel_Flood_Drought_Indicator_Sig_half==4)=NaN;
Pop_Flood_inc_Drought_dec_Sig_half=Population_HYDE_CIESIN_2015; % Population experiencing increased Flood and decresed Drought (indicator 2)
Pop_Flood_inc_Drought_dec_Sig_half(isnan(Multimodel_Flood_Drought_Indicator_Sig_half) | Multimodel_Flood_Drought_Indicator_Sig_half==1 | Multimodel_Flood_Drought_Indicator_Sig_half==3 | Multimodel_Flood_Drought_Indicator_Sig_half==4)=NaN;
Pop_Flood_dec_Drought_inc_Sig_half=Population_HYDE_CIESIN_2015; % Population experiencing decresed Flood and increased Drought (indicator 3)
Pop_Flood_dec_Drought_inc_Sig_half(isnan(Multimodel_Flood_Drought_Indicator_Sig_half) | Multimodel_Flood_Drought_Indicator_Sig_half==1 | Multimodel_Flood_Drought_Indicator_Sig_half==2 | Multimodel_Flood_Drought_Indicator_Sig_half==4)=NaN;
Pop_Flood_dec_Drought_dec_Sig_half=Population_HYDE_CIESIN_2015; % Population experiencing decresed Flood and decresed Drought (indicator 4)
Pop_Flood_dec_Drought_dec_Sig_half(isnan(Multimodel_Flood_Drought_Indicator_Sig_half) | Multimodel_Flood_Drought_Indicator_Sig_half==1 | Multimodel_Flood_Drought_Indicator_Sig_half==2 | Multimodel_Flood_Drought_Indicator_Sig_half==3)=NaN;

xlswrite(excel_out_name_1, nansum(nansum(Pop_Flood_inc_Drought_inc_Sig_half)), 'Population_Area_Affected', 'E11')
xlswrite(excel_out_name_1, nansum(nansum(Pop_Flood_inc_Drought_dec_Sig_half)), 'Population_Area_Affected', 'D11')
xlswrite(excel_out_name_1, nansum(nansum(Pop_Flood_dec_Drought_inc_Sig_half)), 'Population_Area_Affected', 'F11')
xlswrite(excel_out_name_1, nansum(nansum(Pop_Flood_dec_Drought_dec_Sig_half)), 'Population_Area_Affected', 'G11')

Area_Flood_inc_Drought_inc_Sig_half=GridCell_Areas; % Land Area experiencing increased Flood and increased Drought (indicator 1)
Area_Flood_inc_Drought_inc_Sig_half(isnan(Multimodel_Flood_Drought_Indicator_Sig_half) | Multimodel_Flood_Drought_Indicator_Sig_half==2 | Multimodel_Flood_Drought_Indicator_Sig_half==3 | Multimodel_Flood_Drought_Indicator_Sig_half==4)=NaN;
Area_Flood_inc_Drought_dec_Sig_half=GridCell_Areas; % Land Area experiencing increased Flood and decresed Drought (indicator 2)
Area_Flood_inc_Drought_dec_Sig_half(isnan(Multimodel_Flood_Drought_Indicator_Sig_half) | Multimodel_Flood_Drought_Indicator_Sig_half==1 | Multimodel_Flood_Drought_Indicator_Sig_half==3 | Multimodel_Flood_Drought_Indicator_Sig_half==4)=NaN;
Area_Flood_dec_Drought_inc_Sig_half=GridCell_Areas; % Land Area experiencing decresed Flood and increased Drought (indicator 3)
Area_Flood_dec_Drought_inc_Sig_half(isnan(Multimodel_Flood_Drought_Indicator_Sig_half) | Multimodel_Flood_Drought_Indicator_Sig_half==1 | Multimodel_Flood_Drought_Indicator_Sig_half==2 | Multimodel_Flood_Drought_Indicator_Sig_half==4)=NaN;
Area_Flood_dec_Drought_dec_Sig_half=GridCell_Areas; % Land Area experiencing decresed Flood and decresed Drought (indicator 4)
Area_Flood_dec_Drought_dec_Sig_half(isnan(Multimodel_Flood_Drought_Indicator_Sig_half) | Multimodel_Flood_Drought_Indicator_Sig_half==1 | Multimodel_Flood_Drought_Indicator_Sig_half==2 | Multimodel_Flood_Drought_Indicator_Sig_half==3)=NaN;

xlswrite(excel_out_name_1, nansum(nansum(Area_Flood_inc_Drought_inc_Sig_half)), 'Population_Area_Affected', 'E13')
xlswrite(excel_out_name_1, nansum(nansum(Area_Flood_inc_Drought_dec_Sig_half)), 'Population_Area_Affected', 'D13')
xlswrite(excel_out_name_1, nansum(nansum(Area_Flood_dec_Drought_inc_Sig_half)), 'Population_Area_Affected', 'F13')
xlswrite(excel_out_name_1, nansum(nansum(Area_Flood_dec_Drought_dec_Sig_half)), 'Population_Area_Affected', 'G13')

%%%% Population and Land Area - Cells with agreement in half of the models on trend direction  %%%
Pop_Flood_inc_Drought_inc_Dir_half=Population_HYDE_CIESIN_2015; % Population experiencing increased Flood and increased Drought (indicator 1)
Pop_Flood_inc_Drought_inc_Dir_half(isnan(Multimodel_Flood_Drought_Indicator_Dir_half) | Multimodel_Flood_Drought_Indicator_Dir_half==2 | Multimodel_Flood_Drought_Indicator_Dir_half==3 | Multimodel_Flood_Drought_Indicator_Dir_half==4)=NaN;
Pop_Flood_inc_Drought_dec_Dir_half=Population_HYDE_CIESIN_2015; % Population experiencing increased Flood and decresed Drought (indicator 2)
Pop_Flood_inc_Drought_dec_Dir_half(isnan(Multimodel_Flood_Drought_Indicator_Dir_half) | Multimodel_Flood_Drought_Indicator_Dir_half==1 | Multimodel_Flood_Drought_Indicator_Dir_half==3 | Multimodel_Flood_Drought_Indicator_Dir_half==4)=NaN;
Pop_Flood_dec_Drought_inc_Dir_half=Population_HYDE_CIESIN_2015; % Population experiencing decresed Flood and increased Drought (indicator 3)
Pop_Flood_dec_Drought_inc_Dir_half(isnan(Multimodel_Flood_Drought_Indicator_Dir_half) | Multimodel_Flood_Drought_Indicator_Dir_half==1 | Multimodel_Flood_Drought_Indicator_Dir_half==2 | Multimodel_Flood_Drought_Indicator_Dir_half==4)=NaN;
Pop_Flood_dec_Drought_dec_Dir_half=Population_HYDE_CIESIN_2015; % Population experiencing decresed Flood and decresed Drought (indicator 4)
Pop_Flood_dec_Drought_dec_Dir_half(isnan(Multimodel_Flood_Drought_Indicator_Dir_half) | Multimodel_Flood_Drought_Indicator_Dir_half==1 | Multimodel_Flood_Drought_Indicator_Dir_half==2 | Multimodel_Flood_Drought_Indicator_Dir_half==3)=NaN;

xlswrite(excel_out_name_1, nansum(nansum(Pop_Flood_inc_Drought_inc_Dir_half)), 'Population_Area_Affected', 'E15')
xlswrite(excel_out_name_1, nansum(nansum(Pop_Flood_inc_Drought_dec_Dir_half)), 'Population_Area_Affected', 'D15')
xlswrite(excel_out_name_1, nansum(nansum(Pop_Flood_dec_Drought_inc_Dir_half)), 'Population_Area_Affected', 'F15')
xlswrite(excel_out_name_1, nansum(nansum(Pop_Flood_dec_Drought_dec_Dir_half)), 'Population_Area_Affected', 'G15')

Area_Flood_inc_Drought_inc_Dir_half=GridCell_Areas; % Land Area experiencing increased Flood and increased Drought (indicator 1)
Area_Flood_inc_Drought_inc_Dir_half(isnan(Multimodel_Flood_Drought_Indicator_Dir_half) | Multimodel_Flood_Drought_Indicator_Dir_half==2 | Multimodel_Flood_Drought_Indicator_Dir_half==3 | Multimodel_Flood_Drought_Indicator_Dir_half==4)=NaN;
Area_Flood_inc_Drought_dec_Dir_half=GridCell_Areas; % Land Area experiencing increased Flood and decresed Drought (indicator 2)
Area_Flood_inc_Drought_dec_Dir_half(isnan(Multimodel_Flood_Drought_Indicator_Dir_half) | Multimodel_Flood_Drought_Indicator_Dir_half==1 | Multimodel_Flood_Drought_Indicator_Dir_half==3 | Multimodel_Flood_Drought_Indicator_Dir_half==4)=NaN;
Area_Flood_dec_Drought_inc_Dir_half=GridCell_Areas; % Land Area experiencing decresed Flood and increased Drought (indicator 3)
Area_Flood_dec_Drought_inc_Dir_half(isnan(Multimodel_Flood_Drought_Indicator_Dir_half) | Multimodel_Flood_Drought_Indicator_Dir_half==1 | Multimodel_Flood_Drought_Indicator_Dir_half==2 | Multimodel_Flood_Drought_Indicator_Dir_half==4)=NaN;
Area_Flood_dec_Drought_dec_Dir_half=GridCell_Areas; % Land Area experiencing decresed Flood and decresed Drought (indicator 4)
Area_Flood_dec_Drought_dec_Dir_half(isnan(Multimodel_Flood_Drought_Indicator_Dir_half) | Multimodel_Flood_Drought_Indicator_Dir_half==1 | Multimodel_Flood_Drought_Indicator_Dir_half==2 | Multimodel_Flood_Drought_Indicator_Dir_half==3)=NaN;

xlswrite(excel_out_name_1, nansum(nansum(Area_Flood_inc_Drought_inc_Dir_half)), 'Population_Area_Affected', 'E17')
xlswrite(excel_out_name_1, nansum(nansum(Area_Flood_inc_Drought_dec_Dir_half)), 'Population_Area_Affected', 'D17')
xlswrite(excel_out_name_1, nansum(nansum(Area_Flood_dec_Drought_inc_Dir_half)), 'Population_Area_Affected', 'F17')
xlswrite(excel_out_name_1, nansum(nansum(Area_Flood_dec_Drought_dec_Dir_half)), 'Population_Area_Affected', 'G17')

%%%% Population and Land Area - Cells with agreement in two-thirds of the models on trend direction  %%%
Pop_Flood_inc_Drought_inc_Dir_twothirds=Population_HYDE_CIESIN_2015; % Population experiencing increased Flood and increased Drought (indicator 1)
Pop_Flood_inc_Drought_inc_Dir_twothirds(isnan(Multimodel_Flood_Drought_Indicator_Dir_twothirds) | Multimodel_Flood_Drought_Indicator_Dir_twothirds==2 | Multimodel_Flood_Drought_Indicator_Dir_twothirds==3 | Multimodel_Flood_Drought_Indicator_Dir_twothirds==4)=NaN;
Pop_Flood_inc_Drought_dec_Dir_twothirds=Population_HYDE_CIESIN_2015; % Population experiencing increased Flood and decresed Drought (indicator 2)
Pop_Flood_inc_Drought_dec_Dir_twothirds(isnan(Multimodel_Flood_Drought_Indicator_Dir_twothirds) | Multimodel_Flood_Drought_Indicator_Dir_twothirds==1 | Multimodel_Flood_Drought_Indicator_Dir_twothirds==3 | Multimodel_Flood_Drought_Indicator_Dir_twothirds==4)=NaN;
Pop_Flood_dec_Drought_inc_Dir_twothirds=Population_HYDE_CIESIN_2015; % Population experiencing decresed Flood and increased Drought (indicator 3)
Pop_Flood_dec_Drought_inc_Dir_twothirds(isnan(Multimodel_Flood_Drought_Indicator_Dir_twothirds) | Multimodel_Flood_Drought_Indicator_Dir_twothirds==1 | Multimodel_Flood_Drought_Indicator_Dir_twothirds==2 | Multimodel_Flood_Drought_Indicator_Dir_twothirds==4)=NaN;
Pop_Flood_dec_Drought_dec_Dir_twothirds=Population_HYDE_CIESIN_2015; % Population experiencing decresed Flood and decresed Drought (indicator 4)
Pop_Flood_dec_Drought_dec_Dir_twothirds(isnan(Multimodel_Flood_Drought_Indicator_Dir_twothirds) | Multimodel_Flood_Drought_Indicator_Dir_twothirds==1 | Multimodel_Flood_Drought_Indicator_Dir_twothirds==2 | Multimodel_Flood_Drought_Indicator_Dir_twothirds==3)=NaN;

xlswrite(excel_out_name_1, nansum(nansum(Pop_Flood_inc_Drought_inc_Dir_twothirds)), 'Population_Area_Affected', 'E19')
xlswrite(excel_out_name_1, nansum(nansum(Pop_Flood_inc_Drought_dec_Dir_twothirds)), 'Population_Area_Affected', 'D19')
xlswrite(excel_out_name_1, nansum(nansum(Pop_Flood_dec_Drought_inc_Dir_twothirds)), 'Population_Area_Affected', 'F19')
xlswrite(excel_out_name_1, nansum(nansum(Pop_Flood_dec_Drought_dec_Dir_twothirds)), 'Population_Area_Affected', 'G19')

Area_Flood_inc_Drought_inc_Dir_twothirds=GridCell_Areas; % Land Area experiencing increased Flood and increased Drought (indicator 1)
Area_Flood_inc_Drought_inc_Dir_twothirds(isnan(Multimodel_Flood_Drought_Indicator_Dir_twothirds) | Multimodel_Flood_Drought_Indicator_Dir_twothirds==2 | Multimodel_Flood_Drought_Indicator_Dir_twothirds==3 | Multimodel_Flood_Drought_Indicator_Dir_twothirds==4)=NaN;
Area_Flood_inc_Drought_dec_Dir_twothirds=GridCell_Areas; % Land Area experiencing increased Flood and decresed Drought (indicator 2)
Area_Flood_inc_Drought_dec_Dir_twothirds(isnan(Multimodel_Flood_Drought_Indicator_Dir_twothirds) | Multimodel_Flood_Drought_Indicator_Dir_twothirds==1 | Multimodel_Flood_Drought_Indicator_Dir_twothirds==3 | Multimodel_Flood_Drought_Indicator_Dir_twothirds==4)=NaN;
Area_Flood_dec_Drought_inc_Dir_twothirds=GridCell_Areas; % Land Area experiencing decresed Flood and increased Drought (indicator 3)
Area_Flood_dec_Drought_inc_Dir_twothirds(isnan(Multimodel_Flood_Drought_Indicator_Dir_twothirds) | Multimodel_Flood_Drought_Indicator_Dir_twothirds==1 | Multimodel_Flood_Drought_Indicator_Dir_twothirds==2 | Multimodel_Flood_Drought_Indicator_Dir_twothirds==4)=NaN;
Area_Flood_dec_Drought_dec_Dir_twothirds=GridCell_Areas; % Land Area experiencing decresed Flood and decresed Drought (indicator 4)
Area_Flood_dec_Drought_dec_Dir_twothirds(isnan(Multimodel_Flood_Drought_Indicator_Dir_twothirds) | Multimodel_Flood_Drought_Indicator_Dir_twothirds==1 | Multimodel_Flood_Drought_Indicator_Dir_twothirds==2 | Multimodel_Flood_Drought_Indicator_Dir_twothirds==3)=NaN;

xlswrite(excel_out_name_1, nansum(nansum(Area_Flood_inc_Drought_inc_Dir_twothirds)), 'Population_Area_Affected', 'E21')
xlswrite(excel_out_name_1, nansum(nansum(Area_Flood_inc_Drought_dec_Dir_twothirds)), 'Population_Area_Affected', 'D21')
xlswrite(excel_out_name_1, nansum(nansum(Area_Flood_dec_Drought_inc_Dir_twothirds)), 'Population_Area_Affected', 'F21')
xlswrite(excel_out_name_1, nansum(nansum(Area_Flood_dec_Drought_dec_Dir_twothirds)), 'Population_Area_Affected', 'G21')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Normalized Change Averaged by Quadrants %%%

%[Quadrant 1: P95 incrs,P05 incrs]  [Quadrant 2: P95 incrs,P05 decrs]  [Quadrant 3: P95 decrs, P05 incrs] [Quadrant 4: P95 decrs, P05 decrs]
Multimodel_3_Quadrant_Ave_P95_P05_hist_rcp8p5_norm=NaN(1,8); % [columns 1-4] = ave. norm. change in P95 in quadrants 1-4   ,   [columns 5-8] = ave. norm. change in P05 in quadrants 1-4
Var_P05=nanmean(Multimodel_Q_Change_P05_hist_rcp8p5_norm,3);
Var_P95=nanmean(Multimodel_Q_Change_P95_hist_rcp8p5_norm,3);

%%%% All cells %%%
Var_P05_help=Var_P05; Var_P95_help=Var_P95; GridCell_Areas_help=GridCell_Areas; % The Quadrant experiencing increased Flood and increased Drought (indicator 1)
Var_P05_help(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==2 | Multimodel_Flood_Drought_Indicator_All==3 | Multimodel_Flood_Drought_Indicator_All==4)=NaN; % excluding the cells that do not fit in this indicator
Var_P95_help(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==2 | Multimodel_Flood_Drought_Indicator_All==3 | Multimodel_Flood_Drought_Indicator_All==4)=NaN;
GridCell_Areas_help(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==2 | Multimodel_Flood_Drought_Indicator_All==3 | Multimodel_Flood_Drought_Indicator_All==4)=NaN;
Multimodel_3_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1,1)=nansum(nansum( Var_P95_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help )); % Cell-Area-Weighted averaging
Multimodel_3_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1,5)=nansum(nansum( Var_P05_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help ));

Var_P05_help=Var_P05; Var_P95_help=Var_P95; GridCell_Areas_help=GridCell_Areas; % The Quadrant experiencing increased Flood and decresed Drought (indicator 2)
Var_P05_help(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==1 | Multimodel_Flood_Drought_Indicator_All==3 | Multimodel_Flood_Drought_Indicator_All==4)=NaN; % excluding the cells that do not fit in this indicator
Var_P95_help(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==1 | Multimodel_Flood_Drought_Indicator_All==3 | Multimodel_Flood_Drought_Indicator_All==4)=NaN;
GridCell_Areas_help(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==1 | Multimodel_Flood_Drought_Indicator_All==3 | Multimodel_Flood_Drought_Indicator_All==4)=NaN;
Multimodel_3_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1,2)=nansum(nansum( Var_P95_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help )); % Cell-Area-Weighted averaging
Multimodel_3_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1,6)=nansum(nansum( Var_P05_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help ));

Var_P05_help=Var_P05; Var_P95_help=Var_P95; GridCell_Areas_help=GridCell_Areas; % The Quadrant experiencing decresed Flood and increased Drought (indicator 3)
Var_P05_help(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==1 | Multimodel_Flood_Drought_Indicator_All==2 | Multimodel_Flood_Drought_Indicator_All==4)=NaN; % excluding the cells that do not fit in this indicator
Var_P95_help(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==1 | Multimodel_Flood_Drought_Indicator_All==2 | Multimodel_Flood_Drought_Indicator_All==4)=NaN;
GridCell_Areas_help(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==1 | Multimodel_Flood_Drought_Indicator_All==2 | Multimodel_Flood_Drought_Indicator_All==4)=NaN;
Multimodel_3_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1,3)=nansum(nansum( Var_P95_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help )); % Cell-Area-Weighted averaging
Multimodel_3_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1,7)=nansum(nansum( Var_P05_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help ));

Var_P05_help=Var_P05; Var_P95_help=Var_P95; GridCell_Areas_help=GridCell_Areas; % The Quadrant experiencing decresed Flood and decresed Drought (indicator 4)
Var_P05_help(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==1 | Multimodel_Flood_Drought_Indicator_All==2 | Multimodel_Flood_Drought_Indicator_All==3)=NaN; % excluding the cells that do not fit in this indicator
Var_P95_help(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==1 | Multimodel_Flood_Drought_Indicator_All==2 | Multimodel_Flood_Drought_Indicator_All==3)=NaN;
GridCell_Areas_help(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==1 | Multimodel_Flood_Drought_Indicator_All==2 | Multimodel_Flood_Drought_Indicator_All==3)=NaN;
Multimodel_3_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1,4)=nansum(nansum( Var_P95_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help )); % Cell-Area-Weighted averaging
Multimodel_3_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1,8)=nansum(nansum( Var_P05_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help ));

xlswrite(excel_out_name_1, Multimodel_3_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1,1), 'Multimodel_Quadrants', 'E39')
xlswrite(excel_out_name_1, Multimodel_3_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1,2), 'Multimodel_Quadrants', 'D39')
xlswrite(excel_out_name_1, Multimodel_3_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1,3), 'Multimodel_Quadrants', 'F39')
xlswrite(excel_out_name_1, Multimodel_3_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1,4), 'Multimodel_Quadrants', 'G39')
xlswrite(excel_out_name_1, Multimodel_3_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1,5) * -1, 'Multimodel_Quadrants', 'J39')
xlswrite(excel_out_name_1, Multimodel_3_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1,6) * -1, 'Multimodel_Quadrants', 'I39')
xlswrite(excel_out_name_1, Multimodel_3_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1,7) * -1, 'Multimodel_Quadrants', 'K39')
xlswrite(excel_out_name_1, Multimodel_3_Quadrant_Ave_P95_P05_hist_rcp8p5_norm(1,8) * -1, 'Multimodel_Quadrants', 'L39')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Normalized Change Averaged by Increase/Decrease %%%

%[Column 1: P95 incrs]  [Column 2: P95 decrs]  [Column 3: P05 incrs]  [Column 4: P05 decrs]
% Multimodel_3_IndDec_Ave_P95_P05_hist_rcp8p5_norm : column 1 = P95 change average in increasing cells, column 2 = P95 change average in decreasing cells, column 3 = P05 change average in increasing cells, column 4 = P05 change average in decreasing cells
Multimodel_3_IndDec_Ave_P95_P05_hist_rcp8p5_norm=NaN(1,4);
Var_P05=nanmean(Multimodel_Q_Change_P05_hist_rcp8p5_norm,3);
Var_P95=nanmean(Multimodel_Q_Change_P95_hist_rcp8p5_norm,3);

%%%% All cells %%%
Var_P95_help=Var_P95; GridCell_Areas_help=GridCell_Areas; % Cells with Increased Flood (Quads 1 and 2)(Indicators 1 and 2)
Var_P95_help(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==3 | Multimodel_Flood_Drought_Indicator_All==4)=NaN; % excluding the cells that do not fit in this indicator
GridCell_Areas_help(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==3 | Multimodel_Flood_Drought_Indicator_All==4)=NaN;
Multimodel_3_IndDec_Ave_P95_P05_hist_rcp8p5_norm(1,1)=nansum(nansum( Var_P95_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help )); % Cell-Area-Weighted averaging

Var_P95_help=Var_P95; GridCell_Areas_help=GridCell_Areas; % Cells with Decreased Flood (Quads 3 and 4)(Indicators 3 and 4)
Var_P95_help(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==1 | Multimodel_Flood_Drought_Indicator_All==2)=NaN; % excluding the cells that do not fit in this indicator
GridCell_Areas_help(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==1 | Multimodel_Flood_Drought_Indicator_All==2)=NaN;
Multimodel_3_IndDec_Ave_P95_P05_hist_rcp8p5_norm(1,2)=nansum(nansum( Var_P95_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help )); % Cell-Area-Weighted averaging

Var_P05_help=Var_P05;GridCell_Areas_help=GridCell_Areas; % Cells with Increased Drought (Quads 1 and 3)(Indicators 1 and 3)
Var_P05_help(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==2 | Multimodel_Flood_Drought_Indicator_All==4)=NaN; % excluding the cells that do not fit in this indicator
GridCell_Areas_help(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==2 | Multimodel_Flood_Drought_Indicator_All==4)=NaN;
Multimodel_3_IndDec_Ave_P95_P05_hist_rcp8p5_norm(1,3)=nansum(nansum( Var_P05_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help )); % Cell-Area-Weighted averaging

Var_P05_help=Var_P05; GridCell_Areas_help=GridCell_Areas; % Cells with Decreased Drought (Quads 2 and 4)(Indicators 2 and 4)
Var_P05_help(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==1 | Multimodel_Flood_Drought_Indicator_All==3)=NaN; % excluding the cells that do not fit in this indicator
GridCell_Areas_help(isnan(Multimodel_Flood_Drought_Indicator_All) | Multimodel_Flood_Drought_Indicator_All==1 | Multimodel_Flood_Drought_Indicator_All==3)=NaN;
Multimodel_3_IndDec_Ave_P95_P05_hist_rcp8p5_norm(1,4)=nansum(nansum( Var_P05_help .* GridCell_Areas_help )) / nansum(nansum( GridCell_Areas_help )); % Cell-Area-Weighted averaging

xlswrite(excel_out_name_1, Multimodel_3_IndDec_Ave_P95_P05_hist_rcp8p5_norm(1,1), 'Multimodel_IncDec', 'D39')
xlswrite(excel_out_name_1, Multimodel_3_IndDec_Ave_P95_P05_hist_rcp8p5_norm(1,2), 'Multimodel_IncDec', 'E39')
xlswrite(excel_out_name_1, Multimodel_3_IndDec_Ave_P95_P05_hist_rcp8p5_norm(1,3) * -1, 'Multimodel_IncDec', 'F39')
xlswrite(excel_out_name_1, Multimodel_3_IndDec_Ave_P95_P05_hist_rcp8p5_norm(1,4) * -1, 'Multimodel_IncDec', 'G39')

toc;
