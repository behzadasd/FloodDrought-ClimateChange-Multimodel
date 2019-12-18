function [b_slope, Data_ave, P_value, Reg_Significance, NO_st_d]=func_LinReg1_Trend_grid_P(Data, Lat_n, Lon_n, yrs_n, min_NO_st_d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Behzad Asadieh   , Ph.D. Candidate                  %%%
%%% Civil Engineering Department - Water Resources      %%%
%%% The City College of The City University of New York %%%
%%% basadie00@citymail.cuny.edu                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%        Linear Regression Trend Analysis         %%%
%%%           P-Value included in outputs           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Strt_yr=zeros(Lat_n,Lon_n); % The count of year that the data is available from, for each station
NO_st_d=zeros(Lat_n,Lon_n); % Number of the data available for each station

b_slope=NaN(Lat_n,Lon_n);
Reg_Significance=NaN(Lat_n,Lon_n,2); % Linear Regression trend Significance at 1% and 5% level (1st and 2nd layer, respectively)
Data_ave=NaN(Lat_n,Lon_n); % Average Data of each station
P_value=NaN(Lat_n,Lon_n); % P-value of the linear regression - Represents the significance of the trend
% % bint_05=NaN(Lat_n,Lon_n,2,2); % Intervals of b at 95% confidence (5% Significance Level)
% % bint_01=NaN(Lat_n,Lon_n,2,2); % Intervals of b at 99% confidence (1% Significance Level)
% % All_stast=NaN(Lat_n,Lon_n,4);
% % All_Z1_05=NaN(Lat_n,Lon_n); % Z1= b / delta(b) at 5% Level
% % All_Z1_01=NaN(Lat_n,Lon_n); % Z1= b / delta(b) at 1% Level

for Lt=1:Lat_n
    for Ln=1:Lon_n
        
        for i=1:yrs_n
            if ~isnan(Data(Lt,Ln,i))
                Strt_yr(Lt,Ln)=i; % Gives the count of year that the data of this grid has been started from
                break
            end
        end
        
        if Strt_yr(Lt,Ln)>0 % That means this grid has at least 1 data available
            counter=0;
            for i=Strt_yr(Lt,Ln):yrs_n
                if ~isnan(Data(Lt,Ln,i))
                    counter=counter+1;
                end
            end
            NO_st_d(Lt,Ln)=counter; % Number of the data available for this grid
        end
        
        if NO_st_d(Lt,Ln)>= 1 % That means this station has at least 1 data available
            Reg_Significance(Lt,Ln,:)=0;
        end
        
        if NO_st_d(Lt,Ln)>= min_NO_st_d %%% Minimum number of available data for the station to have a reliable calculation %%%
            
            t_vec=zeros(NO_st_d(Lt,Ln),1);
            counter_a=1;
            counter_b=1;
            for i=1:yrs_n
                if ~isnan(Data(Lt,Ln,i))
                    t_vec(counter_b,1)=counter_a;
                    counter_b=counter_b+1;
                end
                counter_a=counter_a+1;
            end
            
            help_Data=zeros(NO_st_d(Lt,Ln),1); % Gathers all non-NaN values of the current station in a vector
            counter=1;
            for i=1:yrs_n
                if ~isnan(Data(Lt,Ln,i))
                    help_Data(counter,1)=Data(Lt,Ln,i);
                    counter=counter+1;
                end
            end
            
            x_vec=[ones(NO_st_d(Lt,Ln),1) t_vec];
                       
            %%% 5% Level Calculations %%%
% %             [b,bint,~,~,stats] = regress(help_Data,x_vec);
            [b,~,~,~,stats] = regress(help_Data,x_vec);
            
            b_slope(Lt,Ln)=b(2,1);
            P_value(Lt,Ln)=stats(1,3);
% %             bint_05(Lt,Ln,:,1)=bint(1,:);
% %             bint_05(Lt,Ln,:,2)=bint(2,:);
% %             All_stast(Lt,Ln,:)=stats;
            Data_ave(Lt,Ln)=mean(help_Data);
            
% %             All_Z1_05(Lt,Ln)=b(2,1)/((abs(bint(2,2)-bint(2,1)))/2);
            
            if ( stats(1,3) < 0.05 && b(2,1)>=0 )
                Reg_Significance(Lt,Ln,2)=1;
            elseif ( stats(1,3) < 0.05 && b(2,1)<0 )
                Reg_Significance(Lt,Ln,2)=-1;
            end
            
            %%% 1% Level Calculations %%%
% %             [b,bint,~,~,stats] = regress(help_Data,x_vec,0.01);
            [b,~,~,~,stats] = regress(help_Data,x_vec,0.01);
            
% %             bint_01(Lt,Ln,:,1)=bint(1,:);
% %             bint_01(Lt,Ln,:,2)=bint(2,:);          
% %             All_Z1_01(Lt,Ln)=b(2,1)/((abs(bint(2,2)-bint(2,1)))/2);
            
            if ( stats(1,3) < 0.01 && b(2,1)>=0 )
                Reg_Significance(Lt,Ln,1)=1;
            elseif ( stats(1,3) < 0.01 && b(2,1)<0 )
                Reg_Significance(Lt,Ln,1)=-1;
            end
            
        end
        
    end
end


end

