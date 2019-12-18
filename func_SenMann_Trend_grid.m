function [Q_sen_med, Z_mk, Significance_Sens_slope, Significance_Mann_Kendal, NO_st_d]=func_SenMann_Trend_grid(Data, Lat_n, Lon_n, yrs_n, min_NO_st_d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Behzad Asadieh   , Ph.D. Candidate                  %%%
%%% Civil Engineering Department - Water Resources      %%%
%%% The City College of The City University of New York %%%
%%% basadie00@citymail.cuny.edu                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Sen's slope estimator and Mann-Kendall trend test   %%%

Z_01=2.575; % Z for 1% Significancy level
Z_05=1.96; % Z for 5% Significancy level

Strt_yr=zeros(Lat_n,Lon_n); % The count of year that the data is available from, for each grid
NO_st_d=zeros(Lat_n,Lon_n); % Number of the data available for each grid
S_mk=zeros(Lat_n,Lon_n); % The Mann-Kendall test statistic S
Z_mk=NaN(Lat_n,Lon_n); % The Mann-Kendall test statistic Z
Var_s=zeros(Lat_n,Lon_n); % The Variance for Mann-Kendall method
Q_sen_med=NaN(Lat_n,Lon_n); % Sen's slope estimatr Q_med
Q_min_01=zeros(Lat_n,Lon_n); % Sen's slope estimatr Q_min 1%
Q_max_01=zeros(Lat_n,Lon_n); % Sen's slope estimatr Q_max 1%
Q_min_05=zeros(Lat_n,Lon_n); % Sen's slope estimatr Q_min 5%
Q_max_05=zeros(Lat_n,Lon_n); % Sen's slope estimatr Q_max 5%

Significance_Mann_Kendal=NaN(Lat_n,Lon_n,2); % Mann Kendal Significance - 1st layer is for 1% and 2nd layer is for 5% significancy level - Number 1 means that trend is significant
Significance_Sens_slope=NaN(Lat_n,Lon_n,2); % Sen's slope estimatr Significance - 1st layer is for 1% and 2nd layer is for 5% significancy level - Number 1 means that trend is significant

for Lt=1:Lat_n
    for Ln=1:Lon_n
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%  Mann-Kendall trend test   %%%
        
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
        
        if NO_st_d(Lt,Ln)>= 1 % That means this grid has at least 1 data available
            Significance_Mann_Kendal(Lt,Ln,:)=0;
            Significance_Sens_slope(Lt,Ln,:)=0;
        end
        
        if NO_st_d(Lt,Ln)>= min_NO_st_d %%% Minimum number of available data for the grid to have a reliable calculation %%%
            
            % The Mann-Kendall test statistic S
            for i=1:yrs_n-1
                if ~isnan(Data(Lt,Ln,i))
                    for j=i+1:yrs_n
                        if ~isnan(Data(Lt,Ln,j))
                            if (Data(Lt,Ln,j)-Data(Lt,Ln,i)) > 0
                                S_mk(Lt,Ln)=S_mk(Lt,Ln)+1;
                            elseif (Data(Lt,Ln,j)-Data(Lt,Ln,i)) < 0
                                S_mk(Lt,Ln)=S_mk(Lt,Ln)-1;
                            end
                        end
                    end
                end
            end
            
            Var_s(Lt,Ln)=(NO_st_d(Lt,Ln)*(NO_st_d(Lt,Ln)-1)*(2*NO_st_d(Lt,Ln)+5))/18; % Variance(s) of the data
            
            % The Mann-Kendall test statistic Z
            if S_mk(Lt,Ln)>=0
                Z_mk(Lt,Ln)=(S_mk(Lt,Ln)-1)/(sqrt(Var_s(Lt,Ln)));
            elseif S_mk(Lt,Ln)<0
                Z_mk(Lt,Ln)=(S_mk(Lt,Ln)+1)/(sqrt(Var_s(Lt,Ln)));
            end
            
            %%% Mann-Kendall Significance Assessment %%%
            if abs(Z_mk(Lt,Ln)) > Z_01
                if Z_mk(Lt,Ln) > Z_01
                    Significance_Mann_Kendal(Lt,Ln,1)=1;
                end
                if Z_mk(Lt,Ln) < (-1*Z_01)
                    Significance_Mann_Kendal(Lt,Ln,1)=-1;
                    
                end
            end
            if abs(Z_mk(Lt,Ln)) > Z_05
                if Z_mk(Lt,Ln) > Z_05
                    Significance_Mann_Kendal(Lt,Ln,2)=1;
                end
                if Z_mk(Lt,Ln) < (-1*Z_05)
                    Significance_Mann_Kendal(Lt,Ln,2)=-1;
                    
                end
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Sen's slope estimator %%%
            
            n_sen=(NO_st_d(Lt,Ln)*(NO_st_d(Lt,Ln)-1)/2);
            q_sen_0=zeros(n_sen,1); % Un-Sorted Q_sen matrix
            counter=1;
            for i=1:yrs_n-1
                if ~isnan(Data(Lt,Ln,i))
                    for j=i+1:yrs_n
                        if ~isnan(Data(Lt,Ln,j))
                            q_sen_0(counter,1)=(Data(Lt,Ln,j)-Data(Lt,Ln,i))/(j-i);
                            counter=counter+1;
                        end
                    end
                end
            end
            
            q_sen=sort(q_sen_0); % Sorted Q_sen matrix
            Q_sen_med(Lt,Ln)=median(q_sen);
            
            %%% alpha=0.01 Significance %%%
            c_alpha01=Z_01*sqrt(Var_s(Lt,Ln)); % Confidence Interval (Z(1-a/2)*sqrt(Var(s)))
            m_1_01=(n_sen - c_alpha01)/2;
            m_2_01=(n_sen + c_alpha01)/2;
            Q_min_01(Lt,Ln)=q_sen((n_sen - ceil(m_1_01)),1); % Q_min is (M1)th largest of N-ordered slope estimates
            Q_max_01(Lt,Ln)=q_sen((n_sen - ceil(m_2_01)+1),1); % Q_max is (M2+1)th largest of N-ordered slope estimates
            
            if Q_sen_med(Lt,Ln) > 0 % That means there is a positive trend in data
                if (Q_min_01(Lt,Ln)*Q_max_01(Lt,Ln)) > 0 % that means Q_min and Q_max have the same sign and the trend is significant
                    Significance_Sens_slope(Lt,Ln,1)=1;
                end
            elseif Q_sen_med(Lt,Ln) < 0 % That means there is a negetive trend in data
                if (Q_min_01(Lt,Ln)*Q_max_01(Lt,Ln)) > 0 % that means Q_min and Q_max have the same sign and the trend is significant
                    Significance_Sens_slope(Lt,Ln,1)=-1;
                end
            end
            
            %%% alpha=0.05 Significance %%%
            c_alpha05=Z_05*sqrt(Var_s(Lt,Ln)); % Confidence Interval (Z(1-a/2)*sqrt(Var(s)))
            m_1_05=(n_sen - c_alpha05)/2;
            m_2_05=(n_sen + c_alpha05)/2;
            Q_min_05(Lt,Ln)=q_sen((n_sen - ceil(m_1_05)),1); % Q_min is (M1)th largest of N-ordered slope estimates
            Q_max_05(Lt,Ln)=q_sen((n_sen - ceil(m_2_05)+1),1); % Q_max is (M2+1)th largest of N-ordered slope estimates
            
            if Q_sen_med(Lt,Ln) > 0 % That means there is a positive trend in data
                if (Q_min_05(Lt,Ln)*Q_max_05(Lt,Ln)) > 0 % that means Q_min and Q_max have the same sign and the trend is significant
                    Significance_Sens_slope(Lt,Ln,2)=1;
                end
            elseif Q_sen_med(Lt,Ln) < 0 % That means there is a negetive trend in data
                if (Q_min_05(Lt,Ln)*Q_max_05(Lt,Ln)) > 0 % that means Q_min and Q_max have the same sign and the trend is significant
                    Significance_Sens_slope(Lt,Ln,2)=-1;
                end
            end
            
            
        end
        
    end
end


end

