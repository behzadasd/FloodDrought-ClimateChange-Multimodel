function [Color_Value]=func_FloodDroughtRisk_indicator(Var_P05, Var_P95)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Behzad Asadieh   , Ph.D. Candidate                  %%%
%%% Civil Engineering Department - Water Resources      %%%
%%% The City College of The City University of New York %%%
%%% basadie00@citymail.cuny.edu                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lat_n=size(Var_P05,1);
Lon_n=size(Var_P05,2);

Color_Value=NaN(Lat_n, Lon_n);

for lt=1:Lat_n
    for ln=1:Lon_n
        if ~isnan(Var_P05(lt,ln)) && ~isnan(Var_P95(lt,ln))
            
            %%% Flood and Drought risk increase %%%
            if (Var_P95(lt,ln)>=0 && Var_P95(lt,ln)<= 1) && (Var_P05(lt,ln)<= 0 && Var_P05(lt,ln)>= -1) % Line 1 ColorMap - Purple
                Color_Value(lt,ln) = 1;
            end
            
            %%% Flood risk increase, Drought risk decrease %%%
            if (Var_P95(lt,ln)>=0 && Var_P95(lt,ln)<= 1) && (Var_P05(lt,ln)<= 1 && Var_P05(lt,ln)> 0) % Line 2 ColorMap - Blue
                Color_Value(lt,ln) = 2;
            end
            
            %%% Drought risk increase, Flood risk decrease %%%
            if (Var_P95(lt,ln)>= -1 && Var_P95(lt,ln)< 0) && (Var_P05(lt,ln)<= 0 && Var_P05(lt,ln)>= -1) % Line 3 ColorMap - Red
                Color_Value(lt,ln) = 3;
            end
            
            %%% Flood and Drought risk decrease %%%
            if (Var_P95(lt,ln)>= -1 && Var_P95(lt,ln)< 0) && (Var_P05(lt,ln)<= 1 && Var_P05(lt,ln)> 0) % Line 4 ColorMap - Green
                Color_Value(lt,ln) = 4;
            end
            
        end
    end
end


end

