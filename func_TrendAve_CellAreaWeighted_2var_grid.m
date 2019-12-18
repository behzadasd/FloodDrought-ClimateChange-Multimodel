function [C_Ave_Var1, C_Ave_Var2]=func_TrendAve_CellAreaWeighted_2var_grid(Var1, Var2, Lat, Lon, Lat_bound, Lon_bound, NO_st_d, min_NO_st_d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Behzad Asadieh                                      %%%
%%% Civil Engineering Department - Water Resources      %%%
%%% The City College of The City University of New York %%%
%%% basadie00@citymail.cuny.edu                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Trend test results Averaging - Cell Area Weighted - Globally and by Continent   %%%

Continents_LatLong_Bound=[12,90,-180,-52.5;-60,11.999,-100,-30;-50,-10,110,180;42,90,-52.499,60;36,41.999,-52.499,26;42,90,60.001,180;36,41.999,26.001,180;30,35.999,34,73.999;30,35.999,79.5,180; ...
    11.999,29.999,34,70;-10,30,88.5,180;30,35.999,74,79.499;0,29.999,70,88.499;-50,35.999,-30,34;-50,11.999,34.001,60];

% Area = R^2 * ( Lon2 - Lon1 ) * ( Sin(Lat2) - Sin(Lat1) )
earth_R = 6378; % Earth Radius - Unit is kilometer (km)
lon_diff_miltiplier = ( Lon_bound(:,2) - Lon_bound(:,1) )' * (pi/180) ; % ( Lon2 - Lon1 ) in the formula
lat_sin_multiplier =  sind(Lat_bound(:,1)) - sind(Lat_bound(:,2))  ;

GridCell_Areas1=NaN(size(Lat_bound,1), size(Lon_bound,1));
for i=1:size(Lat_bound,1)
    for j=1:size(Lon_bound,1)
        
        GridCell_Areas1 (i,j) = abs( (earth_R^2) * lon_diff_miltiplier (1,j) * lat_sin_multiplier (i,1) );
        
    end
end
GridCell_Areas2=GridCell_Areas1;

if nargin == 6 % Number of input arguments is less than 7 - That means the NO_st_d and min_NO_st_d are not given and so do not play important role in the process
    NO_st_d=ones(size(Var1,1), size(Var1,2)) *10;
    NO_st_d(isnan(Var1))=NaN;
    min_NO_st_d=1;
end

% Columns are: Continents are North America, South America, Oceania, Europe1, Europe2, Asia 1-6, India1, India2, Africa 1, Africa 2
% India will be included in Asia average as well
% Rows are: Lat_Down, Lat_Up, Lon_Left, Lon_Right

GridCell_Areas1(isnan(Var1))=NaN;
GridCell_Areas2(isnan(Var2))=NaN;

sum_Var1=zeros(size(Continents_LatLong_Bound,1),1);
sum_Var2=zeros(size(Continents_LatLong_Bound,1),1);
area_accumulator=zeros(size(Continents_LatLong_Bound,1),1); % number of stations for each continent

for Lt=1:size(Lat,1)
    for Ln=1:size(Lon,1)
        
        if NO_st_d(Lt,Ln)>= min_NO_st_d
            
            for i=1:size(sum_Var1,1)
                
                if (Lat(Lt,1) > Continents_LatLong_Bound(i,1)) && (Lat(Lt,1) <= Continents_LatLong_Bound(i,2)) && (Lon(Ln,1) > Continents_LatLong_Bound(i,3)) && (Lon(Ln,1) <= Continents_LatLong_Bound(i,4))
                    
                    if ~isnan(Var1(Lt,Ln))
                        area_accumulator(i,1)=area_accumulator(i,1)+GridCell_Areas1 (Lt,Ln);
                        sum_Var1(i,1)=sum_Var1(i,1)+Var1(Lt,Ln) * GridCell_Areas1 (Lt,Ln);
                    end
                    if ~isnan(Var2(Lt,Ln))
                        sum_Var2(i,1)=sum_Var2(i,1)+Var2(Lt,Ln) * GridCell_Areas2 (Lt,Ln);
                    end
                    
                end
                
            end
            
        end
        
    end
end

% Global, North America, South America, Europe, Oceania, Africa, Asia, India
C_Ave_Var1=zeros(8,1);

C_Ave_Var1(1,1)=nansum(nansum( Var1 .* GridCell_Areas1 )) / nansum(nansum( GridCell_Areas1 )); % Global
C_Ave_Var1(2,1)=sum_Var1(1,1)/area_accumulator(1,1); % North America
C_Ave_Var1(3,1)=sum_Var1(2,1)/area_accumulator(2,1); % South America
C_Ave_Var1(4,1)=sum(sum_Var1(4:5,1))/sum(area_accumulator(4:5,1)); % Europe
C_Ave_Var1(5,1)=sum_Var1(3,1)/area_accumulator(3,1); % Oceania
C_Ave_Var1(6,1)=sum(sum_Var1(14:15,1))/sum(area_accumulator(14:15,1)); % Africa
C_Ave_Var1(7,1)=sum(sum_Var1(6:13,1))/sum(area_accumulator(6:13,1)); % Asia
C_Ave_Var1(8,1)=sum(sum_Var1(12:13,1))/sum(area_accumulator(12:13,1)); % India

C_Ave_Var2=zeros(8,1);

C_Ave_Var2(1,1)=nansum(nansum( Var2 .* GridCell_Areas2 )) / nansum(nansum( GridCell_Areas2 )); % Global
C_Ave_Var2(2,1)=sum_Var2(1,1)/area_accumulator(1,1); % North America
C_Ave_Var2(3,1)=sum_Var2(2,1)/area_accumulator(2,1); % South America
C_Ave_Var2(4,1)=sum(sum_Var2(4:5,1))/sum(area_accumulator(4:5,1)); % Europe
C_Ave_Var2(5,1)=sum_Var2(3,1)/area_accumulator(3,1); % Oceania
C_Ave_Var2(6,1)=sum(sum_Var2(14:15,1))/sum(area_accumulator(14:15,1)); % Africa
C_Ave_Var2(7,1)=sum(sum_Var2(6:13,1))/sum(area_accumulator(6:13,1)); % Asia
C_Ave_Var2(8,1)=sum(sum_Var2(12:13,1))/sum(area_accumulator(12:13,1)); % India

end
