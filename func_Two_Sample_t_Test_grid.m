function [T_test_value, T_Significance]=func_Two_Sample_t_Test_grid(Data_a, Data_b, min_NO_st_d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Behzad Asadieh   , Ph.D. Candidate                  %%%
%%% Civil Engineering Department - Water Resources      %%%
%%% The City College of The City University of New York %%%
%%% basadie00@citymail.cuny.edu                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Two Sample t-test - Data_a compared to Data_b %%%

if nargin == 2 % Number of input arguments is less than 7 - That means the NO_st_d and min_NO_st_d are not given and so do not play important role in the process
    min_NO_st_d=2;
end

lat_nn=size(Data_a,1); % Number of Latitudinal elements
lon_nn=size(Data_a,2); % Number of Longitudinal elements

NO_st_d_a=sum(~isnan(Data_a),3); NO_st_d_a(NO_st_d_a==0)=NaN;
NO_st_d_b=sum(~isnan(Data_b),3); NO_st_d_b(NO_st_d_b==0)=NaN;

N=floor( nanmean(nanmean(NO_st_d_a)) );
alpha=0.05; T_05=tinv(1-(alpha/2),N);

T_test_value=NaN(lat_nn,lon_nn);  % T value from Two Sample t-test 

for Lt=1:lat_nn
    for Ln=1:lon_nn
               
        if NO_st_d_a(Lt,Ln)>= min_NO_st_d %%% Minimum number of available data for the grid to have a reliable calculation %%%
            
            aa=reshape(Data_a(Lt,Ln,:), [size(Data_a,3),1]);
            bb=reshape(Data_b(Lt,Ln,:), [size(Data_b,3),1]);
            
            M_a=nanmean(aa,1); % Average of Sample_a
            M_b=nanmean(bb,1); % Average of Sample_b
            
            N_a=NO_st_d_a(Lt,Ln); % Size of Sample_a
            N_b=NO_st_d_b(Lt,Ln); % Size of Sample_b
            
            Var_a=1/(N_a-1) * sum((aa - M_a).^2); % Sample Variance of Sample_a
            Var_b=1/(N_b-1) * sum((bb - M_b).^2); % Sample Variance of Sample_a
            
            T_test_value(Lt,Ln)=( M_a - M_b ) / ( sqrt ( (Var_a/N_a) + (Var_b/N_b) ) );
            
            
        end
        
    end
end


T_Significance=T_test_value; % Two Sample t-test Significance at 95% confidence level - Number 1 or -1 means that trend is significant
T_Significance(~isnan(T_test_value))=0;
T_Significance(T_test_value>T_05)=1;
T_Significance(T_test_value<-T_05)=-1;

% % T_Significance=NaN(lat_nn,lon_nn,3); % Two Sample t-test Significance - 1st layer is for 99%, 2nd layer is for 95% and 3rd layer is for 90% significancy level - Number 1 or -1 means that trend is significant
% % 
% % N=floor( nanmean(nanmean(NO_st_d_a)) );
% % alpha=0.01; T_01=tinv(1-(alpha/2),N);
% % alpha=0.05; T_05=tinv(1-(alpha/2),N);
% % alpha=0.10; T_10=tinv(1-(alpha/2),N);
% % 
% % T_S01=T_test_value;
% % T_S01(~isnan(T_test_value))=0;
% % T_S01(T_test_value>T_01)=1;
% % T_S01(T_test_value<-T_01)=-1;
% % 
% % T_S05=T_test_value;
% % T_S05(~isnan(T_test_value))=0;
% % T_S05(T_test_value>T_05)=1;
% % T_S05(T_test_value<-T_05)=-1;
% % 
% % T_S10=T_test_value;
% % T_S10(~isnan(T_test_value))=0;
% % T_S10(T_test_value>T_10)=1;
% % T_S10(T_test_value<-T_10)=-1;
% % 
% % T_Significance(:,:,1)=T_S01;
% % T_Significance(:,:,2)=T_S05;
% % T_Significance(:,:,3)=T_S10;

end

