%this funciton is used to load time and temepraturte data with given years 
%years is a vector 
%output of t_date and temp are cell array 

function [t_date,temp,p]=load_data(years)

n_year = length(years); %find how many years of data we are pulling out 
min_year= min(years); %find the data starting year 
date0=datenum(1970,1,1); 

%t_total = []; %prelocate time data 
%temp_total=[]; %prelocate temperature data
for ii = 1:n_year
    k = years(ii); % create index that corr sponds to years 
    filein = ['scripps_pier-' num2str(k) '.nc'];

    t{ii} = ncread(filein, 'time');
    t_date{ii} = double(t{ii})/3600/24+date0; %convert time into date 
    t_date{ii} = datetime(t_date{ii},'ConvertFrom','datenum'); % convert datetime 
    %t_total = [t_total;t{ii}];
    
    temp{ii}= ncread(filein, 'temperature');
    maxT_ind = find(temp{ii}>50); 
    temp{ii}(maxT_ind) = nan;% filter out data that is impossibly high
    %temp_total = [temp_total;temp{ii}];
    maxT_ind = []; %clear temp index for next loop 

    p{ii}= ncread(filein, 'pressure');
    %maxT_ind = find(temp{ii}>50); 
    %temp{ii}(maxT_ind) = nan;% filter out data that is impossibly high
    %temp_total = [temp_total;temp{ii}];
    %maxT_ind = []; %clear temp index for next loop 
end 

end 