function pier=MultiYearPierData(path, yrs)
%function for loading multiyear pier data from the Scripps pier
%path should have all the scripps pier data from 2005-2021
%takes vector with year values as input and gives a struct with timeUTC, temp,
%pressure as output
%Ishwari 10/10/22

%if path is not specified
if nargin < 2
filePath='C:\Users\Ishwari\OneDrive\Desktop\SIOC_221A\Data\';
end

%creating an array of filenames
filePath = path;                        %path for data
fileName = strings(1,length(yrs));      %pre-allocating 
fileIn = strings(1,length(yrs));        %pre-allocating
timeUTC = datetime([], [], []);         %pre-allocating
temp = [];                              %pre-allocating
pressure = [];                          %pre-allocating

%loop for creating a time and a temp array for the 2005-2021 period
for i = 1:length(yrs)
    %creating string array of filenames
    fileName(i) = string(['scripps_pier-' num2str(yrs(i)) '.nc']);
    fileIn(i) = strcat(filePath,fileName(i)); 
    
    %reading in temperature for each year
    temp1 = ncread(fileIn(i), 'temperature'); 
    n = length(temp);
    %appending temperature values for each year to 'temp' array
    temp(n+1:n+length(temp1)) = temp1; 
    
    %reading in pressure data for each year
    p1 = ncread(fileIn(i), 'pressure'); 
    n = length(pressure);
    %appending pressure values for each year to 'temp' array
    pressure(n+1:n+length(p1)) = p1; 
    
    %reading in time data for each year
    time = ncread(fileIn(i), 'time'); 
    %converting to datetime
    date0=datenum(1970,1,1); 
    dnum=double(time)/(3600*24) + date0;
    n = length(timeUTC);
    %appending datetime values for each year to 'timeUTC' array
    timeUTC(n+1:n+length(time))  = dnum - min(dnum) + datetime(yrs(i), 1, 1, 00, 00, 00); 
end

%removing anomalously hightemp data above 40 C
temp(find(temp >= 40)) = NaN; 

%assigning variables in structure
pier.time = timeUTC; 
pier.temp = temp;
pier.pressure = pressure;
end