function [out] = getpierdata_kk(file)
%input -- file: file, including full location, as a string
%output -- time, temperature, pressure
    out.time=ncread(file,'time'); %from meta: time given in seconds from 1/1/1970
    out.date0=datenum(1970,1,1); %defined start time
    out.dnum=double(time)/3600/24 + date0; %convert time to datenum form
    out.temperature=ncread(file,'temperature'); %read temperature data
    out.pressure=ncread(file,'pressure'); %read pressure data
end

