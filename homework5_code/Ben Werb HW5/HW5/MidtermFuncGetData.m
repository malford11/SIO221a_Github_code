function [NonChunked,Chunked] = MidtermFuncGetData()
%Get and process data for midterm func
df = Structify('OS_T8S110W_DM134A-20150425_D_WIND_10min.nc' ...
    ,'OS_T8S110W_DM183A-20160321_D_WIND_10min.nc', ...
    'OS_T8S110W_DM231A-20170606_D_WIND_10min.nc'); %all data
NonChunked = FillGaps(df); %data with gaps filled
Chunked = Chunkify(NonChunked);
end

