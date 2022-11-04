function [Chunks] = Chunkify(df)
%CHUNKIFY Summary of this function goes here
%   Detailed explanation goes here
Chunks = struct;
%10 minute samples in 60 days:
%WSPD1
n = (24*60*60)/10;
length(df.WSPD1)/n;
days3090 = n*5.5;
days60 = n*5;
test = df.WSPD1(1:days60);
x1 = reshape(test,n,[]);
test2 = df.WSPD1((n/2)+1:days3090);
x2 = reshape(test2,8640,[]);
Chunks.WSPD1 = [x1 x2];
%VWND1
test = df.VWND1(1:days60);
x1 = reshape(test,n,[]);
test2 = df.VWND1((n/2)+1:days3090);
x2 = reshape(test2,8640,[]);
Chunks.VWND1 = [x1 x2];
%UWND1
test = df.UWND1(1:days60);
x1 = reshape(test,n,[]);
test2 = df.UWND1((n/2)+1:days3090);
x2 = reshape(test2,8640,[]);
Chunks.UWND1 = [x1 x2];
%WSPD2
length(df.WSPD2)/n;
days60 = n*7;
days3090=n*6;
test = df.WSPD2(1:days60);
x1 = reshape(test,n,[]);
test2 = df.WSPD2((n/2)+1:days3090+(n/2))
x2 = reshape(test2,8640,[]);
Chunks.WSPD2 = [x1 x2];
%VWND2
test = df.VWND2(1:days60);
x1 = reshape(test,n,[]);
test2 = df.VWND2((n/2)+1:days3090+(n/2));
x2 = reshape(test2,8640,[]);
Chunks.VWND2 = [x1 x2];
%VWND2
test = df.UWND2(1:days60);
x1 = reshape(test,n,[]);
test2 = df.UWND2((n/2)+1:days3090+(n/2));
x2 = reshape(test2,8640,[]);
Chunks.UWND2 = [x1 x2];
%WSPD3
length(df.WSPD3)/n;
days60 = n*5;
days3090=n*4;
test = df.WSPD3(1:days60);
x1 = reshape(test,n,[]);
test2 = df.WSPD3((n/2)+1:days3090+(n/2));
x2 = reshape(test2,8640,[]);
Chunks.WSPD3 = [x1 x2];
%VWND3
test = df.VWND3(1:days60);
x1 = reshape(test,n,[]);
test2 = df.VWND3((n/2)+1:days3090+(n/2));
x2 = reshape(test2,8640,[]);
Chunks.VWND3 = [x1 x2];
%VWND3
test = df.UWND3(1:days60);
x1 = reshape(test,n,[]);
test2 = df.UWND3((n/2)+1:days3090+(n/2));
x2 = reshape(test2,8640,[]);
Chunks.UWND3 = [x1 x2];

%Multi-Year Chunks
Chunks.WSPD = [Chunks.WSPD1 Chunks.WSPD2 Chunks.WSPD3];
Chunks.VWND = [Chunks.VWND1 Chunks.VWND2 Chunks.VWND3];
Chunks.UWND = [Chunks.UWND1 Chunks.UWND2 Chunks.UWND3];
end

