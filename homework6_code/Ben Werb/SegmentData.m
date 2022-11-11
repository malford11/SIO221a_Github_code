function [SegmentedData] = SegmentData(Data,SegmentLength)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
SegmentedData = struct;

DataNames = fieldnames(Data);
n = length(DataNames) - 1; %Subtract Data.Names

for i = 1:n
    ncols = length(Data.(char(DataNames(i)))) / SegmentLength;
    x1 = reshape(Data.(char(DataNames(i))),SegmentLength,ncols);
    x2 = reshape(Data.(char(DataNames(i)))(SegmentLength/2 + 1:end - SegmentLength/2),SegmentLength,ncols-1);
    SegmentedData.(char(DataNames(i))) = [x1 x2];
end

end

