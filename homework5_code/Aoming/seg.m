function segments =seg(data1,data2,data3)
%input segmented data of three years
%This function combines the three years segmented data

n = 6*24*60;

s1 = buffer(data1,n,0.5*n);
s1 = s1(:,2:end);

s2 = buffer(data2,n,0.5*n);
s2 = s2(:,2:end);

s3 = buffer(data3,n,0.5*n);
s3 = s3(:,2:end);

segments = [s1 s2 s3];
end

