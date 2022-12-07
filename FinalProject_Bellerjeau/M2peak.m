function [peakmax] = M2peak(specmat,k,n)

% M2 Peak
a = find(k>1.52,1,'first');
b = find(k>2.44,1,'first');
kM2 = k(a:b);
for i = 1:n
    peakchunk(:,i) = specmat(a:b,i);
    peakmax(i) = max(specmat(a:b,i));
end

% %%%%verifying that we found the correct peak locations
% figure()
% for j = 1:length(mp3.z)
%     plot3(kM2,mp3.z(j)*ones(length(kM2),1),log10(M2all(:,j)),'Color',DepthSpectra(j,:))
%     plot3(1.95,mp3.z(j),log10(M2peaks(j)),'ko')
%     hold on
% end
% 
% figure()
% plot(M2peaks,mp3.z,'-k')
% axis ij
% ylabel('Depth')
% xlabel('Power')
% title('M2 Peak over Depth')

