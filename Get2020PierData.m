function pier=Get2020PierData(date_start,date_end,file)
%function pier=Get2020PierData(date_start,date_end,file)
%Return temperature and pressure and time (datenum format) from the 2020
%pier record.  
%
%
%ALB change 09/22/2022
%MHA SIO211a
%8/3/2022
% MHA: v2: include option to specify file path and name.

if nargin < 3
file='/Users/malford/GoogleDrive/Work/Projects/Teaching/sio221a/FromSarah/homework_2021/data/scripps_pier-2020.nc';
end


%
%file='/Users/malford/GoogleDrive/Work/Projects/Teaching/sio221a/FromSarah/homework_2021/data/scripps_pier-2020.nc';
time=ncread(file,'time');
temperature=ncread(file,'temperature');
pressure=ncread(file,'pressure');

date0=datenum(1970,1,1);
dnum=double(time)/3600/24+date0;


if nargin < 1
    date_start=min(dnum);
    date_end=max(dnum);
end

%i1=find(dnum > datenum(2020,6,1,0,0,0) & dnum < datenum(2020,6,7,1,1,1));
i1=find(dnum > date_start & dnum < date_end);

pier.dnum=dnum(i1);
pier.temperature=temperature(i1);
pier.pressure=pressure(i1);
pier.readme='2020 Pier data, SIO221a, function Get2020PierData.m';


%plot(dnum(i1),pressure(i1),'LineWidth',1)
%datetick('x')
%set(gca,'FontSize',16)
%xlabel('months (of 2020)','FontSize',16)
%label('Pressure (decibars)', 'FontSize',16)

%pier will then be returned.