%% randomly produce samples
clear;clc;
[MT,Rmt] = geotiffread('J:\xinyue\MOUNTAIN\MT_re.tif');
fs=dir('J:\xinyue\HANSEN\v1.6\*.tif');     
lat_lon_yr_save=[];
for pp =1:length(fs); 
    tic;
    disp(['begin pixel ', num2str(pp)])
    fn = [char(fs(pp).folder),'/',char(fs(pp).name)];   
    if fn(end-9) == 'N' 
        lat = str2double(fn(end-11:end-10));
    elseif fn(end-9) == 'S'
        lat = -str2double(fn(end-11:end-10));
    end    
    if fn(end-4) == 'E'
        lon = str2double(fn(end-7:end-5));
    elseif fn(end-4) == 'W'
        lon = -str2double(fn(end-7:end-5));
    end 
    lat_start = (Rmt.LatitudeLimits(2)-lat)/Rmt.CellExtentInLatitude+1;
    lat_end = (Rmt.LatitudeLimits(2)-(lat-10))/Rmt.CellExtentInLatitude;
    lon_start = (lon-Rmt.LongitudeLimits(1))/Rmt.CellExtentInLongitude+1;
    lon_end = (lon+10-Rmt.LongitudeLimits(1))/Rmt.CellExtentInLongitude;
    mts = MT(lat_start:lat_end, lon_start:lon_end);
    mtmax = max(mts,[],'all');  
    if mtmax == 0  
        continue
    end   
    [m,n] = size(mts);
    re_mt = uint16(zeros(40000));
    for ii=1:m
        for jj=1:n
            re_mt((ii-1)*40+1:40*ii, (jj-1)*40+1:jj*40) = uint16(zeros(40)) + mts(ii,jj);
        end
    end    
    [LOSS, Rloss] = geotiffread(fn);
    LOSS(re_mt==0) = 0;     
    lat = Rloss.LatitudeLimits(2); %left corner lat
    lon = Rloss.LongitudeLimits(1); 
    [lons,lats] = meshgrid(lon+10/40000/2:10/40000:lon+10-10/40000/2,lat-10/40000/2:-10/40000:lat-10+10/40000/2);
    loss = [lats(LOSS>0),lons(LOSS>0),single(LOSS(LOSS>0))];  
    [mm,~] = size(loss);
    II = randperm(mm,ceil(mm/1000));   %¡¡1% random number
    lat_lon_yr_temp = [];
    lat_lon_yr_temp = loss(II,:);
    lat_lon_yr_save = [lat_lon_yr_save;lat_lon_yr_temp];
    t = toc;
    disp(['pixel ', num2str(pp),' costs ', num2str(t/60), ' min'])
end   
loss_yr = sortrows(lat_lon_yr_save,3); 
loss_all_save = [];
[nn,~]=size(lat_lon_yr_save);
JJ = randperm(nn,5000);   
loss_all_temp = lat_lon_yr_save(JJ,:);
loss_all_save = [loss_all_save;loss_all_temp];
loss_all_yr = sortrows(loss_all_save,3); 
xlswrite('sample_random.xlsx',loss_all_yr);

%% summary after visual interpretation
data = xlsread('F:\onedrive_xinyue\SupplementaryData.xlsx');
numel(find(isnan(data(:,6))))
numel(find(data(:,6)==1))
no_gain=data(data(:,6)==0,:);
gain=data(data(:,6)==1,:);
xlswrite('sample_nogain.xls',no_gain);
xlswrite('sample_gain.xls',gain);

% read elevation of forest gain 
gain = xlsread('sample_gain.xls');
lat = gain(:,2); lon = gain(:,3); 
fs = dir('J:\xinyue\ASTER GDEM v3\*.tif');
for pp = 1:length(gain)
    tic;
    disp(['begin point ', num2str(pp)])
    lon_start = floor(lon(pp)/10)*10; 
    lat_start = ceil(lat(pp)/10)*10;
    lon_read =num2str(abs(lon_start),'%03d'); 
    lat_read =num2str(abs(lat_start),'%02d');
    if lon_start <0 && lat_start <0
        fn = [char(fs(1).folder),'\ASTER_GDEM_V003_',lat_read,'S_',lon_read,'W.tif'];
    elseif lon_start >0 && lat_start <0
        fn = [char(fs(1).folder),'\ASTER_GDEM_V003_',lat_read,'S_',lon_read,'E.tif'];
    elseif lon_start <0 && lat_start >0
        fn = [char(fs(1).folder),'\ASTER_GDEM_V003_',lat_read,'N_',lon_read,'W.tif'];  
    elseif lon_start >0 && lat_start >0
        fn = [char(fs(1).folder),'\ASTER_GDEM_V003_',lat_read,'N_',lon_read,'E.tif'];
    end    
    [DEM Rdem] = geotiffread(fn);
    row = ceil(4000*(lat_start-lat(pp)));
    col = ceil(4000*abs(lon_start-lon(pp)));
    gain(pp,7) = DEM(row,col); 
    t = toc;
    disp(['point ', num2str(pp),' costs ', num2str(t/60), ' min'])
end 
xlswrite('sample_gain_elevation.xls',gain);
% p = histogram(gain(:,7),'Normalization','probability'); 
% n = histogram(gain(:,7));
