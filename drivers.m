[DRIVER,Rdriver] = geotiffread('F:\onedrive_xinyue\data\Curtis_foest divers\Goode_FinalClassification_19_05pcnt_prj.tif');
drivers = geotiffinfo('F:\onedrive_xinyue\data\Curtis_foest divers\Goode_FinalClassification_19_05pcnt_prj.tif');
%% driver for each country 
[MT,Rmt] = geotiffread('J:\xinyue\MOUNTAIN\country_mt2_re.tif'); % value by country
%[MT,Rmt] = geotiffread('J:\xinyue\MOUNTAIN\hotspot\revise_landgrid\country_rsr_re.tif'); % RSR_threatened hotspots by country 
%[MT,Rmt] = geotiffread('J:\xinyue\MOUNTAIN\hotspot\revise_landgrid\PAs\country_pa_rsrthr_re.tif'); % RSR(PAs) by country
lat_lon_yr_ct = [];
fs = dir('J:\xinyue\HANSEN\v1.6\*.tif'); 
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
    %clip global mt to 10*10 degree (1000*1000)
    lat_start = (Rmt.LatitudeLimits(2)-lat)/Rmt.CellExtentInLatitude+1;
    lat_end = (Rmt.LatitudeLimits(2)-(lat-10))/Rmt.CellExtentInLatitude;
    lon_start = (lon-Rmt.LongitudeLimits(1))/Rmt.CellExtentInLongitude+1;
    lon_end = (lon+10-Rmt.LongitudeLimits(1))/Rmt.CellExtentInLongitude;
    mts = MT(lat_start:lat_end, lon_start:lon_end);
    mtmax = max(mts,[],'all');  % the largest no. of mt
    if mtmax == 0  % if there is no mt in the region
        continue
    end   
    %resample mt to 30m (40000*40000)
    [m,n] = size(mts);
    re_mt = uint8(zeros(40000));
    for ii = 1:m
        for jj = 1:n
            re_mt((ii-1)*40+1:40*ii, (jj-1)*40+1:jj*40) = uint8(zeros(40)) + mts(ii,jj);
        end
    end    
    [LOSS, Rloss] = geotiffread(fn);
    LOSS(re_mt==0) = 0;    
    [lons,lats] = meshgrid(lon+10/40000/2:10/40000:lon+10-10/40000/2,lat-10/40000/2:-10/40000:lat-10+10/40000/2);
    loss = [lats(LOSS>0),lons(LOSS>0),single(LOSS(LOSS>0)),single(re_mt(LOSS>0))];  
    lat_lon_yr_ct = [lat_lon_yr_ct;loss];
    t = toc;
    disp(['pixel ', num2str(pp),' costs ', num2str(t/60), ' min'])
end   
fn = ['lat_lon_yr_ct.mat'];
save(fn,'lat_lon_yr_ct')

ct_drivers = nan(141,5);
for mt=1:141
    if  any(lat_lon_yr_ct(:,4)==mt)==0
        continue
    end
    bycountry = lat_lon_yr_ct(lat_lon_yr_ct(:,4)==mt,:);
    lat =  bycountry(:,1);
    lon =  bycountry(:,2);
    for pp = 1:length(bycountry)
        row = ceil((89.4992225958-lat(pp))/drivers.PixelScale(1));
        col = ceil((abs(179.999984049+lon(pp))/drivers.PixelScale(1)));
        bycountry(pp,5) = DRIVER(row,col);
    end
    tbl = tabulate(bycountry(:,5)); %calculate the percentage
    for driver = 1:5
        if  any(tbl(:,1)==driver)==0
            continue
        end
        ct_drivers(mt,driver)=tbl(tbl(:,1)==driver,3);
    end
end

% summary for all mountains
clear;
load('lat_lon_yr_ct.mat');
lat_lon_driver = lat_lon_yr_ct(:,1:2);
lat =  lat_lon_yr_ct(:,1);
lon = lat_lon_yr_ct(:,2);
for pp = 1:length(lat_lon_driver)
    row = ceil((89.4992225958-lat(pp))/drivers.PixelScale(1));
    col = ceil((abs(179.999984049+lon(pp))/drivers.PixelScale(1)));
    lat_lon_driver(pp,3) = DRIVER(row,col);
end
uni_drivers = unique(lat_lon_driver(:,3));
tbl = tabulate(lat_lon_driver(:,3)); 
tbl(8,3) = 100-tbl(3,3)-tbl(4,3)-tbl(5,3)-tbl(6,3)-tbl(7,3);

% summary for tropical, temperate, boreal regions
load('lat_lon_driver.mat');
lat =  lat_lon_driver(:,1);
ii = find(lat>=50);%boreal
lat_lon_boreal = (lat_lon_driver(ii,:));
tbl_boreal = tabulate(lat_lon_boreal(:,3)); 
tbl_boreal(8,3) = 100-tbl_boreal(3,3)-tbl_boreal(4,3)-tbl_boreal(5,3)-tbl_boreal(6,3)-tbl_boreal(7,3);
jj = find(lat<24 & lat>-24);%tropical
lat_lon_tropical = (lat_lon_driver(jj,:));
tbl_tropical = tabulate(lat_lon_tropical(:,3));
tbl_tropical(8,3) = 100-tbl_tropical(3,3)-tbl_tropical(4,3)-tbl_tropical(5,3)-tbl_tropical(6,3)-tbl_tropical(7,3);
kk = find(lat<50 & lat >=24 | lat <=-24);%temperate
lat_lon_temperate = (lat_lon_driver(kk,:));
tbl_temperate = tabulate(lat_lon_temperate(:,3)); 
tbl_temperate(8,3) = 100-tbl_temperate(3,3)-tbl_temperate(4,3)-tbl_temperate(5,3)-tbl_temperate(6,3)-tbl_temperate(7,3);

%% bar plot (figure 3a)
C = [tbl(3,3),tbl_tropical(3,3),tbl_temperate(3,3),tbl_boreal(3,3)];     
S = [tbl(4,3),tbl_tropical(4,3),tbl_temperate(4,3),tbl_boreal(4,3)];
F = [tbl(5,3),tbl_tropical(5,3),tbl_temperate(5,3),tbl_boreal(5,3)];
W = [tbl(6,3),tbl_tropical(6,3),tbl_temperate(6,3),tbl_boreal(6,3)];
U = [tbl(7,3),tbl_tropical(7,3),tbl_temperate(7,3),tbl_boreal(7,3)];
O = [tbl(8,3),tbl_tropical(8,3),tbl_temperate(8,3),tbl_boreal(8,3)];
y = [C;S;F;W;U;O];
b = bar(y,1);
set(b,'edgecolor','none');
h = legend('Global','Topical', 'Temperate', 'Boreal');
set(h,'Box','off');
xticklabels({'C','S','F','W','U','O'});
ylabel('Percentage of various drivers (%)')
%% bar plot (figure 3b)
data = xlsread('F:\onedrive_xinyue\mountainforests_backup\results\matlab\revise\drivers\summary drivers.xlsx');
driver_rsr = data(3:8,7);driver_rsr_pa = data(3:8,11);
C = [driver_rsr(1,1),driver_rsr_pa(1,1)];     
S = [driver_rsr(2,1),driver_rsr_pa(2,1)];
F = [driver_rsr(3,1),driver_rsr_pa(3,1)];
W = [driver_rsr(4,1),driver_rsr_pa(4,1)];
U = [driver_rsr(5,1),driver_rsr_pa(5,1)];
O = [driver_rsr(6,1),driver_rsr_pa(6,1)];
y = [C;S;F;W;U;O];
b1 = bar(y,1)
set(b1,'edgecolor','none');
map = [200 67 64;96 57 18]/255;
for ki = 1:2
    set(b1(ki), 'Facecolor',map(ki,:))  
end
h1 = legend( 'RSR', 'RSR (PAs)');
set(h1,'Box','off');
xticklabels({'C','S','F','W','U','O'});
ylabel('Percentage of various drivers (%)')
