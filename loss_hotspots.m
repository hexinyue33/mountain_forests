%% Mountain forest loss within the RSR(threatened) biodiversity hotspots
clear;clc;
[MT,Rmt] = geotiffread('J:\xinyue\MOUNTAIN\MT_re.tif');
[RSR,Rrsr]=geotiffread('J:\xinyue\BIODIVERSITY\RSR_clip.tif'); 
fs1=dir('J:\xinyue\HANSEN\v1.6\*.tif'); 
for pp =1:length(fs1); 
    tic;
    disp(['begin pixel ', num2str(pp)])
    fn1 = [char(fs1(pp).folder),'/',char(fs1(pp).name)];
    if fn1(end-9) == 'N' 
        lat = str2double(fn1(end-11:end-10));
    elseif fn1(end-9) == 'S'
        lat = -str2double(fn1(end-11:end-10));
    end    
    if fn1(end-4) == 'E'
        lon = str2double(fn1(end-7:end-5));
    elseif fn1(end-4) == 'W'
        lon = -str2double(fn1(end-7:end-5));
    end         
    fn = ['loss_rsrthr',fn1(end-11:end-3),'mat'];       
    if exist(fn,'file') == 2
        continue
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
    rsr_new = RSR(lat_start:lat_end, lon_start:lon_end);
    rsrmax = max(rsr_new,[],'all');    
    if rsrmax < 0.00012 % Outside of hotspots
        continue
    end 
    [LOSS, Rloss] = geotiffread(fn1); 
    [m,n] = size(mts);
    re_mt = uint16(zeros(40000));
    for ii=1:m
        for jj=1:n
            re_mt((ii-1)*40+1:40*ii, (jj-1)*40+1:jj*40) = uint16(zeros(40)) + mts(ii,jj);
        end
    end   
    re_rsr=single(zeros(40000));
     for ii=1:m
        for jj=1:n
            re_rsr((ii-1)*40+1:40*ii, (jj-1)*40+1:jj*40) = single(zeros(40)) + rsr_new(ii,jj);
        end
    end  
    rersrmax = max(re_rsr,[],'all');
    rersrmin = min(re_rsr(re_rsr>0)); 
    mts = unique(mts(mts>0));   %%display the same value once     
    disp(['mt range', num2str(length(mts))])
    loss_hotspot_mt = nan(1048,18);  
    for kk = 1:1:length(mts)   
        mt = mts(kk); 
        for yr = 1:18           
        loss_hotspot_mt(mt,yr) = sum(re_mt==mt & LOSS==yr & re_rsr>=0.00012,'all')*cosd(lat)*900/10000; 
        end         
    end        
    save(fn,'loss_hotspot_mt')
    t = toc;
    disp(['pixel ', num2str(pp),' costs ', num2str(t/60), ' min'])
end
%% stat.
clear;
fs=dir('J:\xinyue\mobaxterm\revise\LANDGRID\loss_rsrthr\*.mat'); fs=struct2cell(fs); fs=fs(1,:)';
for yr=1:18;
    for pp = 1:length(fs); 
        display(pp)
        fn = char(fs(pp,1)); 
        load(fn);
        sum_loss_pixel(pp,yr)=nansum(loss_hotspot_mt(:,yr));
    end
end
sum_loss=sum(sum_loss_pixel); 
sum(sum_loss)

% Theil-sen estimates 
loss=(sum_loss)';
N = size(loss,1); % [col, 1]
x = [2001:2018]'; % [col, 1]
bls = regress(loss,[ones(N,1) x]);
Comb = combnk(1:N,2);
theil = diff(loss(Comb),1,2)./diff(x(Comb),1,2);
b = median(theil);% b is the slope, i.e., changes per step. 

% MK significant test 
y = sum_loss;
n = length( y );
dt=1;
% calculate statistic
s = 0;
for k = 1:n-1,
    for j = k+1:n,
        s = s + sign( y(j) - y(k) );
    end;
end;
% variance ( assuming no tied groups )
v = ( n * ( n - 1 ) * ( 2 * n + 5 ) ) / 18;
% test statistic
if s == 0,
z = 0;
elseif s > 0,
z = ( s - 1 ) / sqrt( v );
else
z = ( s + 1 ) / sqrt( v );
end; 
% should calculate Normal value here
nor = 1.96;
% results
disp( [ ' n = ' num2str( n ) ] );
disp( [ ' Mean Value = ' num2str( mean( y ) ) ] );
disp( [ ' Z statistic = ' num2str( z ) ] );
if abs( z ) < nor,
disp( 'No significant trend' );
z = 0;
elseif z > 0,
disp( 'Upward trend detected' );
else
disp( 'Downward trend detected' );
end
