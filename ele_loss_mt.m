%% forest loss in elevation bin (50 m)
clear;clc;
fs1 = dir('J:\xinyue\ASTER GDEM v3\*.tif'); 
fs = dir('J:\xinyue\HANSEN\v1.6\*.tif'); 
[MT,Rmt] = geotiffread('/nfs/a336/eexhe/matlab/MT_revalue_0.01/MT_re.tif');
for pp =1:length(fs) 
    tic;
    disp(['begin pixel ', num2str(pp)])
    fn1 = [char(fs(pp).folder),'\',char(fs(pp).name)];
    fn2 = [char(fs1(1).folder),'\ASTER_GDEM_V003_',fn1(end-11:end-4),'.tif']; 
    % pixels lat lon
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
    fn = ['bin_loss_',fn1(end-11:end-3),'mat'];       
    if exist(fn,'file') == 2
       continue
    end   
    lat_start = (Rmt.LatitudeLimits(2)-lat)/Rmt.CellExtentInLatitude+1;
    lat_end = (Rmt.LatitudeLimits(2)-(lat-10))/Rmt.CellExtentInLatitude;
    lon_start = (lon-Rmt.LongitudeLimits(1))/Rmt.CellExtentInLongitude+1;
    lon_end = (lon+10-Rmt.LongitudeLimits(1))/Rmt.CellExtentInLongitude;
    mts = MT(lat_start:lat_end, lon_start:lon_end);
    mtmax = max(mts,[],'all');  
    if mtmax == 0  % if there is no mt in the region
        continue
    end  
    [LOSS, Rloss] = geotiffread(fn1);   
    % resample mt data
    [m,n] = size(mts);
    re_mt = int16(zeros(40000));
    for ii=1:m
        for jj=1:n
            re_mt((ii-1)*40+1:40*ii, (jj-1)*40+1:jj*40) = int16(zeros(40)) + mts(ii,jj);
        end
    end   
    [DEM,Rdem] = geotiffread(fn2);
    loss_ele = nan(20,20,18,71);
    for ii = 1:20 
        disp(['begin loop ', num2str(ii)])
        for jj = 1:20 
            grid_loss = LOSS((ii-1)*2000+1:ii*2000,(jj-1)*2000+1:jj*2000);  
%             if sum(sum(grid_loss>0),'all') < 100  % ignore pixels with small loss
%                 continue
%             end    
            grid_mt = re_mt((ii-1)*2000+1:(ii-1)*2000+2000,(jj-1)*2000+1:jj*2000);  
            gridmtmax = max(grid_mt,[],'all');  
            gridmtmin = min(grid_mt(grid_mt>0)); 
            if gridmtmax == 0  % if there is no mt in this grid
                 continue
            end     
            grid_dem = DEM((ii-1)*2000+1:(ii-1)*2000+2000,(jj-1)*2000+1:jj*2000);         
            elemax = ceil(double(max(grid_dem(grid_dem>0)))/50)*50;
            elemin = ceil(double(min(grid_dem(grid_dem>0)))/50)*50;           
            for yr = 1:18       
                for ele = elemin:50:min(elemax,3500)
                    loss_ele(ii,jj,yr,ele/50) = sum(grid_dem>=ele-50 & grid_dem<ele & grid_loss==yr & grid_mt>0,'all')*cosd(lat)*900/10000;
                end                
                if elemax >= 3500
                    loss_ele(ii,jj,yr,71) = sum(grid_dem>=3500 & grid_loss==yr & grid_mt>0,'all')*cosd(lat)*900/10000;
                end                  
            end           
        end
    end
    save(fn,'loss_ele')
    t = toc;
    disp(['pixel ', num2str(pp),' costs ', num2str(t/60), ' min'])
end
%% stat.
clear;clc;
fs = dir('J:\xinyue\mobaxterm\revise\LANDGRID\bin_loss\*.mat'); fs = struct2cell(fs); fs = fs(1,:)';
for pp = 1:length(fs); 
    display(pp)
    fn = char(fs(pp,1)); 
    load(fn);  
    for yr = 1:18; 
        for ele = 1:71;
            loss_grid = loss_ele(:,:,yr,ele);
            sum_grid(ele,yr,pp) = nansum(nansum(loss_grid));
        end
    end   
end
loss = sum(sum_grid,3);
loss_elevation = sum(loss,2);
% ele = [0:50:3500];
% plot(ele,loss_elevation)
