clear;clc;
fs=dir('J:\xinyue\HANSEN\v1.6\*.tif'); 
[MT,Rmt] = geotiffread('J:\xinyue\MOUNTAIN\MT_re.tif');
loss_mt = nan(360,720,18);
for pp =1:length(fs)
    tic;
    disp(['begin pixel ', num2str(pp)])
    fn1 = [char(fs(pp).folder),'\',char(fs(pp).name)];
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
    lat_start = (Rmt.LatitudeLimits(2)-lat)/Rmt.CellExtentInLatitude+1;
    lat_end = (Rmt.LatitudeLimits(2)-(lat-10))/Rmt.CellExtentInLatitude;
    lon_start = (lon-Rmt.LongitudeLimits(1))/Rmt.CellExtentInLongitude+1;
    lon_end = (lon+10-Rmt.LongitudeLimits(1))/Rmt.CellExtentInLongitude;
    mts = MT(lat_start:lat_end, lon_start:lon_end);  %1000*1000 
    mtmax = max(mts,[],'all'); 
    if mtmax == 0  % if there is no mt in the region
        continue
    end    
    % resample mt data, to 40000*40000
    [m,n] = size(mts);
    re_mt = uint8(zeros(40000));
    for ii=1:m
        for jj=1:n
            re_mt((ii-1)*40+1:40*ii, (jj-1)*40+1:jj*40) = uint8(zeros(40)) + mts(ii,jj);
        end
    end   
    [LOSS, Rloss] = geotiffread(fn1);   
    if sum(LOSS>0,'all') == 0
       continue
    end      
    % 0.5 degree(20*20) in 10 degree
    loss = nan(20,20,18); 
    for ii=1:20 
        disp(['begin loop ', num2str(ii)])
        for jj=1:20 
            grid_loss = LOSS((ii-1)*2000+1:ii*2000,(jj-1)*2000+1:jj*2000);  
            grid_mt = re_mt((ii-1)*2000+1:(ii-1)*2000+2000,(jj-1)*2000+1:jj*2000);       
            gridmtmax = max(grid_mt,[],'all');  
            if gridmtmax == 0  % if there is no mt in this grid
                 continue
            end  
            for yr = 1:18       
                  loss(ii,jj,yr) = sum(grid_mt>0 & grid_loss==yr,'all')*cosd(lat)*900/10000;   %unit: ha 
            end
        end  
    end
    row = (90-lat)*2+1;
    col = (180+lon)*2+1;
    loss_mt(row:row+19,col:col+19,:,:) = loss; 
    t = toc;
    disp(['pixel ', num2str(pp),' costs ', num2str(t/60), ' min'])
end
fn = ['loss_mt.mat'];
save(fn,'loss_mt')
%% total loss and loss trend in each 0.5 degree grid
clear;
load('loss_mt.mat');
total_loss = zeros(360,720); trend_loss = nan(360,720); 
for aa = 1:360 
    for bb = 1:720 
        grid = squeeze(loss_mt(aa,bb,:,:));
        loss = grid(:,1);
        ele = grid(:,2);   
        total_loss(aa,bb)=nansum(grid(:,1));
          % Theil-sen estimator
         N=size(loss,1); 
         x=(1:N)'; 
         bls = regress(loss,[ones(N,1) x]);
         Comb = combnk(1:N,2);
         theil=diff(loss(Comb),1,2)./diff(x(Comb),1,2);
         b_loss=median(theil);% b is the slope, i.e., changes per step.
         trend_loss(aa,bb)=b_loss;     
    end
end
fn1 = ['total_loss.mat'];
save(fn1,'total_loss')
fn2 = ['trend_loss.mat'];
save(fn2,'trend_loss')
%% spatial pattern (figure 2)
load('total_loss.mat');load('trend_loss.mat')
MT = geotiffread('F:\onedrive_xinyue\mountainforests_backup\datasets\MT_revalues_0.5\mt_point5.tif');
mt = shaperead('F:\onedrive_xinyue\mountainforests_backup\datasets\GMBA mountain inventory V1.2(entire world)_0.01\GMBAdissolve.shp','UseGeoCoords',true);
load('F:\onedrive_xinyue\phd_study\code\matlab practice\colorbar_matlab\MPL_RdYlGn_r.mat');
MPL_RdYlGn_r(64:65,:) = repmat([220,220,220]/255,[2,1]);
for kk = 1:3
    MPL_RdYlGn_r(65:72,kk) = MPL_RdYlGn_r(65,kk):(MPL_RdYlGn_r(72,kk)-MPL_RdYlGn_r(65,kk))/7:MPL_RdYlGn_r(72,kk);
    MPL_RdYlGn_r(57:64,kk) = MPL_RdYlGn_r(57,kk):(MPL_RdYlGn_r(64,kk)-MPL_RdYlGn_r(57,kk))/7:MPL_RdYlGn_r(64,kk);
end
load coast
long(lat<-56)=nan;
lat(lat<-56)=nan;
[LON,LAT]=meshgrid(-180:0.5:179.5,90:-0.5:-89.5);
trend_loss(isnan(trend_loss) & MT>0) = 0; %b
total_loss(MT==0)=nan;%a

fig = figure('units','centimeters','position',[5,1.5,18,19]);
ax2 = axes;
set(ax2,'position',[0 0 1 0.5]);
axesm('MapProjection','robinson','MapLatLimit',[-90 90],'MapLonLimit',[-180 180],'frame','off','grid','on','MLineLocation',20,...
    'MeridianLabel','off','MLabelParallel','south', 'PLineLocation',10,'ParallelLabel','off','PLabelMeridian','east','FontSize',20)
hold on;
h2 = framem;
set(h2,'LineWidth',0.5)
pcolorm(LAT,LON, trend_loss)
caxis([-60 60])
colormap(ax2,MPL_RdYlGn_r)
set(gca,'Visible','off')
plotm([mt.Lat],[mt.Lon],'-','LineWidth',0.1,'Color',[128 128 128]/255)
plotm(lat,long,'-k','LineWidth',0.5)
h = colorbar('location','northoutside');
set(h,'position',[0.38 0.08 0.4 0.015],'xtick',-60:30:60,'xticklabel',{'<-60','-30','0','30','>60'},'FontSize',9);
h1 = text(0,-1.25,'Accerelation in mountain forest loss (ha yr^-^2)','FontSize',10);
set(h1,'Position',[-0.65,-0.85])
text(-2.5,1.3,'b','FontSize',12,'FontWeight','bold')

ax1 = axes;
set(ax1,'position',[0 0.5 1 0.5]);
axesm('MapProjection','robinson','MapLatLimit',[-90 90],'MapLonLimit',[-180 180],'frame','off','grid','on','MLineLocation',20,...
    'MeridianLabel','off','MLabelParallel','south', 'PLineLocation',10,'ParallelLabel','off','PLabelMeridian','east','FontSize',20)
hold on;
h2 = framem;
set(h2,'LineWidth',0.5)
pcolorm(LAT,LON, total_loss)
caxis([0 50000])
colormap(ax1,MPL_RdYlGn_r(65:end,:))
set(gca,'Visible','off')
plotm([mt.Lat],[mt.Lon],'-','LineWidth',0.1,'Color',[128 128 128]/255)
plotm(lat,long,'-k','LineWidth',0.5)
h = colorbar('location','southoutside');
set(h,'position',[0.38 0.58 0.4 0.015],'xtick',0:10000:50000,'xticklabel',{'0','1','2','3','4','>5'},'FontSize',9);
h1 = text(0,0.465,'Total mountain forest loss (10^4 ha)','FontSize',10);
set(h1,'Position',[-0.45,-0.87])
text(-2.5,1.3,'a','FontSize',12,'FontWeight','bold')

set(gcf, 'renderer', 'painters');
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig.Position(3:4),'PaperPosition',[0,0,1,1]) 
