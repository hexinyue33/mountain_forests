clear;clc;
[MT,Rmt] = geotiffread('J:\xinyue\MOUNTAIN\MT_re.tif');
fs1=dir('J:\xinyue\HANSEN\v1.6\*.tif'); 
mt_fml_tile = nan(1048,14,504);
lat_fml = [80.000496;44.9995039365;9.99851187302;-25.0024801905];
lon_fml = [-180.000496;-135.000496;-90.000496;-45.000496;-0.000496000000595;44.999504;89.999504;134.999504];
for pp = 1:length(fs1)
    tic;
    disp(['begin tile ', num2str(pp)])
    fn1 = [char(fs1(pp).folder),'/',char(fs1(pp).name)];
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
     % 1000*1000
    lat_start = (Rmt.LatitudeLimits(2)-lat)/Rmt.CellExtentInLatitude+1;
    lat_end = (Rmt.LatitudeLimits(2)-(lat-10))/Rmt.CellExtentInLatitude;
    lon_start = (lon-Rmt.LongitudeLimits(1))/Rmt.CellExtentInLongitude+1;
    lon_end = (lon+10-Rmt.LongitudeLimits(1))/Rmt.CellExtentInLongitude;
    mts = MT(lat_start:lat_end, lon_start:lon_end);  
    mtmax = max(mts,[],'all');  % the largest no. of mt  
    if mtmax == 0  % if there is no mt in the region
        continue
    end    
   [LOSS, Rloss] = geotiffread(fn1);
   LOSS_info = geotiffinfo(fn1);
    % resample mt data
    [m,n] = size(mts);
    re_mt = uint16(zeros(40000));
    for ii=1:m
        for jj=1:n
            re_mt((ii-1)*40+1:40*ii, (jj-1)*40+1:jj*40) = uint16(zeros(40)) + mts(ii,jj);
        end
    end     
    mts = unique(mts(mts>0));   %%display the same value once
    disp(['mt total = ', num2str(length(mts))])
    loss_mt_fml = nan(1048,14); 
    for kk = 1:1:length(mts)    %%from 1 to 1048
        ROW=[]; COL=[];  new_class=[];group=[];
        mt = mts(kk); 
        disp(['mt range ', num2str(mt)])
        [a,b] = find(re_mt==mt & LOSS==15);%% loss in 2015
        if length(a) == 0
            continue 
        end
        [LAT,LON] = pix2latlon(LOSS_info.RefMatrix,a,b); %center
        for aa = 1:length(LAT)
%             disp(['begin ', num2str(aa),'/',num2str(length(LAT))])
            if LAT(aa) >lat_fml(2,1) && LON(aa) < lon_fml(2,1)
               ROW(aa) = 1; COL(aa) = 1;
            elseif  LAT(aa) <lat_fml(2,1) && LAT(aa) >lat_fml(3,1) && LON(aa) < lon_fml(2,1)
               ROW(aa) = 2; COL(aa) = 1;        
            elseif  LAT(aa) <lat_fml(3,1) && LAT(aa) >lat_fml(4,1) && LON(aa) < lon_fml(2,1)
               ROW(aa) = 3; COL(aa) = 1;          
            elseif  LAT(aa) <lat_fml(4,1) && LON(aa) <lon_fml(2,1)
               ROW(aa) = 4; COL(aa) = 1;     
               
            elseif LAT(aa) >lat_fml(2,1) && LON(aa) > lon_fml(2,1) && LON(aa)<lon_fml(3,1)
               ROW(aa) = 1; COL(aa) = 2;
            elseif  LAT(aa) <lat_fml(2,1) && LAT(aa) >lat_fml(3,1) && LON(aa) > lon_fml(2,1) && LON(aa)<lon_fml(3,1)
               ROW(aa) = 2; COL(aa) = 2;        
            elseif LAT(aa) <lat_fml(3,1) && LAT(aa) >lat_fml(4,1) && LON(aa) > lon_fml(2,1) && LON(aa)<lon_fml(3,1)
               ROW(aa) = 3; COL(aa) = 2;          
            elseif  LAT(aa) <lat_fml(4,1) && LON(aa) > lon_fml(2,1) && LON(aa)<lon_fml(3,1)
               ROW(aa) = 4; COL(aa) = 2;   
               
             elseif LAT(aa) >lat_fml(2,1) && LON(aa) > lon_fml(3,1) && LON(aa)<lon_fml(4,1)
               ROW(aa) = 1; COL(aa) = 3;
            elseif  LAT(aa) <lat_fml(2,1) && LAT(aa) >lat_fml(3,1) && LON(aa) > lon_fml(3,1) && LON(aa)<lon_fml(4,1)
               ROW(aa) = 2; COL(aa) = 3;        
            elseif LAT(aa) <lat_fml(3,1) && LAT(aa) >lat_fml(4,1) && LON(aa) > lon_fml(3,1) && LON(aa)<lon_fml(4,1)
               ROW(aa) = 3; COL(aa) = 3;          
            elseif  LAT(aa) <lat_fml(4,1) && LON(aa) > lon_fml(3,1) && LON(aa)<lon_fml(4,1)
               ROW(aa) = 4; COL(aa) = 3;   
               
            elseif LAT(aa) >lat_fml(2,1) && LON(aa) > lon_fml(4,1) && LON(aa)<lon_fml(5,1)
               ROW(aa) = 1; COL(aa) = 4;
            elseif  LAT(aa) <lat_fml(2,1) && LAT(aa) >lat_fml(3,1) && LON(aa) > lon_fml(4,1) && LON(aa)<lon_fml(5,1)
               ROW(aa) = 2; COL(aa) = 4;        
            elseif LAT(aa) <lat_fml(3,1) && LAT(aa) >lat_fml(4,1) && LON(aa) > lon_fml(4,1) && LON(aa)<lon_fml(5,1)
               ROW(aa) = 3; COL(aa) = 4;          
            elseif  LAT(aa) <lat_fml(4,1) && LON(aa) > lon_fml(4,1) && LON(aa)<lon_fml(5,1)
               ROW(aa) = 4; COL(aa) = 4;   
               
             elseif LAT(aa) >lat_fml(2,1) && LON(aa) > lon_fml(5,1) && LON(aa)<lon_fml(6,1)
               ROW(aa) = 1; COL(aa) = 5;
            elseif  LAT(aa) <lat_fml(2,1) && LAT(aa) >lat_fml(3,1) && LON(aa) > lon_fml(5,1) && LON(aa)<lon_fml(6,1)
               ROW(aa) = 2; COL(aa) = 5;        
            elseif LAT(aa) <lat_fml(3,1) && LAT(aa) >lat_fml(4,1) && LON(aa) > lon_fml(5,1) && LON(aa)<lon_fml(6,1)
               ROW(aa) = 3; COL(aa) = 5;          
            elseif  LAT(aa) <lat_fml(4,1) && LON(aa) > lon_fml(5,1) && LON(aa)<lon_fml(6,1)
               ROW(aa) = 4; COL(aa) = 5;         
               
             elseif LAT(aa) >lat_fml(2,1) && LON(aa) > lon_fml(6,1) && LON(aa)<lon_fml(7,1)
               ROW(aa) = 1; COL(aa) = 6;
            elseif  LAT(aa) <lat_fml(2,1) && LAT(aa) >lat_fml(3,1) && LON(aa) > lon_fml(6,1) && LON(aa)<lon_fml(7,1)
               ROW(aa) = 2; COL(aa) = 6;        
            elseif LAT(aa) <lat_fml(3,1) && LAT(aa) >lat_fml(4,1) && LON(aa) > lon_fml(6,1) && LON(aa)<lon_fml(7,1)
               ROW(aa) = 3; COL(aa) = 6;          
            elseif  LAT(aa) <lat_fml(4,1) && LON(aa) > lon_fml(6,1) && LON(aa)<lon_fml(7,1)
               ROW(aa) = 4; COL(aa) = 6;   
               
            elseif LAT(aa) >lat_fml(2,1) && LON(aa) > lon_fml(7,1) && LON(aa)<lon_fml(8,1)
               ROW(aa) = 1; COL(aa) = 7;
            elseif  LAT(aa) <lat_fml(2,1) && LAT(aa) >lat_fml(3,1) && LON(aa) > lon_fml(7,1) && LON(aa)<lon_fml(8,1)
               ROW(aa) = 2; COL(aa) = 7;        
            elseif LAT(aa) <lat_fml(3,1) && LAT(aa) >lat_fml(4,1) && LON(aa) > lon_fml(7,1) && LON(aa)<lon_fml(8,1)
               ROW(aa) = 3; COL(aa) = 7;          
            elseif  LAT(aa) <lat_fml(4,1) && LON(aa) > lon_fml(7,1) && LON(aa)<lon_fml(8,1)
               ROW(aa) = 4; COL(aa) = 7;   
               
             elseif LAT(aa) >lat_fml(2,1) && LON(aa) > lon_fml(8,1) 
               ROW(aa) = 1; COL(aa) = 8;
            elseif  LAT(aa) <lat_fml(2,1) && LAT(aa) >lat_fml(3,1) && LON(aa) > lon_fml(8,1) 
               ROW(aa) = 2; COL(aa) = 8;        
            elseif LAT(aa) <lat_fml(3,1) && LAT(aa) >lat_fml(4,1) && LON(aa) > lon_fml(8,1) 
               ROW(aa) = 3; COL(aa) = 8;          
            elseif  LAT(aa) <lat_fml(4,1) && LON(aa) > lon_fml(8,1) 
               ROW(aa) = 4; COL(aa) = 8;      
            end
        end  
            group(:,1) = ROW';
            group(:,2) = COL';
            uni_group = unique(group,'rows');
            s_unigroup = size(uni_group);
        for bb = 1:s_unigroup(1,1)
            cc = find(ROW==uni_group(bb,1) & COL == uni_group(bb,2));
            class =[];
            for dd = 1:length(cc)        
                fn = ['J:\xinyue\GlobalForestManagement\tiles v3.2\tiles v3.2\FML_v3-2_with-colorbar_',num2str(ROW(cc(dd))),'_',num2str(COL(cc(dd))),'.tif'];
                if dd == 1
                [FML,Rfml] = geotiffread(fn);
                lat_fml_start = lat_fml(ROW(cc(dd)));
                lon_fml_start = lon_fml(COL(cc(dd)));
                row_fml = ceil((lat_fml_start-LAT(cc(dd)))/((lat_fml(1,1)-lat_fml(2,1))/35281));
                col_fml = ceil(abs(lon_fml_start-LON(cc(dd)))/((lon_fml(2,1)-lon_fml(1,1))/45360));
                class(dd,1) = LAT(cc(dd));
                class(dd,2) = LON(cc(dd));
                class(dd,3) = FML(row_fml,col_fml);
                else 
                lat_fml_start = lat_fml(ROW(cc(dd)));
                lon_fml_start = lon_fml(COL(cc(dd)));
                row_fml = ceil((lat_fml_start-LAT(cc(dd)))/((lat_fml(1,1)-lat_fml(2,1))/35281));
                col_fml = ceil(abs(lon_fml_start-LON(cc(dd)))/((lon_fml(2,1)-lon_fml(1,1))/45360));
                class(dd,1) = LAT(cc(dd));
                class(dd,2) = LON(cc(dd));
                class(dd,3) = FML(row_fml,col_fml);
                end
            end
            new_class(length(new_class)*(bb-1)+1:length(new_class)*(bb-1)+length(cc),:) = class(:,:);
        end       
        tab = tabulate(new_class(:,3));     
        tab1= tab(:,1);
        if length(find(tab1==11))>0
            loss_mt_fml(mt,1:2) = tab(find(tab1==11),2:3);
        end
        if length(find(tab1==20))>0
            loss_mt_fml(mt,3:4) = tab(find(tab1==20),2:3);
        end
        if length(find(tab1==31))>0
            loss_mt_fml(mt,5:6) = tab(find(tab1==31),2:3);
        end
        if length(find(tab1==32))>0
            loss_mt_fml(mt,7:8) = tab(find(tab1==32),2:3);
        end
        if length(find(tab1==40))>0
           loss_mt_fml(mt,9:10) = tab(find(tab1==40),2:3);
        end
        if length(find(tab1==53))>0
           loss_mt_fml(mt,11:12) = tab(find(tab1==53),2:3);
        end
        if length(find(tab1==128))>0
           loss_mt_fml(mt,13:14) = tab(find(tab1==128),2:3);
        end
    end
    mt_fml_tile(:,:,pp)= loss_mt_fml;
    t = toc;
    disp(['tile ', num2str(pp),' costs ', num2str(t/60), ' min'])
end
%% stat.
sumup = nansum(mt_fml_tile,3);
stats = zeros(7,2);
stats(1,1) = sum(sumup(:,1));
stats(2,1) = sum(sumup(:,3));
stats(3,1) = sum(sumup(:,5));
stats(4,1) = sum(sumup(:,7));
stats(5,1) = sum(sumup(:,9));
stats(6,1) = sum(sumup(:,11));
stats(7,1) = sum(sumup(:,13));
stats(1,2) = stats(1,1)/sum(stats(:,1))*100;
stats(2,2) = stats(2,1)/sum(stats(:,1))*100;
stats(3,2) = stats(3,1)/sum(stats(:,1))*100;
stats(4,2) = stats(4,1)/sum(stats(:,1))*100;
stats(5,2) = stats(5,1)/sum(stats(:,1))*100;
stats(6,2) = stats(6,1)/sum(stats(:,1))*100;
stats(7,2) = stats(7,1)/sum(stats(:,1))*100;

% seven to three classes
data = xlsread('J:\xinyue\1_mobaxterm\revise\ForestManagement\NEW\loss in 2015\STATS_FML.xlsx');
globe(1,1) = sum(data(1:2,6));
globe(2,1) = sum(data(3:6,6));
globe(3,1) = sum(data(7,6));
tro(1,1) = sum(data(1:2,9));
tro(2,1) = sum(data(3:6,9));
tro(3,1) = sum(data(7,9));
temp(1,1) = sum(data(1:2,12));
temp(2,1) = sum(data(3:6,12));
temp(3,1) = sum(data(7,12));
bor(1,1) = sum(data(1:2,15));
bor(2,1) = sum(data(3:6,15));
bor(3,1) = sum(data(7,15));
combine = [globe,tro,temp,bor];
b=bar(combine,1)
set(b,'edgecolor','none');
h = legend( 'Global', 'Tropical','Temperate','Boreal');
set(h,'Box','off');
ylim([0 100])
xticklabels({'Natural regenerating forests','Plantations','No data'});
ylabel('Percentage of forest management in 2015 (%)')
