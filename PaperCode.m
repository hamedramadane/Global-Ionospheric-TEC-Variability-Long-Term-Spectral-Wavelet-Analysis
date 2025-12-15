%% --------------------------Code for the Paper Plots---------------------%

clc;
clear all;
%% -------------------------Fifure One------------------------------------%

day23 = cdfread('gps_tec15min_igs_20230303_v01.cdf');
lat = -87.5:2.5:87.5;    
lon = -180:5:180;

for k = 1:96

    mp = day23{k,2};
    
end

mp1 = day23{89,2};
mp2 = day23{89,2};        % Choose k={1,10,17,25,33,41,49,57,65,73,81,89}                                             
pcolor(lon,lat,mp2) 
xlabel ('\bf Longitude','FontSize',25);
ylabel ('\bf Latitude','FontSize',25); 
title('\bf March, 03/2023  22:00 UT','FontSize',25);
fontsize(gca, 15,'points') 
set(gca,'FontWeight','bold','FontSize',25);
c=colorbar('Fontsize',15);
c.Title.String = 'TECU';
grid on
shading interp 
load coastlines
geoshow(coastlat,coastlon,'color','k')

%----------------------Adding Terminators lines---------------------------%
TEC_data = mp2; 
t = datetime(2023,3,3,22,0,0,'TimeZone','UTC');  

plot_TEC_with_daynight2(lon, lat, TEC_data, t);

fontsize(gca, 25,'points') 
xlabel ('\bf Longitude','FontSize',22);
ylabel ('\bf Latitude','FontSize',20); 
title('\bf   22:00 UT','FontSize',25);
set(gca,'FontWeight','bold','FontSize',25);
c=colorbar('Fontsize',20);
c.Title.String = 'TECU';
grid on
shading interp 
load coastlines
geoshow(coastlat,coastlon,'color','k')
hold on

%% ----------------------------Figure Two---------------------------------%

%----------Extracting the Global TEC (GTEC) from the IGS .nc data-------

addpath('D:\Projets Master\Tec_map\TEC_2024');
FolderName = 'D:\Projets Master\Tec_map\TEC_2024\';
fileNames = dir(fullfile(FolderName, '*.cdf'));

GTS = zeros(366, 96);                                % Output: 365 days × 96 time steps

for day = 1:366
    currName = fileNames(day).name;
    currFileName = fullfile(FolderName, currName);
    data = cdfread(currFileName);

    for t = 1:96
        curr_map = data{t, 2};                       % 71×73 TEC map
        global_mean = mean(curr_map(:), 'omitnan');  
        GTS(day, t) = global_mean;
    end
end

%% --------------------//----------------//----------------//------------------//----------------//
%---------------Stack plot of GTEC, Dst, and Kp within 2010-2024----------%

filtered_data  = readtable("D:\Second year\First semester\Statistics\Stat_coding\Global_Dst.csv");

an24 =  filtered_data.GTEC24;
an23 =  filtered_data.GTEC23;
an22 =  filtered_data.GTEC22;
an21 =  filtered_data.GTEC21;
an20 =  filtered_data.GTEC20;
an19 =  filtered_data.GTEC19;
an18 =  filtered_data.GTEC18;
an17 =  filtered_data.GTEC17;
an16 =  filtered_data.GTEC16;
an15 =  filtered_data.GTEC15;
an14 =  filtered_data.GTEC14;
an13 =  filtered_data.GTEC13;
an12 =  filtered_data.GTEC12;
an11 =  filtered_data.GTEC11;
an10 =  filtered_data.GTEC10;

mydata = [an10,an11,an12,an13,an14,an15,an16,an17,an18,an19,an20,an21,an22,an23,an24];
mydata = fillmissing(mydata, 'spline'); 

Meantec = mean(mydata,1);
Meantec(Meantec==0)=NaN;
mydata2 = mean(mydata,2);
smooth = smoothdata(mydata2, 'movmean', 15);

figure(1)
rect = [ 0.2  0.2  0.7  0.7 ];
axes( 'position' , rect);
scalex = 1:15;
scaley2 = 1:12;
scaley = linspace(1,366,366);    % modify here (1,365,365)
scaley1 = linspace(datetime('Jan 01, 2010 00:00:00'),datetime('Dec 31, 2024 00:00:00'),length(mydata));
%-------------II--------------%

h =imagesc(scalex, scaley, mydata);

axis xy;
title('\bf b) Dst distribution from 2010 to 2024','FontSize',22);
set(gcf,'color','w');
set(gca,'FontWeight','bold','color', 'w','FontSize',18);
set(gca, 'xtick',[ ] );
set(gca, 'ytick',[ ]);
pos_upper=get(gca,'pos');
colormap jet;
C=colorbar;
C.Limits = [-300 25];   % Kp = [0 5], GTEC = [0 50]
set(C,'FontSize',26);
ylabel(C,'\bf Dst (nT)','FontWeight','bold','FontSize',18);
rect = [ 0.1  0.2  0.1  0.7 ];
axes( 'position' , rect );

% left graph

plot(mydata2, scaley,'Linewidth',2);
hold on
plot(smooth, scaley,'r','Linewidth',2);
month_days = [1 32 60 91 121 152 182 213 244 274 305 336];


month_labels = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
ylim([1 366]);            % modify here depending of the year (365/366)
yticks(month_days);
yticklabels(month_labels);

set(gca,'FontWeight','bold','FontSize',18);
set(gca, 'xdir','reverse','FontWeight','bold','FontSize',18);
ylabel('\bf Months of the year','FontWeight','bold','FontSize',18);
xlabel('\bf Dst (nT)','FontWeight','bold','FontSize',18);
% yticklabels([]);
grid on;
set(gca,'XMinorGrid','on');
set(gca,'YMinorGrid','on');

% the botom graph

scale = scalex;
rect = [ 0.2 , 0.099 , 0.662 , 0.1 ];
axes( 'position' , rect );
plot(scale, Meantec,'bo-','Linewidth',2);
set(gca,'FontWeight','bold','FontSize',18);
axis tight;
set(gca, 'ytick',[ ] );
xticks(1:1:15);
xticklabels({'2010','2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018', '2019', '2020', '2021', '2022','2023','2024'});
xlabel('\bf Time (Years)','FontWeight','bold','FontSize',18);
grid on;
set(gca,'XMinorGrid','on');
set(gca,'YMinorGrid','on');

%% ---------------------------Figure three -----------------------------%

filtered_data  = readtable("D:\Second year\First semester\Statistics\Stat_coding\Storms\cycle.csv");

seriesname = {'GTEC' 'Dst' 'Kp'};

d3 = filtered_data.Kp;
d2 = filtered_data.Dst;
d1 = filtered_data.GTEC;
d1= fillmissing(d1, 'spline');

filtered_data1  = readtable("D:\Second year\First semester\Statistics\Stat_coding\Storms\puissance.csv");                                         %filtered_data.days';
GTECq =  filtered_data1.p1;
kpInd =  filtered_data1.p2;
DstInd =  filtered_data1.p3;
T =  filtered_data1.t;

%------------------Continous Wavelet Transform (Grinted et Al. Software)----------------------------%

subplot(3,2,1);
wt(d1);
fontsize(gca, 15,'points') 
set(gca,'FontWeight','bold','FontSize',15);
title('a) Wavelet scalogram: GTEC','FontSize',15);
xticklabels([])
ylabel('Period (days)')
colormap('jet');
c=colorbar();

subplot(3,2,2)
plot(T, GTECq,LineStyle="-",LineWidth=2,Color='b');
axis tight 
ylabel('Norm power','FontSize',15);
set(gca,'FontWeight','bold','FontSize',15);
%title('Geomagnetic storm of 03-08/04/2010','FontSize',15);
title('b) Wavelet normalized power spectrum: GTEC','FontSize',15);
xticklabels([])
hold off;
grid on

subplot(3,2,3);
wt(d2);
fontsize(gca, 15,'points') 
set(gca,'FontWeight','bold','FontSize',15);
% set(gca,'XLim')
% xticks_years = 2010:2024;
% xticks_days = (xticks_years - 2010) * 365.25;
% set(gca, 'XTick', xticks_days);
% set(gca, 'XTickLabel', string(xticks_years));
title('c) Wavelet scalogram: Kp','FontSize',15);
xticklabels([])
ylabel('Period (days)')
colormap('jet');
c=colorbar();

subplot(3,2,4)
plot(T, kpInd,LineStyle="-",LineWidth=2,Color='b');
axis tight 
ylabel('Norm power','FontSize',15);
set(gca,'FontWeight','bold','FontSize',15);
%title('Geomagnetic storm of 03-08/04/2010','FontSize',15);
title('d) Wavelet normalized power spectrum: Kp','FontSize',15);
xticklabels([])
hold off;
grid on

subplot(3,2,5);
wt(d3);
fontsize(gca, 15,'points') 
set(gca,'FontWeight','bold','FontSize',15);
set(gca,'XLim')
xticks_years = 2010:2024;
xticks_days = (xticks_years - 2010) * 365.25;
set(gca, 'XTick', xticks_days);
set(gca, 'XTickLabel', string(xticks_years));
title('e) Wavelet scalogram: Dst','FontSize',15);
ylabel('Period (days)')
xlabel('Time (years)')
colormap('jet');
c=colorbar();

subplot(3,2,6)
plot(T, DstInd, LineStyle="-", LineWidth=2, Color='b');
axis tight 
ylabel('Norm power','FontSize',15);
set(gca,'FontWeight','bold','FontSize',15);
title('f) Wavelet normalized power spectrum: Dst','FontSize',15);
xlabel('Period (days)')
hold off;
grid on

%--------------Correlation of GTEC and Dst/Kp normalized power plot------------%

filtered_data  = readtable("D:\Second year\First semester\Statistics\Stat_coding\Storms\power.csv");                                         %filtered_data.days';
GTEC =  filtered_data.PGTEC;
Dst =  filtered_data.PDst;
Kp =  filtered_data.PKp;

figure(1);
subplot(1,2,1)
scatter(Dst, GTEC,'black');
xlim([0 20])
set(gca,'FontWeight','bold','FontSize',15);
ylabel('Normalized GTEC-power ','FontSize',15);
xlabel('Normalized Dst-power ','FontSize',15);
grid on
box on;

subplot(1,2,2)
scatter(Kp,GTEC,'black');
yticklabels([])
set(gca,'FontWeight','bold','FontSize',15);
xlabel('Normalized Kp-power','FontSize',15);
xlim([0 20])
grid on
box on;

R = corrcoef(GTEC, Dst);
correlation = R(1,2);
disp(['Correlation = ', num2str(correlation)]);

%% -----------------------------Figure four-------------------------------%

%--------------------------Latitudinal Daily Mean----------------------%

addpath('D:\Projets Master\Tec_map\TEC_2010');
FolderName = 'D:\Projets Master\Tec_map\TEC_2010\';
fileNames = dir(fullfile(FolderName, '*.cdf'));  
latitudes = linspace(87.5, -87.5, 71);           
target_lats = 80:-10:-80;                        
n_days = 365; n_epochs = 96;                     

GTECE = zeros(n_days, n_epochs, length(target_lats));  % Day × Epoch × Latband

for i = 1:n_days
    currFileName = fullfile(FolderName, fileNames(i).name);
    data = cdfread(currFileName);

    for j = 1:n_epochs
        epoch_map = data{j, 2};    % 71×73 map at epoch t

        for k = 1:length(target_lats)
            lat_center = target_lats(k);
            lat_band_idx = find(latitudes >= (lat_center - 10) & latitudes <= (lat_center + 10));
            band_data = epoch_map(lat_band_idx, :);                                         % Submatrix of lat band
            GTECE(i, j, k) = mean(band_data(:), 'omitnan');
        end
    end
end

%-------------------Yearly periods (near 27-day) trend 2010-2024----------%

filtered_data  = readtable("D:\Second year\First semester\Statistics\Stat_coding\YearsPeriods.csv");

GTEC =  filtered_data.GTECp;
Low =  filtered_data.Lowp; 
Mid =  filtered_data.Midp; 
High =  filtered_data.Highp; 
kp =  filtered_data.Kp;
dst =  filtered_data.Dst;
x = 1:15;

figure(1);
subplot(2,1,1)
plot(x, kp, 'c-o', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
plot(x, dst, 'm-o', 'LineWidth', 1.5, 'MarkerSize', 8);
xticklabels([])
title('Periods trend over the years (2010-2024)','FontSize',15);
set(gca,'FontWeight','bold','FontSize',15);
ylabel('Periods ','FontSize',15);
legend('Kp index','Dst index')
grid on
%scalex = linspace(datetime('Jan 01, 2010 00:00:00'),datetime('Dec 31, 2024 00:00:00'),length(GTEC));
subplot(2,1,2)
plot(x, GTEC, 'k-s', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
plot(x, Low, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
plot(x, Mid, '-v',  'LineWidth', 1.5, 'MarkerSize', 8, color = '[0.13 0.55 0.13]');
hold on;
plot(x, High, '-*','LineWidth', 1.5, 'MarkerSize', 8, color = '[0.85 0.44 0.10]');
ylabel('Periods','FontSize',15);
set(gca,'FontWeight','bold','FontSize',15);
xlabel('Time (Years)','FontSize',15);
% title('Periods trend over the years','FontSize',15);
box on;
grid on 
legend('GTEC','Low band ([0, 20N])','Middle band ([30N, 50N])','High band ([60N, 80N])')
xticks(1:1:14);
xticklabels({'2010','2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018', '2019', '2020', '2021', '2022','2023','2024'});

%---------------------------Correlation-----------------------------%

filtered_data1  = readtable("D:\Second year\First semester\Statistics\Stat_coding\NBStorms.csv");
nb = filtered_data1.Nstorms;
R = corrcoef(nb, dst);
correlation = R(1,2);
disp(['Correlation = ', num2str(correlation)]);

%-------------------------Bar plot number of storms and classes-----------%
%---Storms count from (https://www.spaceweatherlive.com/, https://catalogs.astro.bas.bg/)

x = 2010:2024;
y = [8 0 0; 10 3 0; 21 6 0; 22 3 0; 18 1 0; 35 9 0; 22 3 0; 9 3 0; 5 1 0; 4 0 0; 4 0 0; 10 1 0; 26 0 0; 36 8 0; 34 12 4];
b = bar(x,y,'stacked');
myColors = [ ...
    1 0.5 0;    
    1 0 0;      
    0.5 0 0;    
];

for i = 1:length(b)
    b(i).FaceColor = myColors(i,:);
end
set(gca,'FontWeight','bold','FontSize',15);
ylabel('Nbr of Storms');
xlabel('Time (years)');
title('Number of storms throughout the period (2010-2024)');
legend({'G3-Storms', 'G4-Storms', 'G5-Storms'});

%% ----------------------------Figure 5, 6, and 7 ------------------------%

%----------27-day runing median to have the enhanced data sets------------%

filtered_data  = readtable("D:\Second year\First semester\Statistics\Stat_coding\GTEC15mins.csv");
GTEC =  filtered_data.GTEC11; 
% filtered_data  = readtable("D:\Second year\First semester\Statistics\Stat_coding\Storms\Indices24.csv");
% GTEC =  filtered_data.Dst; 

% startDate = datetime(2024, 1, 1, 0, 0, 0);                                  
% endDate = datetime(2024, 12, 31, 23, 0, 0); 
startDate = datetime(2011, 1, 1, 0, 0, 0);
endDate = datetime(2011, 12, 31, 23, 45, 0);
% timeVec = (startDate:hours(1):endDate)';
timeVec = (startDate:minutes(15):endDate)';
n = length(timeVec);

%------------ Define storm period (May 9–14)-------------%
stormStart = datetime(2011,10,22,0,0,0);
stormEnd = datetime(2011,10,27,23,45,0);
% stormStart = datetime(2024,05,08,0,0,0);
% stormEnd = datetime(2024,05,13,23,0,0);
stormIdx = timeVec >= stormStart & timeVec <= stormEnd;

%-------------- Define quiet period (27 days before storm)-------%
quietStart = stormStart - days(26);
quietEnd = stormStart - minutes(15);
% quietEnd = stormStart - hours(1);
quietIdx = timeVec >= quietStart & timeVec <= quietEnd;

% Extract time-of-day for matching
quietTimeOfDay = timeofday(timeVec(quietIdx));
stormTimeOfDay = timeofday(timeVec(stormIdx));

% Build 27-day median reference curve
stormTEC = GTEC(stormIdx);
quietTEC = GTEC(quietIdx);

medianCurve = zeros(size(stormTEC));
for i = 1:length(stormTimeOfDay)
    matchIdx = abs(minutes(quietTimeOfDay - stormTimeOfDay(i))) < 1e-6;
    if any(matchIdx)
        medianCurve(i) = median(quietTEC(matchIdx));
    else
        medianCurve(i) = NaN; 
    end
end

% % Plotting
% figure;
% plot(timeVec(stormIdx), stormTEC, 'k', 'LineWidth', 1.5); hold on;
% plot(timeVec(stormIdx), medianCurve, 'r--', 'LineWidth', 1.5);
% legend('TEC','TEC quiet reference');
% ylabel('TEC (TECu)');
% xlabel('Time UTC');
% title('TEC station during storm period');
% grid on;

delta = stormTEC-medianCurve;
data_reshaped = reshape(delta, 4, []);
mean_vector = mean(data_reshaped, 1)';

%--------------------Plot of the enhanced GTEC, Dst, and Kp---------------%

%filtered_data  = readtable("D:\Second year\First semester\Statistics\Stat_coding\Storms\ST17.csv");
filtered_data1  = readtable("D:\Second year\First semester\Statistics\Stat_coding\Storms\STGTEC.csv");
filtered_data  = readtable("D:\Second year\First semester\Statistics\Stat_coding\Storms\Orage24.csv");
Dst =  filtered_data.Dst;
Kp =  filtered_data.Kp;   
% bz =  filtered_data.Bz; 
% MedDst =  filtered_data.MedDst; 
% MedKp =  filtered_data.MedKp; 
Med =  filtered_data1.Med24;
% med = filtered_data1.med11;
GTEC =  filtered_data1.S24;
% SWSpeed =  filtered_data.SWSpeed; 
delta = GTEC-Med;
% delta1 = Kp-MedKp
% delta11 = GTEC-med;
% delta2 = Dst-MedDst;
% ---------------
b2 = datetime("2024-05-10 17:00:00");
a1 = datetime("2024-05-10 22:00:00");    % tec
a2 = datetime("2024-05-11 02:00:00");    % dst
a3 = datetime("2024-05-11 00:00:00");    % kp
t1 =  datetime("2024-05-10 17:00:00"); % SSC time
t2 =  datetime("2024-05-12 08:00:00"); % recovery start

scalex = linspace(datetime('May 08, 2024 00:00:00'),datetime('May 13, 2024 23:00:00'),length(Dst));
scalex1 = linspace(datetime('May 08, 2024 00:00:00'),datetime('May 13, 2024 23:45:00'),length(GTEC));

subplot(3,1,1)

plot(scalex, Kp./10,LineStyle="-",LineWidth=2,Color='b');
axis tight 
ylabel('Kp','FontSize',15);
set(gca,'FontWeight','bold','FontSize',15);
title('Geomagnetic storm of 08-13/05/2024','FontSize',15);
xticklabels([])
hold on;
xline(b2, 'r--', 'LineWidth', 2);
xline(a3, 'r--', 'LineWidth', 2);
yline(0,'--');
hold off;
grid on


subplot(3,1,2)

plot(scalex, Dst,LineStyle="-",LineWidth=2,Color='b');
axis tight 
ylabel('Dst (nT)','FontSize',15);
set(gca,'FontWeight','bold','FontSize',15);
xticklabels([])
hold on;
xline(b2, 'r--', 'LineWidth', 2);
xline(a2, 'r--', 'LineWidth', 2);
yline(0,'--');

box on;
grid on

subplot(3,1,3)

plot(scalex1, delta,LineStyle="-",LineWidth=2,Color='b');
axis tight 
ylabel('GTEC (TECU)','FontSize',15);
set(gca,'FontWeight','bold','FontSize',15);
xlabel('Time (days)','FontSize',15);
hold on;
xline(b2, 'r--', 'LineWidth', 2);
xline(a1, 'r--', 'LineWidth', 2);
yline(0,'--');
box on;
grid on 
 %--------initiate the gray part------------------%
for i = 1:3
    subplot(3,1,i)   
    hold on
    yl = ylim;       
    fill([t1 t2 t2 t1], [yl(1) yl(1) yl(2) yl(2)], ...
         [0.7 0.7 0.7], 'FaceAlpha',0.3, 'EdgeColor','none')
end

%% ------------Wavelet Coherence (WC, Grinted et Al Software)---------------%
%--------------------------Enhanced time series-----------------------%
filtered_data  = readtable("D:\Second year\First semester\Statistics\Stat_coding\Storms\Orage24.csv");

seriesname = {'GTEC' 'Dst' 'Kp'};

z3 = filtered_data.Kp;
z2 = filtered_data.Dst;
z11 = filtered_data.gtec; % Enhanced GTEC

figure('color',[1 1 1])
wtc(z11,z3);   %,'mcc',0);


set(gca,'FontWeight','bold','FontSize',15);

startDate = datetime(2024,5,08);   % adjust year as needed
endDate = datetime(2024,5,13);
xticks = linspace(1, length(z1), 6);  
xticklabels({'May 08', 'May 09', 'May 10', 'May 11', 'May 12', 'May 13'});
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
colormap('jet');
c=colorbar();
fontsize(gca, 15,'points') 

xlabel('Time (days)')
ylabel('Period (hours)')
title(['d) Wavelet coherence: ' seriesname{1} '-' seriesname{3} ] )

%-------------------------Phase average calculation-----------------------%
 
t=(1:1:144)';
X=   z11;                
Y=   z3;               
[Wxy,period,scale,coi,sig95]=xwt([t X], [t Y]);
[mn,rowix]=min(abs(period-27.7879)); 

incoi=(period(:)*(1./coi)>1);
issig=(sig95>=1);
angles=angle(Wxy(rowix,issig(rowix,:)&~incoi(rowix,:)));
[meantheta,anglestrength,sigma]=anglemean(angles);

D1 = rad2deg(meantheta);
D2 = rad2deg(sigma);

 


