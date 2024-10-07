%% Permafrost Area Fraction Plot for introduction

clear
clc

size1 = 30;

% Load the NC file
file_path = 'C:\\Users\\mohamed.ahajjam\\Desktop\\UND\\Defense resiliency platform\\Month23\\Fire\\Maps\\ESACCI-PERMAFROST-L4-PFR-ERA5_MODISLST_BIASCORRECTED-AREA4_PP-2000-fv04.0.nc';

% Read variables
lat = ncread(file_path, 'lat');
lon = ncread(file_path, 'lon');
% pfr = ncread(file_path, 'PFR');
load('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\results\PFR_all_years.mat') %just load this!
% Remove fill values (assuming 255 is the fill value for missing data)
PFR_all_years(PFR_all_years == 255) = NaN;
% Calculate the modal value for each grid cell
pfr = mean(PFR_all_years, 3);
filtered_pfr = pfr';

% Define the Alaska region
lat_min = 50; lat_max = 72;
lon_min = -170; lon_max = -130;

% Convert longitudes to the range [-180, 180]
lon = mod(lon + 180, 360) - 180;

% Filter the data within the Alaska region
lat_idx = find(lat >= lat_min & lat <= lat_max);
lon_idx = find(lon >= lon_min & lon <= lon_max);

filtered_lat = lat(lat_idx);
filtered_lon = lon(lon_idx);
% filtered_pfr = pfr(lon_idx, lat_idx)';

% Mask out values outside of Alaska
masked_pfr = filtered_pfr;
masked_pfr(masked_pfr == 0) = NaN; % Assuming 0 is the value to mask

% Create a meshgrid for the filtered data
[lon_grid, lat_grid] = meshgrid(filtered_lon, filtered_lat);

lat_grid = lat_grid(:,1:2900); %cheat to remove data outside of alaska
lon_grid = lon_grid(:,1:2900); %cheat to remove data outside of alaska
masked_pfr = masked_pfr(:,1:2900); %cheat to remove data outside of alaska

figure('units','normalized','outerposition',[0 0 1 1]);
% Create a geoaxes
gx = geoaxes;
geobasemap(gx, 'streets-light'); % You can choose other basemap options like 'streets', 'topographic', etc.

% Plot the PFR data using geoscatter (point-based) for demonstration
geoscatter(lat_grid(:), lon_grid(:), 5, masked_pfr(:), 'filled');

% Adjust color limits and add colorbar
load('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\results\CustomColormap_ice2.mat');
colormap(CustomColormap_ice2);

c = colorbar;
c.Label.String = "Permafrost Area Fraction";
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
c.TickLabels = strcat(c.TickLabels,'$\%$');
c.FontSize = size1-10;

gx.LongitudeAxis.FontSize = size1-10;
gx.LongitudeAxis.TickLabelInterpreter = 'latex';
gx.LongitudeAxis.TickLabels = strrep(get(gx.LongitudeAxis,'TickLabels'),'°','$^{\circ}$');
gx.LatitudeAxis.FontSize = size1-10;
gx.LatitudeAxis.TickLabelInterpreter = 'latex';
gx.LatitudeAxis.TickLabels = strrep(get(gx.LatitudeAxis,'TickLabels'),'°','$^{\circ}$');

hold on

% Plot boundaries of Alaska
shapefile = 'usastatelo.shp';
states = readgeotable(shapefile);
alaska = states(strcmp(states.Name, 'Alaska'), :);
h = geoplot(alaska, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);


% Plot wildfire events as red dots
cmap = get(0, 'defaultaxescolororder');
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork',"Fire_archive_AlaskaOnly.csv");
opts = detectImportOptions(path1);
MODIS = readtable(path1,opts);
geoscatter(MODIS.latitude(year(MODIS.acq_date)<2022), MODIS.longitude(year(MODIS.acq_date)<2022), 2,cmap(2,:),'filled');

% Create legend
legend1 = legend('FIRMS Wildfire detected events');
set(legend1,'Interpreter','latex','FontSize',size1-10);

%%make sure to change the interpereter of the x and y axis labels manually
%%for now :p

% exportgraphics(gcf,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\PFR_FIRMS_plot_v2.eps'),"ContentType","vector")
%% Run this after making the changes in the figure of the previous section
ExportGEOTIFF(gca, 'PFR_FIRMS_plot_v2')

%% wildfire occurrence prediction as figure - Task1
clc
size1 = 35;
%accuracy
Season_naive_baseline   = [53.07	51.99	48.98	53.31];
Location_naive_baseline	= [49.82	49.82	49.82	49.82];
Season_Location_HGB	    = [59.81	58.60	58.12	59.33];
% LogisticRegression= [50.78	50.54	50.18	49.58];
KNeighbors	      = [62.09	61.97	61.13	57.64];
% DecisionTree 	  = [59.69	58.36	58.24	56.80];
% RandomForest 	  = [62.33	63.18	61.73	57.40];
% LGBM         	  = [94.58	94.34	94.71	93.38];
SVC          	  = [55.96	53.79	53.79	55.84];
MLP          	  = [57.28	58.00	57.40	55.11];
TSiT         	  = [50.42	50.42	50.42	50.66];
ConvTran     	  = [50.18	50.18	50.18	50.54];
TSSequencer  	  = [49.94	49.94	49.94	50.18];

%Average improvement over any baseline: 
% mean(([94.71	94.46	93.50	93.38] - mean([Season_naive_baseline;Location_naive_baseline;Season_Location_HGB;LogisticRegression;KNeighbors;DecisionTree;RandomForest;SVC;MLP;TSiT;ConvTran;TSSequencer]))./mean([Season_naive_baseline;Location_naive_baseline;Season_Location_HGB;LogisticRegression;KNeighbors;DecisionTree;RandomForest;SVC;MLP;TSiT;ConvTran;TSSequencer]))


Proposed = [94.71	94.46	93.50	93.38];

figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact','units','normalized','outerposition',[0 0 1 1]);
nexttile%LogisticRegression;DecisionTree;RandomForest;

b1 = bar(([0,1,2,3]),[KNeighbors;SVC;MLP;TSiT;ConvTran;TSSequencer;Proposed],'EdgeColor','none','FaceAlpha',0.6);
ylabel('Accuracy ($\%$)','FontSize',size1-12,'Interpreter','latex')
hold on,
y1 = yline(mean([Season_naive_baseline;Location_naive_baseline;Season_Location_HGB],'all'),':','Interpreter','latex',...
'LabelHorizontalAlignment','left','LabelVerticalAlignment','top',...
'FontSize',size1-12,'LineWidth',2, 'Color',[0 0 0]);

% yline(mean(Season_Location_HGB,'all'),'-.','Avg. Season-Location HGB','HandleVisibility','off','Interpreter','latex',...
% 'LabelHorizontalAlignment','left','LabelVerticalAlignment','top',...
% 'FontSize',size1-10,'LineWidth',1.5, 'Color',[0 0 0]);

set(gca,'ylim',[40 100])
set(gca, 'XTickLabel', {'+0','+7','+30','+90'},'TickLabelInterpreter','latex')

set(gca,'Box','on','FontSize',size1-12,'TickLabelInterpreter','latex','LineWidth',0.5,'YMinorGrid','on','YMinorTick','on','TickDir','in','YGrid','on');
% legend1 = legend({'LogisticRegression','KNeighbors','DecisionTree','RandomForest','SVC','MLP','TSiT','ConvTran','TSSequencer','Proposed Model'});
% set(legend1,'Interpreter','latex','FontSize',size1,'Orientation','horizontal','Location','northeastoutside','NumColumns',1);
% title(legend1,'Technique');
colororder(gca,"gem12")

xtips1 = b1(end).XEndPoints;
ytips1 = b1(end).YEndPoints-2.5;
labels1 = strcat(string(b1(end).YData),'\%');
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','middle', 'FontSize',size1-20,'Interpreter','latex','Rotation',90)

% exportgraphics(gcf,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\ResultsAccuracy_task1.eps'),"ContentType","vector")
% 
% %% wildfire occurrence prediction as figure - Task1|FAR

clc
size1 = 35;
%accuracy
Season_naive_baseline   = [30.19	42.51	35.27	24.40];
Location_naive_baseline	= [0.00	0.00	0.00	0.00];
Season_Location_HGB	    = [35.49	37.17	38.13	35.25];

% LogisticRegression= [49.88	51.08	50.84	52.04];
KNeighbors	      = [37.89	38.37	39.09	42.45];
% DecisionTree 	  = [38.13	42.93	41.01	44.12];
% RandomForest 	  = [37.17	35.73	37.65	40.05];
% LGBM         	  = [5.28	5.76	5.04	5.52];
SVC          	  = [62.35	52.04	50.84	61.63];
MLP          	  = [41.25	42.45	45.32	44.12];
TSiT         	  = [14.39	14.39	14.39	22.06];
ConvTran     	  = [0.00+0.1	0.00+0.1	0.00+0.1	4.32];
TSSequencer  	  = [0.72	0.72	0.72	6.47];
Proposed = [3.84	3.84	5.04	4.56];


nexttile%LogisticRegression;DecisionTree;RandomForest;
b = bar(([0,1,2,3]),[KNeighbors;SVC;MLP;TSiT;ConvTran;TSSequencer;Proposed],'EdgeColor','none','FaceAlpha',0.6);
ylabel('False Alarm Rate ($\%$)','FontSize',size1-12,'Interpreter','latex')
hold on,
yline(mean([Season_naive_baseline;Location_naive_baseline;Season_Location_HGB],'all'),':','Interpreter','latex',...
'LabelHorizontalAlignment','left','LabelVerticalAlignment','top',...
'FontSize',size1-12,'LineWidth',2, 'Color',[0 0 0]);
% 
% yline(mean(Season_Location_HGB,'all'),'-.','Avg. Season-Location HGB','HandleVisibility','off','Interpreter','latex',...
% 'LabelHorizontalAlignment','left','LabelVerticalAlignment','top',...
% 'FontSize',size1-10,'LineWidth',1.5, 'Color',[0 0 0]);

set(gca,'ylim',[0 70])
set(gca, 'XTickLabel', {'+0','+7','+30','+90'},'TickLabelInterpreter','latex')
%'LogisticRegression','DecisionTree','RandomForest',
set(gca,'Box','on','FontSize',size1-12,'TickLabelInterpreter','latex','LineWidth',0.5,'YMinorGrid','on','YMinorTick','on','TickDir','in','YGrid','on');

for i = length(b):length(b)
    xtips1 = b(i).XEndPoints;
    ytips1 = b(i).YEndPoints+2.6;
    labels1 = strcat(string(b(i).YData),'\%');
    text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','middle', 'FontSize',size1-20,'Interpreter','latex','Rotation',90)
end
xlabel(t,'Prediction horizons (days)','FontSize',size1-10,'Interpreter','latex')


legend1 = legend([y1 b1],{'Avg. baselines', 'KNeighbors','SVC','MLP','TSiT','ConvTran','TSSequencer','Proposed Model'});
set(legend1,'Interpreter','latex','FontSize',size1-10,'Orientation','horizontal');
% title(legend1,'Technique');
colororder(gca,"gem12")
legend1.Layout.Tile = 'North';
% 
exportgraphics(gcf,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Results_task1_v2.eps'),"ContentType","vector")

%% wildfire occurrence prediction as figure - Task2
clc
clear
size1 = 35;


%RMSE
Season_naive_baseline   = [204.83	204.17	204.99	205.30];
Season_Location_HGB	    = [217.90	215.76	215.79	218.86];
% LinearRegression  = [202.86	203.01	204.72	205.51];
% SGD               = [210.93	214.96	211.52	207.95]; 
KNeighbors	      = [187.64	196.87	190.57	200.38];
% RandomForest 	  = [185.32	196.21	193.10	200.46];
SVR          	  = [211.05	210.89	211.19	210.91];
MLP          	  = [214.62	214.48	201.59	212.98];
TSiT         	  = [205.02	205.37 207.24 206.65];
ConvTran     	  = [205.24	205.24 205.52 206.30];
TSSequencer  	  = [207.81	208.22 209.50 208.98];

Proposed = [190.55	195.18	194.51	195.76];

figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact','units','normalized','outerposition',[0 0 1 1]);
nexttile%categorical([0,7,30,90])%LinearRegression;SGD;RandomForest;
b1 = bar([0 1 2 3],[KNeighbors;SVR;MLP;TSiT;ConvTran;TSSequencer;Proposed],'EdgeColor','none','FaceAlpha',0.6);
ylabel('Root Mean Squared Error ($km^2$)','FontSize',size1-10,'Interpreter','latex');
% xlabel('Prediction horizons (days)','FontSize',size1,'Interpreter','latex')
hold on,
y1 = yline(mean([Season_naive_baseline;Season_Location_HGB],'all'),':','Interpreter','latex',...
'LabelHorizontalAlignment','left','LabelVerticalAlignment','top',...
'FontSize',size1-12,'LineWidth',2, 'Color',[0 0 0]);

set(gca,'ylim',[180 216])
set(gca, 'XTickLabel', {'+0','+7','+30','+90'},'TickLabelInterpreter','latex')

set(gca,'Box','on','FontSize',size1-10,'TickLabelInterpreter','latex','LineWidth',0.5,'YMinorGrid','on','YMinorTick','on','TickDir','in','YGrid','on');
% legend1 = legend({'LogisticRegression','KNeighbors','DecisionTree','RandomForest','SVR','MLP','TSiT','ConvTran','TSSequencer','Proposed Model'});
% set(legend1,'Interpreter','latex','FontSize',size1,'Orientation','horizontal','Location','northeastoutside','NumColumns',1);
% title(legend1,'Technique');
colororder(gca,"gem12")

xtips1 = b1(end).XEndPoints;
ytips1 = b1(end).YEndPoints+2.4;
labels1 = strcat(string(b1(end).YData),'$km^2$');
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','middle', 'FontSize',size1-20,'Interpreter','latex','Rotation',90)

% exportgraphics(gcf,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\ResultsRMSE_task2.eps'),"ContentType","vector")
% 
% %% wildfire occurrence prediction as figure - Task1|FAR

clc
size1 = 35;
%MASE

Season_naive_baseline   = [0.81	0.80	0.81	0.82];
Season_Location_HGB	    = [0.84	0.82	0.83	0.85];
% LinearRegression  = [0.81	0.81	0.83	0.82];
% SGD               = [0.54	0.74	0.53	0.63]; 
KNeighbors	      = [0.69	0.71	0.68	0.69];
% RandomForest 	  = [0.70	0.73	0.75	0.74];
SVR          	  = [0.56	0.56	0.56	0.56];
MLP          	  = [0.80	0.72	0.72	0.89];
TSiT         	  = [0.77	0.76	0.69	0.77];
ConvTran     	  = [0.86	0.87	0.76	0.87];
TSSequencer  	  = [0.64	0.62	0.59	0.62];

Proposed = [0.66	0.72	0.74	0.72];

nexttile%categorical([0,7,30,90])%LinearRegression;SGD;RandomForest;
b = bar([0 1 2 3],[KNeighbors;SVR;MLP;TSiT;ConvTran;TSSequencer;Proposed],'EdgeColor','none','FaceAlpha',0.6);
ylabel('Mean Absolute Scaled Error','FontSize',size1-10,'Interpreter','latex')
% xlabel('Prediction horizons (days)','FontSize',size1,'Interpreter','latex')
hold on,
yline(mean([Season_naive_baseline;Season_Location_HGB],'all'),':','Interpreter','latex',...
'LabelHorizontalAlignment','left','LabelVerticalAlignment','top',...
'FontSize',size1-12,'LineWidth',2, 'Color',[0 0 0]);

set(gca,'ylim',[0.4 1])
set(gca, 'XTickLabel', {'+0','+7','+30','+90'},'TickLabelInterpreter','latex')

set(gca,'Box','on','FontSize',size1-10,'TickLabelInterpreter','latex','LineWidth',0.5,'YMinorGrid','on','YMinorTick','on','TickDir','in','YGrid','on');
% legend1 = legend({'LinearRegression','SGD','KNeighbors','RandomForest','SVR','MLP','TSiT','ConvTran','TSSequencer','Proposed Model'});
% set(legend1,'Interpreter','latex','FontSize',size1-10,'Orientation','horizontal','Location','northeastoutside','NumColumns',1);
% title(legend1,'Technique');
% colororder(gca,"gem12")

for i = length(b):length(b)
    xtips1 = b(i).XEndPoints;
    ytips1 = b(i).YEndPoints+.015;
    labels1 = (string(b(i).YData));
    text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','middle', 'FontSize',size1-20,'Interpreter','latex','Rotation',90)
end

xlabel(t,'Prediction horizons (days)','FontSize',size1-10,'Interpreter','latex')

legend1 = legend([y1 b1],{'Avg. baselines', 'KNeighbors','SVR','MLP','TSiT','ConvTran','TSSequencer','Proposed Model'});
set(legend1,'Interpreter','latex','FontSize',size1-10,'Orientation','horizontal');
% title(legend1,'Technique');
colororder(gca,"gem12")
legend1.Layout.Tile = 'North';
% 
exportgraphics(gcf,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Results_task2_v2.eps'),"ContentType","vector")

%% wildfire duration prediction as figure - Task3
clc
clear
size1 = 35;

%RMSE
Season_naive_baseline   = [9.10	9.04	9.11	9.08];
Season_Location_HGB	    = [8.89	9.05	9.08	8.97];
% LinearRegression  = [8.83	8.82	8.91	8.93];
% SGD               = [10.22	9.86	9.51	9.98]; 
KNeighbors	      = [7.28	6.92	6.98	6.85];
% RandomForest 	  = [6.82	6.89	6.92	7.33];
SVR          	  = [9.28	9.46	9.59	9.58];
MLP          	  = [8.48	7.71	7.48	7.98];
TSiT         	  = [8.22	8.37 8.26 8.38];
ConvTran     	  = [8.94	8.91	9.11	9.10];
TSSequencer  	  = [10.20	10.28	9.29	10.28];

Proposed = [6.68	7.14	6.69	7.03];

figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact','units','normalized','outerposition',[0 0 1 1]);
nexttile%categorical([0,7,30,90])%LinearRegression;SGD;RandomForest;
b1 = bar([0 1 2 3],[KNeighbors;SVR;MLP;TSiT;ConvTran;TSSequencer;Proposed],'EdgeColor','none','FaceAlpha',0.6);
ylabel('Root Mean Squared Error (days)','FontSize',size1-10,'Interpreter','latex')
% xlabel('Prediction horizons (days)','FontSize',size1,'Interpreter','latex')
hold on,
y1 = yline(mean([Season_naive_baseline;Season_Location_HGB],'all'),':','Interpreter','latex',...
'LabelHorizontalAlignment','left','LabelVerticalAlignment','top',...
'FontSize',size1-12,'LineWidth',2, 'Color',[0 0 0]);

set(gca,'ylim',[6 10.5])
set(gca, 'XTickLabel', {'+0','+7','+30','+90'},'TickLabelInterpreter','latex')

set(gca,'Box','on','FontSize',size1-10,'TickLabelInterpreter','latex','LineWidth',0.5,'YMinorGrid','on','YMinorTick','on','TickDir','in','YGrid','on');
% legend1 = legend({'LogisticRegression','KNeighbors','DecisionTree','RandomForest','SVC','MLP','TSiT','ConvTran','TSSequencer','Proposed Model'});
% set(legend1,'Interpreter','latex','FontSize',size1,'Orientation','horizontal','Location','northeastoutside','NumColumns',1);
% title(legend1,'Technique');
colororder(gca,"gem12")

xtips1 = b1(end).XEndPoints;
ytips1 = b1(end).YEndPoints+.26;
labels1 = strcat(string(b1(end).YData),' days');
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','middle', 'FontSize',size1-20,'Interpreter','latex','Rotation',90)

% exportgraphics(gcf,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\ResultsRMSE_task3.eps'),"ContentType","vector")
% 
% %% wildfire occurrence prediction as figure - Task1|FAR

clc
size1 = 35;
%MASE

Season_naive_baseline   = [0.80	0.80	0.80	0.80];
Season_Location_HGB	    = [0.73	0.73	0.73	0.73];
% LinearRegression  = [0.78	0.78	0.79	0.78];
% SGD               = [0.65	0.63	0.67	0.67]; 
KNeighbors	      = [0.59	0.58	0.59	0.58];
% RandomForest 	  = [0.57	0.58	0.58	0.60];
SVR          	  = [0.59	0.60	0.61	0.61];
MLP          	  = [0.75	0.73	0.70	0.75];
TSiT         	  = [0.72 0.73 0.76 0.75];
ConvTran     	  = [0.77	0.77	0.80	0.79];
TSSequencer  	  = [0.68	0.68	0.71	0.71];

Proposed = [0.52	0.57	0.53	0.55];

nexttile%categorical([0,7,30,90])%LinearRegression;SGD;RandomForest;
b = bar([0 1 2 3],[KNeighbors;SVR;MLP;TSiT;ConvTran;TSSequencer;Proposed],'EdgeColor','none','FaceAlpha',0.6);
ylabel('Mean Absolute Scaled Error','FontSize',size1-10,'Interpreter','latex')
% xlabel('Prediction horizons (days)','FontSize',size1,'Interpreter','latex')
hold on,
yline(mean([Season_naive_baseline;Season_Location_HGB],'all'),':','Interpreter','latex',...
'LabelHorizontalAlignment','left','LabelVerticalAlignment','top',...
'FontSize',size1-12,'LineWidth',2, 'Color',[0 0 0]);

set(gca,'ylim',[0.4 1])
set(gca, 'XTickLabel', {'+0','+7','+30','+90'},'TickLabelInterpreter','latex')

set(gca,'Box','on','FontSize',size1-10,'TickLabelInterpreter','latex','LineWidth',0.5,'YMinorGrid','on','YMinorTick','on','TickDir','in','YGrid','on');

for i = length(b):length(b)
    xtips1 = b(i).XEndPoints;
    ytips1 = b(i).YEndPoints+.015;
    labels1 = (string(b(i).YData));
    text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','middle', 'FontSize',size1-20,'Interpreter','latex','Rotation',90)
end

xlabel(t,'Prediction horizons (days)','FontSize',size1-10,'Interpreter','latex')

legend1 = legend([y1 b1], {'Avg. baselines','KNeighbors','SVR','MLP','TSiT','ConvTran','TSSequencer','Proposed Model'});
set(legend1,'Interpreter','latex','FontSize',size1-10,'Orientation','horizontal');
% title(legend1,'Technique');
colororder(gca,"gem12")
legend1.Layout.Tile = 'North';

exportgraphics(gcf,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Results_task3_v2.eps'),"ContentType","vector")


%% Fire Detection Prediction Map Plot
clear
clc
size1 = 40;

i=0;

%results (predictions and targets)
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\results',['FireResultsV4_HistGradientBoostingClassifier_hrz' num2str(i) '_GA_testresults.csv']);
opts = detectImportOptions(path1);
testset = readtable(path1,opts);

%target locations
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',['FireDetection_testsetV4_' num2str(i) 'Days.csv']);
opts = detectImportOptions(path1);
targets = readtable(path1,opts);
testset.avgLongitude = targets.avgLongitude;
testset.avgLatitude = targets.avgLatitude;

for j = 1:height(testset)
    switch testset.test_targs(j)
        case 0
            switch isequal(testset.test_targs(j), testset.test_preds(j))
                case 1
                    testset.TN(j) = 1;
                case 0
                    testset.FN(j) = 1;
            end
        case 1
            switch isequal(testset.test_targs(j), testset.test_preds(j))
                case 1
                    testset.TP(j) = 1;
                case 0
                    testset.FP(j) = 1;
            end
    end
end

% % All in one
% size1 = 18;
% figure('units','normalized','outerposition',[0 0 1 1]);
% t = tiledlayout(1,1,'TileSpacing','tight','Padding','compact','units','normalized','outerposition',[0 0 1 1]);
% geoplot(testset.avgLatitude(testset.TP==1), testset.avgLongitude(testset.TP==1), "o", "LineWidth", 1, ...
%     "MarkerSize", 5, "MarkerEdgeColor",[0.8500 0.3250 0.0980], "MarkerFaceColor", [0.8500 0.3250 0.0980], 'DisplayName',strcat('Fire - Correct (',num2str(round(length(find(testset.TP==1)) / [length(find(testset.FP==1))+length(find(testset.TP==1))]*100,2)),'\%)')) % Plot TP
% hold on
% geoplot(testset.avgLatitude(testset.FP==1), testset.avgLongitude(testset.FP==1), "x", "LineWidth", 1, ...
%     "MarkerSize", 10, "MarkerEdgeColor", [0.8500 0.3250 0.0980], "MarkerFaceColor", [0.8500 0.3250 0.0980], 'DisplayName',strcat('Fire - Incorrect (',num2str(round(length(find(testset.FP==1)) / [length(find(testset.FP==1))+length(find(testset.TP==1))]*100,2)),'\%)')) % Plot FP
% geoplot(testset.avgLatitude(testset.TN==1), testset.avgLongitude(testset.TN==1), "o", "LineWidth", 1, ...
%     "MarkerSize", 5, "MarkerEdgeColor", [0 0.4470 0.7410], "MarkerFaceColor", 'None', 'DisplayName',strcat('No fire - Correct (',num2str(round(length(find(testset.TN==1)) / [length(find(testset.FN==1))+length(find(testset.TN==1))]*100,2)),'\%)')) % Plot TN
% geoplot(testset.avgLatitude(testset.FN==1), testset.avgLongitude(testset.FN==1), "x", "LineWidth", 1, ...
%     "MarkerSize", 10, "MarkerEdgeColor", [0 0.4470 0.7410], "MarkerFaceColor", 'None', 'DisplayName',strcat('No fire - Incorrect (',num2str(round(length(find(testset.FN==1)) / [length(find(testset.FN==1))+length(find(testset.TN==1))]*100,2)),'\%)')) % Plot FN
% geobasemap topographic
% legend1 = legend(gca,'show');
% set(legend1,'Interpreter','latex','FontSize',size1);
% set(gca,'Box','on','FontSize',size1,'FontName','helvetica','LineWidth',0.5);

% separate
figure('units','normalized','outerposition',[0 0 1 1]);
t = tiledlayout(1,2,'TileSpacing','tight','Padding','compact','units','normalized','outerposition',[0 0 1 1]);
nexttile
geoplot(testset.avgLatitude(testset.TP==1), testset.avgLongitude(testset.TP==1), "o", "LineWidth", 1,'HandleVisibility','on', ...
    "MarkerSize", 5, "MarkerEdgeColor",[0 0.4470 0.7410], "MarkerFaceColor", [0 0.4470 0.7410], 'DisplayName',strcat('Correct (',num2str(round(length(find(testset.TP==1)) / [length(find(testset.FP==1))+length(find(testset.TP==1))]*100,2)),'\%)')) % Plot TP
hold on
geoplot(testset.avgLatitude(testset.FP==1), testset.avgLongitude(testset.FP==1), "x", "LineWidth", 2,'HandleVisibility','on', ...
    "MarkerSize", 10, "MarkerEdgeColor", [0.8500 0.3250 0.0980], "MarkerFaceColor", [0.8500 0.3250 0.0980], 'DisplayName',strcat('Incorrect (',num2str(round(length(find(testset.FP==1)) / [length(find(testset.FP==1))+length(find(testset.TP==1))]*100,2)),'\%)')) % Plot FP
[limlats,limlons] = geolimits(gca);
legend1 = legend(gca,'show');
set(legend1,'Interpreter','latex','FontSize',size1-10);
geobasemap topographic
title('Fire class','FontSize',size1-15,'Interpreter','latex')
set(gca,'Box','on','FontSize',size1-15,'FontName','helvetica','LineWidth',0.5);
gx = gca;
gx.LongitudeAxis.FontSize = size1-10;
gx.LongitudeAxis.TickLabelInterpreter = 'latex';
gx.LongitudeAxis.TickLabels = strrep(get(gx.LongitudeAxis,'TickLabels'),'°','$^{\circ}$');
gx.LatitudeAxis.FontSize = size1-10;
gx.LatitudeAxis.TickLabelInterpreter = 'latex';
gx.LatitudeAxis.TickLabels = strrep(get(gx.LatitudeAxis,'TickLabels'),'°','$^{\circ}$');

nexttile
geoplot(testset.avgLatitude(testset.TN==1), testset.avgLongitude(testset.TN==1), "o", "LineWidth", 1,'HandleVisibility','on', ...
    "MarkerSize", 5, "MarkerEdgeColor", [0 0.4470 0.7410], "MarkerFaceColor", 'None', 'DisplayName',strcat('Correct (',num2str(round(length(find(testset.TN==1)) / [length(find(testset.FN==1))+length(find(testset.TN==1))]*100,2)),'\%)')) % Plot TN
hold on
geoplot(testset.avgLatitude(testset.FN==1), testset.avgLongitude(testset.FN==1), "x", "LineWidth", 2,'HandleVisibility','on', ...
    "MarkerSize", 10, "MarkerEdgeColor", [0.8500 0.3250 0.0980], "MarkerFaceColor", 'None', 'DisplayName',strcat('Incorrect (',num2str(round(length(find(testset.FN==1)) / [length(find(testset.FN==1))+length(find(testset.TN==1))]*100,2)),'\%)')) % Plot FN
legend1 = legend(gca,'show');
set(legend1,'Interpreter','latex','FontSize',size1-10);
geolimits(limlats,limlons)
geobasemap topographic
title('No fire class','FontSize',size1-15,'Interpreter','latex')
set(gca,'Box','on','FontSize',size1-15,'FontName','helvetica','LineWidth',0.5);
gx2 = gca;
gx2.LongitudeAxis.FontSize = size1-10;
gx2.LongitudeAxis.TickLabelInterpreter = 'latex';
gx2.LongitudeAxis.TickLabels = strrep(get(gx2.LongitudeAxis,'TickLabels'),'°','$^{\circ}$');
gx2.LatitudeAxis.FontSize = size1-10;
gx2.LatitudeAxis.TickLabelInterpreter = 'latex';
gx2.LatitudeAxis.TickLabels = strrep(get(gx2.LatitudeAxis,'TickLabels'),'°','$^{\circ}$');

%manual edit the labels and location of legend first then export
% exportgraphics(gcf,fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2', strcat('Preds_OccurrenceViz','.pdf')),"ContentType","image")
%% Run this after updating figure from previous section
ExportGEOTIFF(gx, 'Preds_OccurrenceViz_1')
ExportGEOTIFF(gx2, 'Preds_OccurrenceViz_2')

%% Fire detection prediction over different horizons
clear
clc
%read data
TASK1_acc = [0.988984088 0.971848225	0.969400245	0.965728274]*100;
TASK1_far = [0.014888337	0.044665012	0.034739454	0.039702233]*100;

size1 = 40;
figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout('flow','TileSpacing','none','Padding','tight','units','normalized','outerposition',[0 0 1 1]);
nexttile,

yyaxis left
b = bar(categorical([0,7,30,90]),TASK1_acc,'EdgeColor','none','FaceAlpha',0.6);
ylabel('Accuracy (\%)','FontSize',size1,'Interpreter','latex','Color','#4698c6')
xlabel('Prediction horizon (days)','FontSize',size1,'Interpreter','latex')

%add text
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = strcat(string(round(b(1).YData,2)),'\%');
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',size1-8,Interpreter='latex',Color='#4698c6')

hold on

yyaxis right
p = plot(categorical([0,7,30,90]),TASK1_far,'LineStyle','-','LineWidth', 3, 'Marker','o', 'MarkerSize',10, 'MarkerFaceColor','w', 'Color','#e8864f');
ylabel('False Alarm Rate (\%)','FontSize',size1,'Interpreter','latex');
set(gca,'Box','on','XTickLabelRotation',0,'FontSize',size1-5,'TickLabelInterpreter','latex','LineWidth',0.5,'YMinorGrid','on','YMinorTick','on','TickDir','in','TickLength',[0.005 0.005],'XGrid','on','YGrid','on');
set(gca,'FontSize',size1,'TickLabelInterpreter','latex');

%add text
labels1 = strcat(string(round(TASK1_far,2)),'\%');
text(categorical([0,7,30,90]),TASK1_far,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',size1-8,Interpreter='latex',Color='w')

ax = gca;
ax.YAxis(1).Limits = [90 100];
ax.YAxis(1).Color = '#4698c6';
ax.YAxis(2).Limits = [0 10];
ax.YAxis(2).Color = '#e8864f';

%export figure
% exportgraphics(t,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month19\Fire\Figure_FireDetectionResults.pdf'),"ContentType","image")

close("all")

%% Fire burnt area prediction over different horizons
clear
clc
%read data
TASK2_rmse = [171.3293202	186.764856	186.7519059	191.344585];
TASK2_mase = [0.633654994	0.709618223	0.6826714	0.667638429];

size1 = 40;
figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout('flow','TileSpacing','none','Padding','tight','units','normalized','outerposition',[0 0 1 1]);
nexttile,

yyaxis left
b = bar(categorical([0,7,30,90]),TASK2_rmse,'EdgeColor','none','FaceAlpha',0.6);
ylabel('RMSE ($km^2$)','FontSize',size1,'Interpreter','latex','Color','#4698c6')
xlabel('Prediction horizon (days)','FontSize',size1,'Interpreter','latex')

%add text
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = strcat(string(round(b(1).YData,2)),'$km^2$');
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',size1-8,Interpreter='latex',Color='#4698c6')

hold on

yyaxis right
p = plot(categorical([0,7,30,90]),TASK2_mase,'LineStyle','-','LineWidth', 3, 'Marker','o', 'MarkerSize',10, 'MarkerFaceColor','w', 'Color','#e8864f');
ylabel('MASE','FontSize',size1,'Interpreter','latex');
set(gca,'Box','on','XTickLabelRotation',0,'FontSize',size1-5,'TickLabelInterpreter','latex','LineWidth',0.5,'YMinorGrid','on','YMinorTick','on','TickDir','in','TickLength',[0.005 0.005],'XGrid','on','YGrid','on');
set(gca,'FontSize',size1,'TickLabelInterpreter','latex');

%add text
labels1 = string(round(TASK2_mase,2));
text(categorical([0,7,30,90]),TASK2_mase,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',size1-8,Interpreter='latex',Color='w')

ax = gca;
ax.YAxis(1).Limits = [150 200];
ax.YAxis(1).Color = '#4698c6';
ax.YAxis(2).Limits = [0.5 1.12];
ax.YAxis(2).Color = '#e8864f';

% %export figure
exportgraphics(t,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month19\Fire\Figure_FireBurntAreaResults.pdf'),"ContentType","image")

close("all")

%% Fire duration prediction over different horizons
clear
clc
%read data
TASK3_rmse = [6.056985184	6.307973973	6.454590846	6.337901246];
TASK3_mase = [0.492857316	0.528951788	0.522964311	0.524003877];

size1 = 40;
figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout('flow','TileSpacing','none','Padding','tight','units','normalized','outerposition',[0 0 1 1]);
nexttile,

yyaxis left
b = bar(categorical([0,7,30,90]),TASK3_rmse,'EdgeColor','none','FaceAlpha',0.6);
ylabel('RMSE (days)','FontSize',size1,'Interpreter','latex','Color','#4698c6')
xlabel('Prediction horizon (days)','FontSize',size1,'Interpreter','latex')

%add text
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = strcat(string(round(b(1).YData,2)),'$days$');
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',size1-8,Interpreter='latex',Color='#4698c6')

hold on

yyaxis right
p = plot(categorical([0,7,30,90]),TASK3_mase,'LineStyle','-','LineWidth', 3, 'Marker','o', 'MarkerSize',10, 'MarkerFaceColor','w', 'Color','#e8864f');
ylabel('MASE','FontSize',size1,'Interpreter','latex');
set(gca,'Box','on','XTickLabelRotation',0,'FontSize',size1-5,'TickLabelInterpreter','latex','LineWidth',0.5,'YMinorGrid','on','YMinorTick','on','TickDir','in','TickLength',[0.005 0.005],'XGrid','on','YGrid','on');
set(gca,'FontSize',size1,'TickLabelInterpreter','latex');

%add text
labels1 = string(round(TASK3_mase,2));
text(categorical([0,7,30,90]),TASK3_mase,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',size1-8,Interpreter='latex',Color='w')

ax = gca;
ax.YAxis(1).Limits = [3 7];
ax.YAxis(1).Color = '#4698c6';
ax.YAxis(2).Limits = [0 0.75];
ax.YAxis(2).Color = '#e8864f';

%export figure
exportgraphics(gca,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month19\Fire\Figure_FireDurationResults.pdf'),"ContentType","image")

close("all")

%% Plot size and duration of wildfires in geobubble
size1 = 40;

path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\STClusteredFires_v2.csv');
opts = detectImportOptions(path1);
data = readtable(path1,opts);
data.Duration = categorical(data.Duration);

% f = figure('units','normalized','outerposition',[0 0 1 1]);
% 
% 
% gb = geobubble(data,'avgLatitude','avgLongitude', ...
%     'SizeVariable','totalArea_km2', 'basemap', 'topographic');
% gb.SourceTable.Durationlevels = discretize(double(data.Duration),[0 5 15 30 60 73],...
%                                  'categorical', {']0-5] days', '[6-15] days', '[16-30] days', '[31-60] days', '> 60 days'});
% gb.ColorVariable = 'Durationlevels';
% gb.SizeLegendTitle = 'Burnt Area (km^2)';
% gb.ColorLegendTitle = 'Duration';
% set(gca,'FontSize',size1-10);
% 
% gb.LongitudeAxis.FontSize = size1-10;
% gb.LongitudeAxis.TickLabelInterpreter = 'latex';
% gb.LongitudeAxis.TickLabels = strrep(get(gx.LongitudeAxis,'TickLabels'),'°','$^{\circ}$');
% gb.LatitudeAxis.FontSize = size1-10;
% gb.LatitudeAxis.TickLabelInterpreter = 'latex';
% gb.LatitudeAxis.TickLabels = strrep(get(gb.LatitudeAxis,'TickLabels'),'°','$^{\circ}$');

clc
close('all')

figure('units','normalized','OuterPosition',[0 0 0.740625 1]);
t = tiledlayout(1,1,"TileSpacing","tight","Padding","tight");
gx = geoaxes;
geobasemap(gx, 'topographic'); % You can choose other basemap options like 'streets', 'topographic', etc.

cats = {']0-5] days', '[6-15] days', '[16-30] days', '[31-60] days', '$>$ 60 days'};
Durationlevels = discretize(double(data.Duration),[0 5 15 30 60 73],...
                                 'categorical', cats);
hold on
for i = 1:length(cats)
    bubblechart(data.avgLatitude(Durationlevels==cats{i}),data.avgLongitude(Durationlevels==cats{i}),data.totalArea_km2(Durationlevels==cats{i}),"MarkerEdgeColor","none");
end
hold off
set(gca,'FontSize',size1-10);

bubblesize([5 20]);
blgd = bubblelegend('Burnt Area ($km^2$)');
set(blgd,'Interpreter','latex','FontSize',size1-20,'Location','eastoutside');
blgd.LimitLabels = {'0.001';'3966'};

lgd = legend(cats);
title(lgd,'Duration');
set(lgd,'Interpreter','latex','FontSize',size1-20,'Location','eastoutside');

lgd.Position = [0.742146375734849 0.696605382950981 0.15425106520849 0.215857817322612];
blgd.Position = [0.728924855859837 0.540763956023519 0.167472585083502 0.134518580238916];

gx.ZoomLevel = 5.020547898922111;
gx.MapCenter = [64.18470805276995,-151.8460677603565];

gx.LongitudeAxis.FontSize = size1-10;
gx.LatitudeAxis.FontSize = size1-10;

gx.LatitudeAxis.TickLabelInterpreter = 'latex';
gx.LatitudeAxis.TickLabels = strrep(get(gx.LatitudeAxis,'TickLabels'),'°','$^{\circ}$');
gx.LongitudeAxis.TickLabelInterpreter = 'latex';
gx.LongitudeAxis.TickLabels = strrep(get(gx.LongitudeAxis,'TickLabels'),'°','$^{\circ}$');
gx.LatitudeLabel.Interpreter = 'latex';
gx.LongitudeLabel.Interpreter = 'latex';

% exportgraphics(gcf,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month19\Fire\Figure_FireClusters.pdf'),"ContentType","image")

ExportGEOTIFF(gca, 'Figure_FireClusters')
%% Figure wildfires and conditioning features maps

clc
clear
size1 = 40;
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\FireDetection_dataset_0Days.csv');
opts = detectImportOptions(path1);
MODIS = readtable(path1,opts);%,"RowTimes",'acq_date');
fireData = (MODIS);
fireData = fireData(fireData.isFire == 1,:);

% NDVI
fireData.NDVI = fireData.NDVI*0.0001;

f = figure('units','normalized','outerposition',[0 0 1 1]);
gx = geoaxes("Basemap","streets-light");
s = geoscatter(fireData,"avgLatitude","avgLongitude","filled");
s.SizeData = 20;
s.ColorVariable = "NDVI";
c = colorbar;
c.Label.String = "NDVI";
c.Label.Interpreter = 'latex';
set(gca,'FontSize',size1-10)
colormap turbo
gx.LatitudeAxis.TickLabelInterpreter  = 'Latex';
gx.LongitudeAxis.TickLabelInterpreter = 'Latex';
gx.LatitudeAxis.TickLabels = strrep(gx.LatitudeAxis.TickLabels, '°', '$^{\circ}$');
gx.LongitudeAxis.TickLabels = strrep(gx.LongitudeAxis.TickLabels, '°', '$^{\circ}$');
gx.LatitudeLabel.Interpreter = 'latex';
gx.LongitudeLabel.Interpreter = 'latex';
exportgraphics(gca,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month19\Fire\Figure_FiresNDVI.pdf'),"ContentType","image")
close('all')

% tmax

f = figure('units','normalized','outerposition',[0 0 1 1]);
gx = geoaxes("Basemap","streets-light");
s = geoscatter(fireData,"avgLatitude","avgLongitude","filled");
s.SizeData = 20;
s.ColorVariable = "tmax";
c = colorbar;
c.Label.String = "Maximum air temperature at 2m ($^{\circ} C$)";
c.Label.Interpreter = 'latex';
set(gca,'FontSize',size1-10)
colormap turbo
gx.LatitudeAxis.TickLabelInterpreter  = 'Latex';
gx.LongitudeAxis.TickLabelInterpreter = 'Latex';
gx.LatitudeAxis.TickLabels = strrep(gx.LatitudeAxis.TickLabels, '°', '$^{\circ}$');
gx.LongitudeAxis.TickLabels = strrep(gx.LongitudeAxis.TickLabels, '°', '$^{\circ}$');
gx.LatitudeLabel.Interpreter = 'latex';
gx.LongitudeLabel.Interpreter = 'latex';
exportgraphics(gca,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month19\Fire\Figure_Firestmax.pdf'),"ContentType","image")
close('all')

% dewpoint
fireData.dewpoint_temperature_2m = fireData.dewpoint_temperature_2m - 273.15;

f = figure('units','normalized','outerposition',[0 0 1 1]);
gx = geoaxes("Basemap","streets-light");
s = geoscatter(fireData,"avgLatitude","avgLongitude","filled");
s.SizeData = 20;
s.ColorVariable = "dewpoint_temperature_2m";
c = colorbar;
c.Label.String = "Dewpoint temperature at 2m ($^{\circ} C$)";
c.Label.Interpreter = 'latex';
set(gca,'FontSize',size1-10)
colormap turbo
gx.LatitudeAxis.TickLabelInterpreter  = 'Latex';
gx.LongitudeAxis.TickLabelInterpreter = 'Latex';
gx.LatitudeAxis.TickLabels = strrep(gx.LatitudeAxis.TickLabels, '°', '$^{\circ}$');
gx.LongitudeAxis.TickLabels = strrep(gx.LongitudeAxis.TickLabels, '°', '$^{\circ}$');
gx.LatitudeLabel.Interpreter = 'latex';
gx.LongitudeLabel.Interpreter = 'latex';
exportgraphics(gca,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month19\Fire\Figure_Firesdewpoint.pdf'),"ContentType","image")
close('all')

% prcp

f = figure('units','normalized','outerposition',[0 0 1 1]);
gx = geoaxes("Basemap","streets-light");
s = geoscatter(fireData,"avgLatitude","avgLongitude","filled");
s.SizeData = 20;
s.ColorVariable = "prcp";
c = colorbar;
c.Label.String = "Precipitation (mm)";
c.Label.Interpreter = 'latex';
set(gca,'FontSize',size1-10)
colormap turbo
gx.LatitudeAxis.TickLabelInterpreter  = 'Latex';
gx.LongitudeAxis.TickLabelInterpreter = 'Latex';
gx.LatitudeAxis.TickLabels = strrep(gx.LatitudeAxis.TickLabels, '°', '$^{\circ}$');
gx.LongitudeAxis.TickLabels = strrep(gx.LongitudeAxis.TickLabels, '°', '$^{\circ}$');
gx.LatitudeLabel.Interpreter = 'latex';
gx.LongitudeLabel.Interpreter = 'latex';
exportgraphics(gca,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month19\Fire\Figure_Firesprcp.pdf'),"ContentType","image")
close('all')

%% Class distribution of select features: fire event
clear
clc
i=0;
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',['FireDetection_dataset_' num2str(i) 'Days.csv']);
opts = detectImportOptions(path1);
data = readtable(path1,opts);
data.NDVI = data.NDVI*0.0001;
data.dewpoint_temperature_2m = data.dewpoint_temperature_2m - 273.15;
data.sur_refl_b07 = data.sur_refl_b07* 0.0001;
data.Albedo_WSA_shortwave = data.Albedo_WSA_shortwave* 0.0001;

features = {'NDVI', 'tmax', 'sur_refl_b07', 'prcp', 'v_component_of_wind_10m', 'Albedo_WSA_shortwave'};
size1 = 40;

figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout('flow','TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);

skwDiff = zeros(length(features),1);
ovlp = zeros(length(features),1);

for j = 1:length(features)
    nexttile,hold on

    feature = features{j};

    % Calculate Density distribution for each class
    xRange = linspace(min(data.(feature)), max(data.(feature)), 1000);
    %class 0
    [f0, xi] = ksdensity(data.(feature)(data.isFire==0),xRange); % Estimate the density
    normalized_f0 = f0 / trapz(xi, f0);
    area(xi, normalized_f0, 'DisplayName','No Fire', FaceAlpha = 0.5, EdgeColor='none'); % Plot the density
    %class 1
    [f1, xi] = ksdensity(data.(feature)(data.isFire==1),xRange); % Estimate the density
    normalized_f1 = f1 / trapz(xi, f1);
    area(xi, normalized_f1, 'DisplayName','Fire', FaceAlpha = 0.5, EdgeColor='none'); % Plot the density
    % ylabel('Density','FontSize',size1,'Interpreter','latex');
    set(gca,'Box','on','XTickLabelRotation',0,'FontSize',size1-15,'TickLabelInterpreter','latex','LineWidth',0.5,'YMinorGrid','on','YMinorTick','on','TickDir','in','TickLength',[0.005 0.005]);
    set(gca,'FontSize',size1-15,'TickLabelInterpreter','latex');

    % Calculate Skewness for each class
    skewnessClass0 = skewness(data.(feature)(data.isFire==0));
    skewnessClass1 = skewness(data.(feature)(data.isFire==1));
    skwDiff(j) = abs(skewnessClass1-skewnessClass0);

    % Calculate Overlap Coefficient (simple method)
    minF = min([normalized_f0; normalized_f1]);
    ovlp(j) = trapz(xRange, minF);
    switch feature
        case 'prcp'
            title(strcat('PRCP (mm)',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);
        case 'NDVI'
            title(strcat('NDVI',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);
        case 'tmax'
            title(strcat('Tmax ($^{\circ} C$)',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);
        case 'sur_refl_b07'
            title(strcat('MIR',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);
        case 'v_component_of_wind_10m'
            title(strcat('Wind\_v (m/s)',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);
        case 'Albedo_WSA_shortwave'
            title(strcat('WSA\_s',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);            
    end
end
ylabel(t, 'Density', 'FontSize', size1, 'Interpreter', 'latex');  % Shared Y-axis label for all plots

legend1 = legend(gca,'show');
legend1.Layout.Tile = 'North'; 
set(legend1,'Orientation','horizontal','Interpreter','latex','FontSize',size1-15);
exportgraphics(gcf,'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month19\Fire\Figure_CDselectFeatures_Task1.pdf',"ContentType","image")
%% Class distribution of select features: fire event V4
clear
clc
i=0;
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',['FireDetection_datasetV4_' num2str(i) 'Days.csv']);
opts = detectImportOptions(path1);
data = readtable(path1,opts); %LC1, Lat, Slope, Aspect, and Snowfall
data.NDVI = data.NDVI*0.0001;
data.dewpoint_temperature_2m = data.dewpoint_temperature_2m - 273.15;
data.sur_refl_b07 = data.sur_refl_b07* 0.0001;
data.Albedo_WSA_shortwave = data.Albedo_WSA_shortwave* 0.0001;
data.LST_Day_1km = data.LST_Day_1km*0.02 - 273.15;

features = {'LC_Type1', 'avgLatitude', 'Slope', 'Aspect', 'snowfall_sum', 'NDVI', 'LST_Day_1km'};
size1 = 30;

figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout(2,4,'TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);

skwDiff = zeros(length(features),1);
ovlp = zeros(length(features),1);

for j = 1:length(features)
    nexttile,hold on

    feature = features{j};

    % Calculate Density distribution for each class
    xRange = linspace(min(data.(feature)), max(data.(feature)), 1000);
    %class 0
    [f0, xi] = ksdensity(data.(feature)(data.isFire==0),xRange); % Estimate the density
    normalized_f0 = f0 / trapz(xi, f0);
    area(xi, normalized_f0, 'DisplayName','No Fire', FaceAlpha = 0.5, EdgeColor='none'); % Plot the density
    %class 1
    [f1, xi] = ksdensity(data.(feature)(data.isFire==1),xRange); % Estimate the density
    normalized_f1 = f1 / trapz(xi, f1);
    area(xi, normalized_f1, 'DisplayName','Fire', FaceAlpha = 0.5, EdgeColor='none'); % Plot the density
    % ylabel('Density','FontSize',size1,'Interpreter','latex');
    set(gca,'Box','on','XTickLabelRotation',0,'FontSize',size1-15,'TickLabelInterpreter','latex','LineWidth',0.5,'YMinorGrid','on','YMinorTick','on','TickDir','in','TickLength',[0.005 0.005]);
    set(gca,'FontSize',size1-15,'TickLabelInterpreter','latex');

    % Calculate Skewness for each class
    skewnessClass0 = skewness(data.(feature)(data.isFire==0));
    skewnessClass1 = skewness(data.(feature)(data.isFire==1));
    skwDiff(j) = abs(skewnessClass1-skewnessClass0);

    % Calculate Overlap Coefficient (simple method)
    minF = min([normalized_f0; normalized_f1]);
    ovlp(j) = trapz(xRange, minF);
    switch feature
        case 'prcp'
            title(strcat('PRCP (mm)',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);
        case 'NDVI'
            title(strcat('NDVI',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);
        case 'tmax'
            title(strcat('Tmax ($^{\circ} C$)',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);
        case 'sur_refl_b07'
            title(strcat('MIR',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);
        case 'v_component_of_wind_10m'
            title(strcat('Wind\_v (m/s)',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);
        case 'Albedo_WSA_shortwave'
            title(strcat('WSA\_s',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);            
        case 'snow_cover'
            title(strcat('Snow\_c',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);            
        case 'LST_Day_1km'
            title(strcat('LST\_d ($^{\circ} C$)',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);            
        case 'snowfall\_sum'
            title(strcat('Snowfall (m = water equ.)',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);            
        case 'LC1'
            title(strcat('LC1',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);            
        case 'avgLatitude'
            title(strcat('Lat',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);            
        case 'Slope'
            title(strcat('Slope',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);            
        case 'Aspect'
            title(strcat('Aspect',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);            
            	
    end
end
ylabel(t, 'Density', 'FontSize', size1, 'Interpreter', 'latex');  % Shared Y-axis label for all plots

legend1 = legend(gca,'show');
legend1.Layout.Tile = 'North'; 
set(legend1,'Orientation','horizontal','Interpreter','latex','FontSize',size1-15);
% exportgraphics(gcf,'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month19\Fire\Figure_CDselectFeatures_V4Task1.pdf',"ContentType","image")

%% Class distribution of select features: burnt area

clear
clc

i=0;
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',['FireBurntArea_dataset_' num2str(i) 'Days.csv']);
opts = detectImportOptions(path1);
data = readtable(path1,opts);
data.NDVI = data.NDVI*0.0001;
data.dewpoint_temperature_2m = data.dewpoint_temperature_2m - 273.15;

% Discretize
[Y,E] = discretize(data.totalArea_km2,[0 100 inf]);

features = {'v_component_of_wind_10m'};%'NDVI', 'tmax', 'dewpoint_temperature_2m', 'prcp', 
size1 = 40;

figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout('flow','TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);

skwDiff = zeros(length(features),1);
ovlp = zeros(length(features),1);

for j = 1:length(features)
    nexttile,hold on

    feature = features{j};
    
    % Calculate Density distribution for each class
    xRange = linspace(min(data.(feature)), max(data.(feature)), 1000);
    %class 0
    [f0, xi] = ksdensity(data.(feature)(Y==1),xRange); % Estimate the density
    normalized_f0 = f0 / trapz(xi, f0);
    area(xi, normalized_f0, 'DisplayName','$]0-100[km^2$', FaceAlpha = 0.5, EdgeColor='none'); % Plot the density
    %class 1
    [f1, xi] = ksdensity(data.(feature)(Y==2),xRange); % Estimate the density
    normalized_f1 = f1 / trapz(xi, f1);
    area(xi, normalized_f1, 'DisplayName','$\ge100km^2$', FaceAlpha = 0.5, EdgeColor='none'); % Plot the density
    ylabel('Density','FontSize',size1,'Interpreter','latex');
    set(gca,'Box','on','XTickLabelRotation',0,'FontSize',size1-15,'TickLabelInterpreter','latex','LineWidth',0.5,'YMinorGrid','on','YMinorTick','on','TickDir','in','TickLength',[0.005 0.005]);
    set(gca,'FontSize',size1-15,'TickLabelInterpreter','latex');

    % Calculate Skewness for each class
    skewnessClass0 = skewness(data.(feature)(Y==1));
    skewnessClass1 = skewness(data.(feature)(Y==2));
    skwDiff(j) = abs(skewnessClass1-skewnessClass0);

    % Calculate Overlap Coefficient (simple method)
    minF = min([normalized_f0; normalized_f1]);
    ovlp(j) = trapz(xRange, minF);
    
    if strcmp(feature, 'prcp')
        title(strcat('Precipitation (mm)',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);
    else
        if strcmp(feature, 'NDVI')
            title(strcat('NDVI',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);
        else
            if strcmp(feature, 'tmax')
                title(strcat('Max. temp. ($^{\circ} C$)',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);
            else
                if strcmp(feature, 'dewpoint_temperature_2m')
                    title(strcat('Dewpoint temp. ($^{\circ} C$)',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);
                end
            end
        end
    end

end
legend1 = legend(gca,'show');
legend1.Layout.Tile = 'North'; % <-- place legend east of tiles
set(legend1,'Orientation','horizontal','Interpreter','latex','FontSize',size1-15);
% exportgraphics(gcf,'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month19\Fire\Figure_CDselectFeatures_Task2.pdf',"ContentType","image")

%% Class distribution of select features: fire duration

clear
clc
i=0;
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',['FireDuration_dataset_' num2str(i) 'Days.csv']);
opts = detectImportOptions(path1);
data = readtable(path1,opts);

% Discretize
[Y,E] = discretize(data.Duration,[0 5 inf]);

features = {'NDVI', 'tmax', 'dewpoint_temperature_2m', 'prcp'};
size1 = 40;

figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout('flow','TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);

skwDiff = zeros(length(features),1);
ovlp = zeros(length(features),1);

for j = 1:length(features)
    nexttile,hold on

    feature = features{j};
    
    % Calculate Density distribution for each class
    xRange = linspace(min(data.(feature)), max(data.(feature)), 1000);
    %class 0
    [f0, xi] = ksdensity(data.(feature)(Y==1),xRange); % Estimate the density
    normalized_f0 = f0 / trapz(xi, f0);
    area(xi, normalized_f0, 'DisplayName','$[1-5[ days$', FaceAlpha = 0.5, EdgeColor='none'); % Plot the density
    %class 1
    [f1, xi] = ksdensity(data.(feature)(Y==2),xRange); % Estimate the density
    normalized_f1 = f1 / trapz(xi, f1);
    area(xi, normalized_f1, 'DisplayName','$\ge5 days$', FaceAlpha = 0.5, EdgeColor='none'); % Plot the density
    ylabel('Density','FontSize',size1,'Interpreter','latex');
    set(gca,'Box','on','XTickLabelRotation',0,'FontSize',size1-15,'TickLabelInterpreter','latex','LineWidth',0.5,'YMinorGrid','on','YMinorTick','on','TickDir','in','TickLength',[0.005 0.005]);
    set(gca,'FontSize',size1-15,'TickLabelInterpreter','latex');

    % Calculate Skewness for each class
    skewnessClass0 = skewness(data.(feature)(Y==1));
    skewnessClass1 = skewness(data.(feature)(Y==2));
    skwDiff(j) = abs(skewnessClass1-skewnessClass0);

    % Calculate Overlap Coefficient (simple method)
    minF = min([normalized_f0; normalized_f1]);
    ovlp(j) = trapz(xRange, minF);
    
    if strcmp(feature, 'prcp')
        title(strcat('Precipitation (mm)',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);
    else
        if strcmp(feature, 'NDVI')
            title(strcat('NDVI',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);
        else
            if strcmp(feature, 'tmax')
                title(strcat('Max. temp. ($^{\circ} C$)',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);
            else
                if strcmp(feature, 'dewpoint_temperature_2m')
                    title(strcat('Dewpoint temp. ($^{\circ} C$)',{', $\Delta s_{0-1}$='},num2str(round(abs(skwDiff(j)),2)),{', OC='},num2str(round(ovlp(j),2))),'Interpreter','latex','FontSize',size1-15);
                end
            end
        end
    end

end
legend1 = legend(gca,'show');
legend1.Layout.Tile = 'North'; % <-- place legend east of tiles
set(legend1,'Orientation','horizontal','Interpreter','latex','FontSize',size1-15);
exportgraphics(gcf,'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month19\Fire\Figure_CDselectFeatures_Task3.pdf',"ContentType","image")
%% (OLD) Correlation analysis among GA-selected wildfire factors - plot
close('all')
clear
i=0;
size1 = 35;
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',['FireBurntArea_dataset_' num2str(i) 'Days.csv']);
opts = detectImportOptions(path1);
data = readtable(path1,opts);
features1 = {'LC_Type1', 'avgLatitude', 'Slope', 'Aspect', 'H', 'surface_pressure', 'snowfall_sum', 'snowmelt_sum', 'snow_cover'};
features2 = {'LC_Type1', 'LC_Type3', 'LC_Type4', 'avgLongitude', 'EVI', 'sur_refl_b02', 'Albedo_BSA_nir', 'Albedo_BSA_shortwave', 'Albedo_WSA_nir', 'dayl', 'LST_Day_1km', 'LST_Night_1km', 'Y', 'M', 'D', 'H', 'v_component_of_wind_10m', 'surface_pressure', 'total_evaporation_sum', 'snowfall_sum', 'snowmelt_sum', 'snow_cover', 'temperature_2m'};
features3 = {'LC_Type1', 'LC_Type2', 'LC_Type3', 'LC_Type4', 'avgLongitude', 'avgLatitude', 'Slope', 'NDVI', 'sur_refl_b07', 'Albedo_WSA_nir', 'Albedo_WSA_shortwave', 'prcp', 'srad', 'tmax', 'Y', 'M', 'D', 'W', 'v_component_of_wind_10m', 'u_component_of_wind_10m', 'runoff_sum', 'snow_depth', 'dewpoint_temperature_2m', 'temperature_2m'};%{'LC_Type1', 'Elevation', 'NDVI', 'sur_refl_b01', 'sur_refl_b07', 'Albedo_BSA_nir', 'Albedo_WSA_shortwave', 'prcp', 'tmax', 'LST_Day_1km', 'M', 'S', 'v_component_of_wind_10m', 'surface_pressure', 'total_evaporation_sum', 'snow_depth', 'snow_cover'};
% featuresV4 = unique([features1 features2 features3]);%{'LC_Type1', 'LC_Type2', 'avgLongitude', 'avgLatitude', 'Elevation', 'Slope', 'sur_refl_b03', 'sur_refl_b07', 'Clear_day_cov', 'Albedo_BSA_nir', 'Albedo_BSA_shortwave', 'Albedo_WSA_nir', 'srad', 'dayl', 'tmax', 'tmin', 'M', 'W', 'H', 'u_component_of_wind_10m', 'surface_pressure', 'runoff_sum', 'snow_density'};

features_cmb = string(unique([categorical(features1)';categorical(features2)';categorical(features3)']));

% Dictionary of misspelled words and their corrections
corrections = {'avgLongitude','Lon';'avgLatitude','Lat';
    'LC_Type1','LC1'; 'LC_Type2','LC2'; 'LC_Type3','LC3'; 'LC_Type4','LC4';...
    'LST_Night_1km','LST\_n'; 'LST_Day_1km','LST\_d';...
    'sur_refl_b07','MIR';'sur_refl_b01','Red';'sur_refl_b02','NIR';'sur_refl_b03','Blue';...
    'Albedo_BSA_nir','BSA\_NIR';'Albedo_WSA_nir','WSA\_NIR';'Albedo_WSA_shortwave','WSA\_s';'Albedo_BSA_shortwave','BSA\_s';...
    'temperature_2m','Temp'; 'tmax', 'Tmax'; 'tmin', 'Tmin';...
    'v_component_of_wind_10m','Wind\_v';'u_component_of_wind_10m','Wind\_u';
    'surface_pressure', 'Pressure';...
    'prcp', 'PRCP'; 'Clear_day_cov','Clear\_d'; 'Clear_night_cov','Clear\_n'; ...
    'total_evaporation_sum', 'Evap';'runoff_sum','Runoff';...
    'snow_density', 'Snow\_de';'snow_depth', 'Snow\_d';'snow_cover', 'Snow\_c'; 'snowmelt_sum', 'Snow\_m'; 'snowfall_sum', 'Snowfall'};

% Number of corrections to perform
numCorrections = size(corrections, 1);
seltask = features2;

% all features at once:
[rho,pval] = corr(data{:,seltask});%featuresV4
figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout(1,1,'TileSpacing','tight','Padding','compact','units','normalized','outerposition',[0 0 1 1]);
nexttile
imagesc(rho); % plot the matrix
num_features = size(data{:,seltask},2);
% Loop through each correction and replace in the cell array
features = []; %#ok<*NASGU>
tmp = seltask;
for i = 1:numCorrections
    % Find indices of the misspelled word
    idx = strcmp(tmp, corrections{i, 1});

    % Replace misspelled words with the correct ones
    tmp(idx) = {corrections{i, 2}};
end
features = tmp;

set(gca, 'XTick', 1:num_features); % center x-axis ticks on bins
set(gca, 'YTick', 1:num_features); % center y-axis ticks on bins
set(gca, 'XTickLabel', features,'TickLabelInterpreter','latex','XTickLabelRotation',90,'FontSize',size1-10); % set x-axis labels
set(gca, 'YTickLabel', features,'TickLabelInterpreter','latex','FontSize',size1-10); % set y-axis labels
c = colorbar;
c.Label.String = "Correlation coefficient ($\rho$)";
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
set(gca,'FontSize',size1-10)
clim([-1 1])
load('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\ColorMapColors_corr.mat')
colormap(CustomColormap)
% Add correlation coefficients as text inside each cell
for i = 1:num_features
    for j = 1:num_features
        if isequal(i,j) || [~isequal(i,j) & abs(rho(i,j))>=0.649]
            if isequal(i,j)
                text(j, i, sprintf('%.0f', rho(i,j)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','Interpreter','latex','FontSize',size1-15, 'Color','w');
            else
                text(j, i, sprintf('%.2f', rho(i,j)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','Interpreter','latex','FontSize',size1-15, 'Color','w');
            end
        % else
        %     if ~isequal(i,j) & abs(rho(i,j))<0.649
        %         text(j, i, sprintf('%.2f', rho(i,j)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','Interpreter','latex','FontSize',size1-15, 'Color','k');
        %     end
        end
    end
end
set(gca, 'TickDir','none')
% exportgraphics(t,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\CorrPlot_GA.pdf'),"ContentType","image")
%% (NEW) Correlation analysis among GA-selected wildfire factors - plot
clc
close('all')
clear
i=0;
size1 = 35;

features1 = {'LC_Type1', 'avgLatitude', 'Slope', 'Aspect', 'H', 'surface_pressure', 'snowfall_sum', 'snowmelt_sum', 'snow_cover'};
features2 = {['LC_Type1', 'avgLongitude', 'EVI', 'Albedo_BSA_nir', 'dayl', 'LST_Night_1km', 'Y', 'M', 'D', 'v_component_of_wind_10m', 'surface_pressure', 'snowmelt_sum']};
features3 = {'LC_Type1', 'LC_Type2', 'LC_Type3', 'LC_Type4', 'avgLongitude', 'avgLatitude', 'Slope', 'NDVI', 'sur_refl_b07', 'Albedo_WSA_nir', 'Albedo_WSA_shortwave', 'prcp', 'srad', 'tmax', 'Y', 'M', 'D', 'W', 'v_component_of_wind_10m', 'u_component_of_wind_10m', 'runoff_sum', 'snow_depth', 'dewpoint_temperature_2m', 'temperature_2m'};%{'LC_Type1', 'Elevation', 'NDVI', 'sur_refl_b01', 'sur_refl_b07', 'Albedo_BSA_nir', 'Albedo_WSA_shortwave', 'prcp', 'tmax', 'LST_Day_1km', 'M', 'S', 'v_component_of_wind_10m', 'surface_pressure', 'total_evaporation_sum', 'snow_depth', 'snow_cover'};


% Dictionary of misspelled words and their corrections
corrections = {'avgLongitude','Lon';'avgLatitude','Lat';
    'LC_Type1','LC1'; 'LC_Type2','LC2'; 'LC_Type3','LC3'; 'LC_Type4','LC4';...
    'LST_Night_1km','LST\_n'; 'LST_Day_1km','LST\_d';...
    'sur_refl_b07','MIR';'sur_refl_b01','Red';'sur_refl_b02','NIR';'sur_refl_b03','Blue';...
    'Albedo_BSA_nir','BSA\_NIR';'Albedo_WSA_nir','WSA\_NIR';'Albedo_WSA_shortwave','WSA\_s';'Albedo_BSA_shortwave','BSA\_s';...
    'temperature_2m','Temp'; 'tmax', 'Tmax'; 'tmin', 'Tmin';...
    'v_component_of_wind_10m','Wind\_v';'u_component_of_wind_10m','Wind\_u';
    'surface_pressure', 'Pressure';...
    'prcp', 'PRCP'; 'Clear_day_cov','Clear\_d'; 'Clear_night_cov','Clear\_n'; ...
    'total_evaporation_sum', 'Evap';'runoff_sum','Runoff'; 'dewpoint_temperature_2m', 'Dewpoint';...
    'snow_density', 'Snow\_de';'snow_depth', 'Snow\_d';'snow_cover', 'Snow\_c'; 'snowmelt_sum', 'Snow\_m'; 'snowfall_sum', 'Snowfall'};

% Number of corrections to perform
numCorrections = size(corrections, 1);

seltask = {'occurrence', 'burnt', 'duration'};
L = 2; %change only this

path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',['correlation_matrix_', seltask{L} ,'.csv']);
opts = detectImportOptions(path1);
rho = readtable(path1,opts);


figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout(1,1,'TileSpacing','tight','Padding','compact','units','normalized','outerposition',[0 0 1 1]);
nexttile
imagesc(rho{:,2:end}); % plot the matrix
num_features = size(rho,2);
% Loop through each correction and replace in the cell array
features = []; %#ok<*NASGU>
tmp = rho{:,1};
for i = 1:numCorrections
    % Find indices of the misspelled word
    idx = strcmp(tmp, corrections{i, 1});

    % Replace misspelled words with the correct ones
    tmp(idx) = {corrections{i, 2}};
end
features = tmp;

set(gca, 'XTick', 1:num_features); % center x-axis ticks on bins
set(gca, 'YTick', 1:num_features); % center y-axis ticks on bins
set(gca, 'XTickLabel', features,'TickLabelInterpreter','latex','XTickLabelRotation',90,'FontSize',size1-10); % set x-axis labels
set(gca, 'YTickLabel', features,'TickLabelInterpreter','latex','FontSize',size1-10); % set y-axis labels
c = colorbar;
set(gca,'FontSize',size1-10)
clim([-1 1])
load('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\ColorMapColors_corr.mat')
colormap(CustomColormap)
colorbar( 'off' )
% Add correlation coefficients as text inside each cell
for i = 1:num_features-1
    rho{i,1};
    for j = 2:num_features
        if isequal(i,j-1) || [~isequal(i,j-1) & abs(rho{i,j})>=0.649]
            if isequal(i,j)
                text(j-1, i, sprintf('%.0f', rho{i,j}), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','Interpreter','latex','FontSize',size1-15, 'Color','w');
            else
                text(j-1, i, sprintf('%.2f', rho{i,j}), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','Interpreter','latex','FontSize',size1-15, 'Color','w');
            end
        % else
        %     if ~isequal(i,j) & abs(rho(i,j))<0.649
        %         text(j, i, sprintf('%.2f', rho(i,j)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','Interpreter','latex','FontSize',size1-15, 'Color','k');
        %     end
        end
    end
end
set(gca, 'TickDir','none')
set(gca, 'XAxisLocation','top','YAxisLocation','right')
box off
% exportgraphics(t,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\CorrPlot_GA_',seltask{L},'.pdf'),"ContentType","image")

%% Density plot of raw fire events in dataset%Geoplot and 
size1 = 40;
figure('units','normalized','outerposition',[0 0 1 1])
gx = geoaxes("Basemap","streets-light");
dp = geodensityplot(MODIS.latitude,MODIS.longitude);
dp.FaceColor = 'interp';
colormap turbo
set(gca,'FontSize',size1)
gx.LatitudeAxis.TickLabelInterpreter  = 'Latex';
gx.LongitudeAxis.TickLabelInterpreter = 'Latex';
gx.LatitudeAxis.TickLabels = strrep(gx.LatitudeAxis.TickLabels, '°', '$^{\circ}$');
gx.LongitudeAxis.TickLabels = strrep(gx.LongitudeAxis.TickLabels, '°', '$^{\circ}$');
gx.LatitudeLabel.Interpreter = 'latex';
gx.LongitudeLabel.Interpreter = 'latex';
% exportgraphics(gca,'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month19\Fire\Figure_density_rawfires.pdf',"ContentType","image")

%% Scatter plot of raw fire events in dataset
size1 = 40;
figure('units','normalized','outerposition',[0 0 1 1])
gx = geoaxes("Basemap","streets-light");
geoscatter(MODIS.latitude,MODIS.longitude, 'filled','o',SizeData=12,MarkerFaceAlpha=.6);
set(gca,'FontSize',size1)
gx.LatitudeAxis.TickLabelInterpreter  = 'Latex';
gx.LongitudeAxis.TickLabelInterpreter = 'Latex';
gx.LatitudeAxis.TickLabels = strrep(gx.LatitudeAxis.TickLabels, '°', '$^{\circ}$');
gx.LongitudeAxis.TickLabels = strrep(gx.LongitudeAxis.TickLabels, '°', '$^{\circ}$');
gx.LatitudeLabel.Interpreter = 'latex';
gx.LongitudeLabel.Interpreter = 'latex';
% exportgraphics(gca,'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month19\Fire\Figure_scatter_rawfires.pdf',"ContentType","image")
%% Plot of wildfire clustered events count
close('all')
clear
clc
i=0;
%fire detection counts per year
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2','Fires_task1_preprocessed_targets.csv');
opts = detectImportOptions(path1);
FireDetectiondata = readtable(path1,opts);
FireDetectiondata = FireDetectiondata(FireDetectiondata.isFire==1,:);
[FireDetectiondata.Year,~,~] = ymd(FireDetectiondata.DetectionDate_first)  ;
years = unique(FireDetectiondata.Year);
fireEventsPerYear = varfun(@sum, FireDetectiondata, 'InputVariables', 'isFire', ...
    'GroupingVariables', 'Year');

size1 = 40;
cmap = get(0, 'defaultaxescolororder');

figure('units','normalized','outerposition',[0 0 1 1])
b = bar((fireEventsPerYear.Year(2:end)),fireEventsPerYear.sum_isFire(2:end),'EdgeColor','none',FaceAlpha= 0.6,HandleVisibility='off');
ylabel('Wildfire occurrence','FontSize',size1,'Interpreter','latex')
% xlabel('Years','FontSize',size1,'Interpreter','latex');
%add text
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(round(b(1).YData,2));
text(xtips1,ytips1,labels1,'HorizontalAlignment','left',...
'VerticalAlignment','middle','FontSize',size1-8,Interpreter='latex',Color=cmap(1,:), Rotation=90)
set(gca,'Box','on','XTickLabelRotation',90,'FontSize',size1-5,'TickLabelInterpreter','latex','LineWidth',0.5,'YMinorGrid','on','YMinorTick','on','TickDir','in','TickLength',[0.005 0.005],'XGrid','on','YGrid','on');
set(gca,'FontSize',size1,'TickLabelInterpreter','latex');
hold on

% Plot the fit
% Get xdata from plot
xdata1 = get(b, 'xdata');
% Get ydata from plot
ydata1 = get(b, 'ydata');
% Make sure data are column vectors
xdata1 = xdata1(:);
ydata1 = ydata1(:);


% Remove NaN values and warn
nanMask1 = isnan(xdata1(:)) | isnan(ydata1(:));
if any(nanMask1)
    warning('GeneratedCode:IgnoringNaNs', ...
        'Data points with NaN coordinates will be ignored.');
    xdata1(nanMask1) = [];
    ydata1(nanMask1) = [];
end

% Find x values for plotting the fit based on xlim
axesLimits1 = xlim(gca);
xplot1 = linspace(axesLimits1(1), axesLimits1(2));

% Preallocate for "Show equations" coefficients
coeffs1 = cell(1,1);

% Find coefficients for polynomial (order = 3)
[fitResults1,S,mu] = polyfit(xdata1,ydata1,5);
R_squared = 1 - (S.normr/norm(ydata1 - mean(ydata1)))^2;

% Evaluate polynomial
yplot1 = polyval(fitResults1,xplot1,[],mu);

% Save type of fit for "Show equations"
fittypesArray1(1) = 2;

% Save coefficients for "Show Equation"
coeffs1{1} = fitResults1;

% % Plot the fit
% fitLine1 = plot(xplot1,yplot1,'XLimInclude','off',...
%     'Tag','linear',...
%     'Parent',gca,...
%     'MarkerSize',6,...
%     'LineWidth',2.4,...
%     'Color',cmap(2,:),'LineStyle',':');%,'FontSize',size1-8,'Interpreter','latex'
% legend1 =legend({strcat('5th order polynomial fit ($R^2=',num2str(round(R_squared,3)),'$)')});
% set(legend1,'Orientation','horizontal','Location','northwest','Interpreter','latex','FontSize',size1-8);
set(gca,'Ylim',[0 300],'XTick',2001:2022, 'XTickLabel',categorical(2001:2022))
yline(mean(fireEventsPerYear.sum_isFire(2:end)),'--',['Annual',newline,'Average'],'Interpreter','latex','LineWidth',2,'LineStyle','--','Color',cmap(2,:),...
    'LabelHorizontalAlignment','left','LabelVerticalAlignment','top','FontSize',size1-20,'HandleVisibility','off')

% exportgraphics(gcf,'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month19\Fire\Figure_counts.pdf',"ContentType","image")

%% Plot of burnt area of wildfire clusters - average
close('all')
clear
clc
i=0;
%fire detection counts per year
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2','Fires_task1_preprocessed_targets.csv');
opts = detectImportOptions(path1);
FireDetectiondata = readtable(path1,opts);
FireDetectiondata = FireDetectiondata(FireDetectiondata.isFire==1,:);
[FireDetectiondata.Year,~,~] = ymd(FireDetectiondata.DetectionDate_first);
years = unique(FireDetectiondata.Year);
fireEventsPerYear = varfun(@mean, FireDetectiondata, 'InputVariables', 'totalArea_km2', ...
    'GroupingVariables', 'Year');

size1 = 40;
cmap = get(0, 'defaultaxescolororder');

figure('units','normalized','outerposition',[0 0 1 1])
b = bar((fireEventsPerYear.Year(2:end)),fireEventsPerYear.mean_totalArea_km2(2:end),'EdgeColor','none',FaceAlpha= 0.6,HandleVisibility='off');
ylabel('Average Burnt Area ($km^2$)','FontSize',size1,'Interpreter','latex')
% xlabel('Years','FontSize',size1,'Interpreter','latex');
%add text
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(round(b(1).YData,1));
text(xtips1,ytips1,labels1,'HorizontalAlignment','left',...
'VerticalAlignment','middle','FontSize',size1-8,Interpreter='latex',Color=cmap(1,:), Rotation=90)
set(gca,'Box','on','XTickLabelRotation',90,'FontSize',size1-5,'TickLabelInterpreter','latex','LineWidth',0.5,'YMinorGrid','on','YMinorTick','on','TickDir','in','TickLength',[0.005 0.005],'XGrid','on','YGrid','on');
set(gca,'FontSize',size1,'TickLabelInterpreter','latex');
hold on

% Plot the fit
% Get xdata from plot
xdata1 = get(b, 'xdata');
% Get ydata from plot
ydata1 = get(b, 'ydata');
% Make sure data are column vectors
xdata1 = xdata1(:);
ydata1 = ydata1(:);


% Remove NaN values and warn
nanMask1 = isnan(xdata1(:)) | isnan(ydata1(:));
if any(nanMask1)
    warning('GeneratedCode:IgnoringNaNs', ...
        'Data points with NaN coordinates will be ignored.');
    xdata1(nanMask1) = [];
    ydata1(nanMask1) = [];
end

% Find x values for plotting the fit based on xlim
axesLimits1 = xlim(gca);
xplot1 = linspace(axesLimits1(1), axesLimits1(2));

% Preallocate for "Show equations" coefficients
coeffs1 = cell(1,1);

% Find coefficients for polynomial (order = 5)
[fitResults1,S,mu] = polyfit(xdata1,ydata1,6);
R_squared = 1 - (S.normr/norm(ydata1 - mean(ydata1)))^2;

% Evaluate polynomial
yplot1 = polyval(fitResults1,xplot1,[],mu);

% Save type of fit for "Show equations"
fittypesArray1(1) = 2;

% Save coefficients for "Show Equation"
coeffs1{1} = fitResults1;

% % Plot the fit
% fitLine1 = plot(xplot1,yplot1,'XLimInclude','off',...
%     'Tag','linear',...
%     'Parent',gca,...
%     'MarkerSize',6,...
%     'LineWidth',2.4,...
%     'Color',cmap(2,:),'LineStyle',':');%,'FontSize',size1-8,'Interpreter','latex'
% legend1 =legend({strcat('6th order polynomial fit ($R^2=',num2str(round(R_squared,3)),'$)')});
% set(legend1,'Orientation','horizontal','Location','northeast','Interpreter','latex','FontSize',size1-8);
set(gca,'Ylim',[0 200],'XTick',2001:2022, 'XTickLabel',categorical(2001:2022))
yline(mean(fireEventsPerYear.mean_totalArea_km2(2:end)),'--',['Annual',newline,'Average'],'Interpreter','latex','LineWidth',2,'LineStyle','--','Color',cmap(2,:),...
    'LabelHorizontalAlignment','left','LabelVerticalAlignment','top','FontSize',size1-20,'HandleVisibility','off')

% exportgraphics(gcf,'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month19\Fire\Figure_AvgArea.pdf',"ContentType","image")

%% Plot of durations of wildfire clusters - average
close('all')
clear
clc
i=0;
%fire detection counts per year
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2','Fires_task1_preprocessed_targets.csv');
opts = detectImportOptions(path1);
FireDetectiondata = readtable(path1,opts);
FireDetectiondata = FireDetectiondata(FireDetectiondata.isFire==1,:);
[FireDetectiondata.Year,~,~] = ymd(FireDetectiondata.DetectionDate_first);
years = unique(FireDetectiondata.Year);
fireEventsPerYear = varfun(@mean, FireDetectiondata, 'InputVariables', 'Duration', ...
    'GroupingVariables', 'Year');

size1 = 40;
cmap = get(0, 'defaultaxescolororder');

figure('units','normalized','outerposition',[0 0 1 1])
b = bar(fireEventsPerYear.Year(2:end),fireEventsPerYear.mean_Duration(2:end),'EdgeColor','none',FaceAlpha= 0.6,HandleVisibility='off');
ylabel(['Average wildfire duration (days)'],'FontSize',size1,'Interpreter','latex')
% xlabel('Years','FontSize',size1,'Interpreter','latex');
%add text
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(round(b(1).YData,2));
text(xtips1,ytips1,labels1,'HorizontalAlignment','left',...
'VerticalAlignment','middle','FontSize',size1-8,Interpreter='latex',Color=cmap(1,:), Rotation=90)
set(gca,'Box','on','XTickLabelRotation',90,'FontSize',size1-5,'TickLabelInterpreter','latex','LineWidth',0.5,'YMinorGrid','on','YMinorTick','on','TickDir','in','TickLength',[0.005 0.005],'XGrid','on','YGrid','on');
set(gca,'FontSize',size1,'TickLabelInterpreter','latex');
hold on

% Plot the fit
% Get xdata from plot
xdata1 = get(b, 'xdata');
% Get ydata from plot
ydata1 = get(b, 'ydata');
% Make sure data are column vectors
xdata1 = xdata1(:);
ydata1 = ydata1(:);


% Remove NaN values and warn
nanMask1 = isnan(xdata1(:)) | isnan(ydata1(:));
if any(nanMask1)
    warning('GeneratedCode:IgnoringNaNs', ...
        'Data points with NaN coordinates will be ignored.');
    xdata1(nanMask1) = [];
    ydata1(nanMask1) = [];
end

% Find x values for plotting the fit based on xlim
axesLimits1 = xlim(gca);
xplot1 = linspace(axesLimits1(1), axesLimits1(2));

% Preallocate for "Show equations" coefficients
coeffs1 = cell(1,1);

% Find coefficients for polynomial (order = 1)
ordr = 7;
[fitResults1,S,mu] = polyfit(xdata1,ydata1,ordr);
R_squared = 1 - (S.normr/norm(ydata1 - mean(ydata1)))^2;

% Evaluate polynomial
yplot1 = polyval(fitResults1,xplot1,[],mu);

% Save type of fit for "Show equations"
fittypesArray1(1) = 2;

% Save coefficients for "Show Equation"
coeffs1{1} = fitResults1;

% Plot the fit
% fitLine1 = plot(xplot1,yplot1,'XLimInclude','off',...
%     'Tag','linear',...
%     'Parent',gca,...
%     'MarkerSize',6,...
%     'LineWidth',2.4,...
%     'Color',cmap(2,:),'LineStyle',':');%,'FontSize',size1-8,'Interpreter','latex'
% legend1 =legend({strcat(num2str(ordr),'th order polynomial fit ($R^2=',num2str(round(R_squared,3)),'$)')});
% set(legend1,'Orientation','horizontal','Location','northeast','Interpreter','latex','FontSize',size1-8);
set(gca,'ylim',[0 20],'XTick',2001:2022, 'XTickLabel',categorical(2001:2022))

yline(mean(fireEventsPerYear.mean_Duration(2:end)),'--',['Annual',newline,'Average'],'Interpreter','latex','LineWidth',2,'LineStyle','--','Color',cmap(2,:),...
    'LabelHorizontalAlignment','left','LabelVerticalAlignment','top','FontSize',size1-20,'HandleVisibility','off')

% exportgraphics(gcf,'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month19\Fire\Figure_AvgDuration.pdf',"ContentType","image")

%% Fire detection 1st instance per year
clear
clc
size1 = 40;
cmap = get(0, 'defaultaxescolororder');

path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2','Fires_task1_preprocessed_targets.csv');
opts = detectImportOptions(path1);
FireDetectiondata = readtable(path1,opts);
FireDetectiondata = FireDetectiondata(FireDetectiondata.isFire==1,:);
[FireDetectiondata.Year,~,~] = ymd(FireDetectiondata.DetectionDate_first);

% Find unique years and the index of their first occurrence
[uniqueYears, firstIdx] = unique(FireDetectiondata.Year,'first');

% Index into the original timetable to get the first instance of each year
firstInstances = FireDetectiondata(firstIdx(2:end), :);
firstInstances(firstInstances.Year==2006,:) = FireDetectiondata(598,:); % remove the volcano -not wildfire- instance (January 2006)

figure('units','normalized','outerposition',[0 0 1 1]);
% plot([2000 2001],[day(FireDetectiondata{1,'DetectionDate_first'}, "dayofyear") day(firstInstances.DetectionDate_first(2), "dayofyear")],'o', 'Color', cmap(1,:), 'MarkerFaceColor', 'w', 'LineStyle', ':', 'MarkerSize', 12, 'LineWidth', 2','HandleVisibility','off')
hold on,
b = plot(firstInstances.Year,day(firstInstances.DetectionDate_first, "dayofyear"),'o', 'Color', cmap(1,:), 'MarkerFaceColor', cmap(1,:), 'LineStyle', '-', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName','Wildfire first occurrence');

set(gca,'ylim',[1 366],'XTick',2001:2022, 'XTickLabel',categorical(2001:2022),'xlim',[2001 2022])
ylabel('Days of year','FontSize',size1,'Interpreter','latex')
set(gca,'Box','on','XTickLabelRotation',90,'FontSize',size1-5,'TickLabelInterpreter','latex','LineWidth',0.5,'YMinorGrid','on','YMinorTick','on','TickDir','in','TickLength',[0.005 0.005],'XGrid','on','YGrid','on');
set(gca,'FontSize',size1,'TickLabelInterpreter','latex');

yregion(60,90,"DisplayName",'March','FaceColor',cmap(1,:),'FaceAlpha',0.2)
yregion(91,120,"DisplayName",'April','FaceColor',cmap(2,:),'FaceAlpha',0.2)
yregion(121,151,"DisplayName",'May','FaceColor',cmap(3,:),'FaceAlpha',0.2)
yregion(152,181,"DisplayName",'June','FaceColor',cmap(4,:),'FaceAlpha',0.2)

legend1 = legend(gca,'show');
set(legend1,'Orientation','horizontal','Location','northoutside','Interpreter','latex','FontSize',size1);

% Create textarrows
annotation(gcf,'textarrow',[0.702205882352941 0.6875],[0.238095238095239 0.304964539007093],'String',{'01-Jun-2017'},...
    'Interpreter','latex', 'FontSize',size1-10);
annotation(gcf,'textarrow',[0.826680672268907 0.829831932773108],[0.2006210739615 0.25531914893617],'String',{'28-Apr-2021'},...
    'Interpreter','latex','FontSize',size1-10);
annotation(gcf,'textarrow',[0.607142857142857 0.608193277310923],[0.279635258358663 0.311043566362716],'String',{'17-May-2015'},...
    'Interpreter','latex','FontSize',size1-10);

% exportgraphics(gcf,'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month22\Fire\Figure_1stOccurrence.pdf',"ContentType","image")

%% Monthly sums plot
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2','Fires_task1_preprocessed_targets.csv');
opts = detectImportOptions(path1);
FireDetectiondata = readtable(path1,opts);
FireDetectiondata = FireDetectiondata(FireDetectiondata.isFire==1,:);
TT = groupsummary(FireDetectiondata,"DetectionDate_first","month", "IncludeEmptyGroups",true);
TT.month_DetectionDate_first = datetime(string(TT.month_DetectionDate_first),'InputFormat','MMM-yyyy','Format','MMM-yyyy');

monthNames = {'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'};
monthsCategorical = categorical(monthNames, monthNames, 'Ordinal', true);

figure('units','normalized','outerposition',[0 0 1 1]);
b = plot(monthsCategorical',TT.GroupCount(find(year(TT.month_DetectionDate_first)==2001)),'o', 'Color', cmap(1,:), 'MarkerFaceColor', cmap(1,:), 'LineStyle', ':', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName','2001');
hold on
for i = 2002:2021
    tmp = TT.GroupCount(find(year(TT.month_DetectionDate_first)==i));
    plot(monthsCategorical,tmp, 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName',num2str(i));
end
linestyleorder("mixedmarkers")

% ylabel('Days of year','FontSize',size1,'Interpreter','latex')
% set(gca,'Box','on','XTickLabelRotation',90,'FontSize',size1-5,'TickLabelInterpreter','latex','LineWidth',0.5,'YMinorGrid','on','YMinorTick','on','TickDir','in','TickLength',[0.005 0.005],'XGrid','on','YGrid','on');
% set(gca,'FontSize',size1,'TickLabelInterpreter','latex');

%% (old) Error boxchart plots per 
clear
clc
dir_name = 'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\results\';
addpath(dir_name)

size1 = 40;
files = dir(fullfile(dir_name, '*.csv'));
trgts_all = {'Fire','Burnt','Duration'};
hrzs = {'hrz0','hrz7','hrz30','hrz90'};
hrzs_names = {'+0','+7','+30','+90'};
figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact','units','normalized','outerposition',[0 0 1 1]);

for j = 2:3
    nexttile
    trgt = trgts_all{j};
    filenames = {files.name}.';
    filenames = filenames(contains(filenames,trgt));

    for i = 1:4%all horizons
        % train_files = filenames(contains(filenames,hrzs{i}) & contains(filenames,'train'));
        test_files = filenames(contains(filenames,hrzs{i}) & contains(filenames,'test'));

        % eval(sprintf('train_data_%s = readtable(train_files{:},detectImportOptions(train_files{:}));',hrzs{i}))
        eval(sprintf('test_data_%s = readtable(test_files{:},detectImportOptions(test_files{:}));',hrzs{i}))

        eval(sprintf('ERROR_%s = (test_data_%s.test_targs-test_data_%s.test_preds);',hrzs{i},hrzs{i},hrzs{i}))
        switch j
            case 1
                eval(sprintf('cats_%s = discretize(test_data_%s.test_targs,[0 0.5 1],''categorical'', {''No fire'', ''Fire''});',hrzs{i},hrzs{i}))
            case 2
                eval(sprintf('cats_%s = discretize(test_data_%s.test_targs,[0 2 100 500 max(test_data_%s.test_targs)],''categorical'', {'']0-2]$km^2$'', '']2-100]$km^2$'', '']100-500]$km^2$'', ''$>$500$km^2$''});',hrzs{i},hrzs{i},hrzs{i}))
            case 3
                eval(sprintf('cats_%s = discretize(test_data_%s.test_targs,[0 2 7 14 30 max(test_data_%s.test_targs)],''categorical'', {'']0-2] days'', '']2-7] days'', '']7-14] days'', '']14-30] days'', ''$>$30 days''});',hrzs{i},hrzs{i},hrzs{i}))
        end
        eval(sprintf('tbl_%s = table(strings(size(cats_%s))+''%s'',cats_%s,ERROR_%s,''VariableNames'',{''Horizons'',''Categories'',''Errors''});',hrzs{i},hrzs{i},hrzs_names{i},hrzs{i},hrzs{i}))

    end
    tbl = [tbl_hrz0;tbl_hrz7;tbl_hrz30;tbl_hrz90];
    tbl.Horizons = categorical(tbl.Horizons,hrzs_names);
             %tbl.Horizons
    boxchart(tbl.Horizons,tbl.Errors,'GroupByColor',tbl.Categories,'JitterOutliers','on','MarkerStyle','none','MarkerSize',size1-20)
    set(gca,'Box','on','FontSize',size1-10,'TickLabelInterpreter','latex','LineWidth',0.5,'YMinorGrid','on','YMinorTick','on','TickDir','in','TickLength',[0.005 0.005],'YGrid','on');
    
    switch j
        case 1
            ylabel('Error','FontSize',size1-10,'Interpreter','latex')
            legend1 = legend(gca,'show');
            set(legend1,'Orientation','horizontal','Location','northoutside','Interpreter','latex','FontSize',size1-20);
            title(legend1,'Wildfire Occurrence')
        case 2
            ylabel('Error ($km^2$)','FontSize',size1-10,'Interpreter','latex')
            set(gca,'ylim',[-560 2000])
            legend1 = legend(gca,'show');
            set(legend1,'Orientation','horizontal','Location','northoutside','Interpreter','latex','FontSize',size1-20);
            title(legend1,'Total burnt area')
        case 3
            ylabel('Error (days)','FontSize',size1-10,'Interpreter','latex')
            set(gca,'ylim',[-10 32])
            legend1 = legend(gca,'show');
            set(legend1,'Orientation','horizontal','Location','northoutside','Interpreter','latex','FontSize',size1-20);
            title(legend1,'Wildfire duration')
    end
end
xlabel(t,'Prediction horizons (days)','FontSize',size1-10,'Interpreter','latex')

%edit location of the second legend before exporting!
% exportgraphics(gcf,fullfile(dir_name, strcat('ErrorsBoxplot_Task2_3.pdf')),"ContentType","image")

%% (new) bar charts with 95% confidence intervals for prediction errors
clear
clc
dir_name = 'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\results\';
addpath(dir_name)

size1 = 40;
files = dir(fullfile(dir_name, '*.csv'));
trgts_all = {'Fire','Burnt','Duration'};
hrzs = {'hrz0','hrz7','hrz30','hrz90'};
hrzs_names = {'+0','+7','+30','+90'};

% Define a custom color map
colors = lines(5); % Use 'lines' colormap

figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact','units','normalized','outerposition',[0 0 1 1]);

for j = 2:3
    nexttile
    trgt = trgts_all{j};
    filenames = {files.name}.';
    filenames = filenames(contains(filenames,trgt));
    
    means = zeros(4, 4); % Assuming 4 horizons and 4 categories
    ci_95 = zeros(4, 4); % Confidence intervals
    
    if j == 2
        all_categories = {']0-2]$km^2$', ']2-100]$km^2$', ']100-500]$km^2$', '$>$500$km^2$'};
        bin_edges = [0 2 100 500 Inf]; % Correct bin edges for burnt area
    else
        all_categories = {']0-2] days', ']2-7] days', ']7-14] days', ']14-30] days', '$>$30 days'};
        bin_edges = [0 2 7 14 30 Inf]; % Correct bin edges for duration
    end

    for i = 1:4 % All horizons
        test_files = filenames(contains(filenames,hrzs{i}) & contains(filenames,'test'));
        eval(sprintf('test_data_%s = readtable(test_files{:},detectImportOptions(test_files{:}));',hrzs{i}));
        eval(sprintf('ERROR_%s = (test_data_%s.test_targs-test_data_%s.test_preds);',hrzs{i},hrzs{i},hrzs{i}));
        
        % Using fixed bin edges that align with your categories
        eval(sprintf('cats_%s = discretize(test_data_%s.test_targs, bin_edges, ''categorical'', all_categories);',hrzs{i},hrzs{i}));

        for k = 1:length(all_categories)
            eval(sprintf('cat_errors = ERROR_%s(cats_%s == all_categories{k});', hrzs{i}, hrzs{i}));
            means(i, k) = mean(cat_errors);
            std_err = std(cat_errors) / sqrt(length(cat_errors)); % Standard error
            t_value = tinv(0.975, length(cat_errors) - 1); % T-value for 95% CI
            ci_95(i, k) = t_value * std_err;
        end
    end
    
    % Create grouped bar chart
    b = bar(1:4, means, 'grouped', 'BarWidth', 0.8, 'EdgeColor','none','FaceAlpha',0.5);
    hold on;
    
    % Color the bars
    for k = 1:length(all_categories)
        b(k).FaceColor = colors(k,:);
    end
    
    % Add error bars to each bar group
    ngroups = size(means, 1);
    nbars = size(means, 2);
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    for i = 1:nbars
        % Align error bars with bars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, means(:, i), ci_95(:, i), 'dk', 'linestyle', 'none', 'CapSize', 10, 'LineWidth', 1.5, 'MarkerFaceColor','white',MarkerEdgeColor='auto');
    end

    % Enhance the grid and axis appearance
    set(gca,'xticklabels', categorical(hrzs_names),'TickLabelInterpreter','latex', 'Box','on','FontSize',size1-10,'TickLabelInterpreter','latex','LineWidth',1,...
        'YMinorGrid','on','XMinorTick','on','TickDir','in','TickLength',[0.01 0.01],...
        'YGrid','off','XGrid','off', 'GridAlpha', 0.6); % Lighter grid lines with 'GridAlpha'

    switch j
        case 2
            ylabel('Error ($km^2$)','FontSize',size1-10,'Interpreter','latex')
            legend1 = legend(all_categories,'Orientation','horizontal','Location','northoutside','Interpreter','latex','FontSize',size1-20);
            title(legend1,'Total burnt area')
        case 3
            ylabel('Error (days)','FontSize',size1-10,'Interpreter','latex')
            legend1 = legend(all_categories,'Orientation','horizontal','Location','northoutside','Interpreter','latex','FontSize',size1-20);
            title(legend1,'Wildfire duration')
    end
end

% Add global xlabel
xlabel(t,'Prediction horizons (days)','FontSize',size1-10,'Interpreter','latex')

% Export % make sure to align the legends first
% exportgraphics(t,fullfile(dir_name, 'Errors95CI_Task2_3.pdf')) % Use vector for best quality

%% Function
function tmp = get_astronomical_season_str(dates)
    sns = ["Spring", "Summer", "Fall", "Winter"];
    % Extract the year from the dates
    yrs = year(dates);
    
    % Create datetime arrays for the equinoxes and solstices
    spring_equinox = datetime(yrs, 3, 20);
    summer_solstice = datetime(yrs, 6, 21);
    fall_equinox = datetime(yrs, 9, 22);
    winter_solstice = datetime(yrs, 12, 21);

    % Preallocate the seasons array
    seasons = zeros(size(dates));

    % Determine the season for each date
    seasons(dates >= spring_equinox & dates < summer_solstice) = 0; % Spring
    seasons(dates >= summer_solstice & dates < fall_equinox) = 1; % Summer
    seasons(dates >= fall_equinox & dates < winter_solstice) = 2; % Fall
    seasons(dates >= winter_solstice | dates < spring_equinox) = 3; % Winter

    tmp = sns(seasons + 1);
end



function ExportGEOTIFF(gx, basename)
    try % Export figure as a high-resolution JPEG
    % basename = 'Figure_RawDataset';
    exportgraphics(gx, [basename '.jpg'], 'Resolution', 300)
    
    % Define geographic reference object for world file
    rasterSize = size(imread([basename '.jpg']));
    [latitudeLimits,longitudeLimits] = geolimits(gx);
    R = georefcells(latitudeLimits, longitudeLimits, rasterSize);
    
    % Write the world file with geographic information
    worldfilename = getworldfilename([basename '.jpg']);
    worldfilewrite(R, worldfilename);
    
    % Read the exported JPEG and world file
    A1 = imread([basename '.jpg']);
    R1 = worldfileread(worldfilename, 'geographic', size(A1));
    
    % Specify the CRS (WGS84, which is EPSG:4326) and other georeferencing information
    geoKeyDirectoryTag.GTModelTypeGeoKey = 2; % Geographic coordinate system
    geoKeyDirectoryTag.GTRasterTypeGeoKey = 1; % Raster pixel is area
    geoKeyDirectoryTag.GeographicTypeGeoKey = 4326; % WGS84 EPSG code
    
    % Write the raster data to a GeoTIFF file
    filename1 = [basename '.tif'];
    geotiffwrite(filename1, A1, R1, 'GeoKeyDirectoryTag', geoKeyDirectoryTag);
    
    % Notify completion
    fprintf('GeoTIFF file has been generated successfully.\n');
    catch
        warning('Error!')
    end
end