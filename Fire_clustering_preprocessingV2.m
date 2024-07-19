%% Validation dataset map part 1
clear
clc
size1 = 28;
cmap = get(0, 'defaultaxescolororder');
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork',"Fire_archive_AlaskaOnly.csv");
opts = detectImportOptions(path1);
MODIS = readtable(path1,opts);

figure('units','normalized','outerposition',[0.5 0 0.5 1])
gx = geoaxes("Basemap","topographic");
geoplot(MODIS.latitude,MODIS.longitude, '.','MarkerFaceColor','auto','MarkerSize',6);%cmap(1,:)
set(gca,'FontSize',size1)
gx.LatitudeAxis.TickLabelInterpreter  = 'Latex';
gx.LongitudeAxis.TickLabelInterpreter = 'Latex';
gx.LatitudeAxis.TickLabels = strrep(gx.LatitudeAxis.TickLabels, '°', '$^{\circ}$');
gx.LongitudeAxis.TickLabels = strrep(gx.LongitudeAxis.TickLabels, '°', '$^{\circ}$');
gx.LatitudeLabel.Interpreter = 'latex';
gx.LongitudeLabel.Interpreter = 'latex';
% exportgraphics(gca,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month19\Fire\Figure_Firesdewpoint.pdf'),"ContentType","image")
[latitudeLimits,longitudeLimits] = geolimits(gx);
zoomlevel = get(gca,'ZoomLevel');
mapcenter = get(gca,'MapCenter');

% exportgraphics(gca,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month22\Fire\Figure_RawDataset.pdf'),"ContentType","image",'BackgroundColor','none')
fprintf('DONE!')

%% Validation dataset map part 2
% Let's try mapping!
% Read in shapefile with perimiters:
% Read as geotable?
T = readgeotable('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\InterAgencyFirePerimeterHistory_All_Years_View\InterAgencyFirePerimeterHistory_All_Years_View.shp');

% Get array of fire years:
Fyears = T(:,9);
Fyears = str2double(table2array(Fyears));
T.FIRE_YEAR = Fyears;

T = T(T.FIRE_YEAR>1999 & T.FIRE_YEAR<2030,:); %remove missing data too!

% Plot:
clc
close('all')
size1 = 28;
figure('units','normalized','outerposition',[0.5 0 0.5 1],'visible','on')
gx = geoaxes("Basemap","topographic");
geoplot(T(startsWith(T.POO_RESP_I,'AK'),:));
geolimits([55,70],[-168,-123])
set(gca,'FontSize',size1)
gx.LatitudeAxis.TickLabelInterpreter  = 'Latex';
gx.LongitudeAxis.TickLabelInterpreter = 'Latex';
gx.LatitudeAxis.TickLabels = strrep(gx.LatitudeAxis.TickLabels, '°', '$^{\circ}$');
gx.LongitudeAxis.TickLabels = strrep(gx.LongitudeAxis.TickLabels, '°', '$^{\circ}$');
gx.LatitudeLabel.Interpreter = 'latex';
gx.LongitudeLabel.Interpreter = 'latex';
geolimits(latitudeLimits,longitudeLimits)
set(gca,'ZoomLevel',zoomlevel,'MapCenter',mapcenter)
% exportgraphics(gcf,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month22\Fire\Figure_ValidationDataset.pdf'),"ContentType","image",'BackgroundColor','none')

fprintf('DONE!')


%% Insel
% Plot:
clc
close('all')
size1 = 28;

path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork',"Fire_archive_AlaskaOnly.csv");
opts = detectImportOptions(path1);
MODIS = readtable(path1,opts);

years = year(MODIS.acq_date);
MODIS(years==max(T.FIRE_YEAR)+1,:) = [];
% MODIS = MODIS(years==2001,:);

figure('units','normalized','outerposition',[0.5 0 0.5 1],'visible','on')
gx = geoaxes("Basemap","topographic");
geoplot(T(startsWith(T.POO_RESP_I,'AK'),:),'DisplayName','Historical records');
hold on
geoplot(MODIS.latitude,MODIS.longitude, '.','MarkerFaceColor','auto','MarkerSize',6,'DisplayName','FIRMS dataset');%cmap(1,:)
geolimits([55,70],[-168,-123])
set(gca,'FontSize',size1)
gx.LatitudeAxis.TickLabelInterpreter  = 'Latex';
gx.LongitudeAxis.TickLabelInterpreter = 'Latex';
gx.LatitudeAxis.TickLabels = strrep(gx.LatitudeAxis.TickLabels, '°', '$^{\circ}$');
gx.LongitudeAxis.TickLabels = strrep(gx.LongitudeAxis.TickLabels, '°', '$^{\circ}$');
gx.LatitudeLabel.Interpreter = 'latex';
gx.LongitudeLabel.Interpreter = 'latex';
legend1 = legend(gca,'show');
set(legend1,'Interpreter','latex','FontSize',size1,'Orientation','horizontal','Location','north');

latitudeLimits = [65.7248   66.7749];
longitudeLimits = [-147.5650 -145.2604];
geolimits(latitudeLimits,longitudeLimits)
set(gca,'ZoomLevel',8.7609,'MapCenter',[66.2553 -146.4127])

% exportgraphics(gcf,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month22\Fire\Figure_ValidationInsel.pdf'),"ContentType","image")

fprintf('DONE!')

%% Validation of clustered events
% Plot:
clc
close('all')
clear
size1 = 28;

latitudeLimits = [65.7248   66.7749];
longitudeLimits = [-147.5650 -145.2604];

T = readgeotable('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\InterAgencyFirePerimeterHistory_All_Years_View\InterAgencyFirePerimeterHistory_All_Years_View.shp');

% Get array of fire years:
Fyears = T(:,9);
Fyears = str2double(table2array(Fyears));
T.FIRE_YEAR = Fyears;

T = T(T.FIRE_YEAR>1999 & T.FIRE_YEAR<2030,:); %remove missing data too!

path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork',"Task_FireDetection_noRS.csv");
opts = detectImportOptions(path1);
opts.VariableTypes(2:3) = repmat({'datetime'},[1, length(opts.VariableNames(2:3))]);%convert all content to datetime
opts = setvaropts(opts, "DetectionDate_first", "InputFormat", "yyyy-dd-MM");
opts = setvaropts(opts, "DetectionDate_last", "InputFormat", "yyyy-dd-MM");
Fires_posClass = readtable(path1,opts);
Fires_posClass(Fires_posClass.isFire==0,:) = [];
years = year(Fires_posClass.DetectionDate_first);
Fires_posClass(years==max(T.FIRE_YEAR)+1,:) = [];
% Fires_posClass = Fires_posClass(years==2005,:);
% T = T(T.FIRE_YEAR == 2005,:);

% Find indices where avgLatitude is within the latitude limits
keep_lat = find(Fires_posClass.avgLatitude >= latitudeLimits(1) & Fires_posClass.avgLatitude <= latitudeLimits(2));
% Find indices where avgLongitude is within the longitude limits
keep_lon = find(Fires_posClass.avgLongitude >= longitudeLimits(1) & Fires_posClass.avgLongitude <= longitudeLimits(2));
%rows of interest 
shared_rows = intersect(keep_lat, keep_lon);

figure('outerposition',[672          -6        1106        1087],'visible','on')
% t = tiledlayout(1,1,'TileSpacing','tight','Padding','compact','units','normalized','outerposition',[0 0 1 1]);
gx = geoaxes("Basemap","streets-light");
geoplot(T(startsWith(T.POO_RESP_I,'AK'),:),'DisplayName','Historical records');
hold on
geoplot(Fires_posClass.avgLatitude(shared_rows),Fires_posClass.avgLongitude(shared_rows), '.','MarkerFaceColor','auto','MarkerSize',15,'DisplayName','Temporal-Spectral clusters');%cmap(1,:)
set(gca,'FontSize',size1)
gx.LatitudeAxis.TickLabelInterpreter  = 'Latex';
gx.LongitudeAxis.TickLabelInterpreter = 'Latex';
gx.LatitudeAxis.TickLabels = strrep(gx.LatitudeAxis.TickLabels, '°', '$^{\circ}$');
gx.LongitudeAxis.TickLabels = strrep(gx.LongitudeAxis.TickLabels, '°', '$^{\circ}$');
gx.LatitudeLabel.Interpreter = 'latex';
gx.LongitudeLabel.Interpreter = 'latex';
legend1 = legend(gca,'show');
set(legend1,'Interpreter','latex','FontSize',size1,'Orientation','horizontal','Location','northoutside');
geolimits(latitudeLimits,longitudeLimits)
set(gca,'ZoomLevel',8.7609,'MapCenter',[66.2553 -146.4127])

exportgraphics(gcf,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month22\Fire\Figure_ValidationClustersInsel.pdf'),"ContentType","image")

%% extract alaska only data
% % Define the latitude and longitude boundaries of Alaska
% % Step 1: Read the shapefile
% alaska_polygon = shaperead('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\Alaska_Polygon/us-state-boundaries.shp');
% 
% % Step 2: Extract the latitude and longitude coordinates of the Alaska boundary
% % Assuming that the polygon is the first (and only) element in the shapefile structure
% alaska_lat = alaska_polygon.Y;
% alaska_lon = alaska_polygon.X;
% 
% % Now you will need the fire events data. Assuming you have two vectors, `fire_lat` and `fire_lon`,
% % which contain the latitude and longitude of the fire events
% 
% % Step 3: Use the inpolygon function to find the points inside the Alaska boundary
% inside = inpolygon(fireData_neg.longitude, fireData_neg.latitude, alaska_lon, alaska_lat);
% 
% % Filter the fire events that are inside the Alaska boundary
% alaska_fires_lat = fireData_neg.latitude(inside);
% alaska_fires_lon = fireData_neg.longitude(inside);
% 
% % If your fire events are in a table, you can use the logical index to extract the rows
% % Assuming your table is named `fireEventsTable` with 'Latitude' and 'Longitude'
% Alaska_only = fireData_neg(inside, :);

%% New temporal-spatial clustering (resulting from the jupyter python scrypt)
clear
clc
%run the script C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\Fire_clustering_preprocessingv2.m
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork',"Stage2_SpatialClusters_DBSCAN_v2.csv");
opts = detectImportOptions(path1);
fireData = readtimetable(path1,opts);%,"RowTimes",'acq_date');
fireData.Properties.VariableNames(end) = "Cluster";
fireData = removevars(fireData, "Var1");
fireData = rmmissing(fireData);


%% visualization of temporal then spatial clusters:
% clc
% % Extract unique cluster IDs
% clusters = unique(fireData.temporalCluster);
% 
% % Select 4 random clusters
% randClusters = clusters(randperm(length(clusters), 4));
% 
% % Create a tiled layout for plots
% t = tiledlayout(2, 2); % 2x2 grid
% t.Padding = 'compact';
% t.TileSpacing = 'compact';
% title(t, 'Random temporal-spatio wildfire clusters');
% 
% % Plot each randomly selected cluster
% for i = 1:length(randClusters)
%     % Extract data for the ith random cluster
%     clusterData = fireData(fireData.temporalCluster == randClusters(i), :);
% 
%     % Create a plot for each cluster in a tile
%     nexttile;
%     s = geoscatter(clusterData,"latitude", "longitude", 'filled');
%     s.ColorVariable = "Cluster";
%     colormap(gca, jet(length(unique(clusterData.Cluster))));
%     title(['Temporal cluster ' num2str(randClusters(i))]);
%     % xlabel('Latitude');
%     % ylabel('Longitude');
% end

%Figure for paper

size1 = 30;

% figure('units','normalized','outerposition',[0 0 1 1])
% Create a tiled layout for plots
% t = tiledlayout(3, 2); % 2x2 grid
% t.Padding = 'compact';
% t.TileSpacing = 'compact';
% title(t, 'Random temporal-spatio wildfire clusters');

% % Plot each randomly selected cluster
% randClusters = [1, 137, 208, 333, 459, 531];
% for i = 1:length(randClusters)
%     % Extract data for the ith random cluster
%     clusterData = fireData(fireData.temporalCluster == randClusters(i), :);
% 
%     % Create a plot for each cluster in a tile
%     % nexttile;
%     figure('units','normalized','outerposition',[0 0 1 1])
%     geoaxes("Basemap","streets-light")
%     s = geoscatter(clusterData,"latitude", "longitude", 'filled');
%     s.ColorVariable = "Cluster";
%     colormap(gca, jet(length(unique(clusterData.Cluster))));
%     set(gca,'FontSize',size1-8);
%     if min(clusterData.DetectionTime) ~= max(clusterData.DetectionTime)
%         if length(unique(clusterData.Cluster)) == 1
%             title(strcat(num2str(day(min(clusterData.DetectionTime))),{' - '},num2str(day(max(clusterData.DetectionTime))),{' '},string(month(min(clusterData.DetectionTime),'shortname')),{' '},string(year(min(clusterData.DetectionTime))),{' ('},num2str(length(unique(clusterData.Cluster))),{' cluster identified)'}),'FontSize',size1,'Interpreter','latex');
%         else
%             title(strcat(num2str(day(min(clusterData.DetectionTime))),{' - '},num2str(day(max(clusterData.DetectionTime))),{' '},string(month(min(clusterData.DetectionTime),'shortname')),{' '},string(year(min(clusterData.DetectionTime))),{' ('},num2str(length(unique(clusterData.Cluster))),{' distinct clusters identified)'}),'FontSize',size1,'Interpreter','latex');
%         end
%     else
%         if length(unique(clusterData.Cluster)) == 1
%             title(strcat(num2str(day(min(clusterData.DetectionTime))),{' '},string(month(min(clusterData.DetectionTime),'shortname')),{' '},string(year(min(clusterData.DetectionTime))),{' ('},num2str(length(unique(clusterData.Cluster))),{' cluster identified)'}),'FontSize',size1,'Interpreter','latex');        
%         else
%             title(strcat(num2str(day(min(clusterData.DetectionTime))),{' '},string(month(min(clusterData.DetectionTime),'shortname')),{' '},string(year(min(clusterData.DetectionTime))),{' ('},num2str(length(unique(clusterData.Cluster))),{' distinct clusters identified)'}),'FontSize',size1,'Interpreter','latex');        
%         end
%     exportgraphics(gcf,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month19\Fire\Figure_SpatioTemporalClusters_examples_',num2str(i),'.pdf'),"ContentType","image")
%     end
% end

% MANUALLY:
i = 3;
clusterData = fireData(fireData.temporalCluster == randClusters(i), :);

% Create a plot for the i = 1 data but zoom out to show the rest of the map
% to look interesting
f = figure('units','normalized','outerposition',[0 0 1 1]);
ax = geoaxes("Basemap","streets-light");
s = geoscatter(clusterData,"latitude", "longitude", 'filled');
s.ColorVariable = "Cluster";
colormap(gca, jet(length(unique(clusterData.Cluster))));
set(gca,'FontSize',size1-8);
if min(clusterData.DetectionTime) ~= max(clusterData.DetectionTime)
    if length(unique(clusterData.Cluster)) == 1
        title(strcat(num2str(day(min(clusterData.DetectionTime))),{' - '},num2str(day(max(clusterData.DetectionTime))),{' '},string(month(min(clusterData.DetectionTime),'shortname')),{' '},string(year(min(clusterData.DetectionTime))),{' ('},num2str(length(unique(clusterData.Cluster))),{' cluster identified)'}),'FontSize',size1,'Interpreter','latex');
    else
        title(strcat(num2str(day(min(clusterData.DetectionTime))),{' - '},num2str(day(max(clusterData.DetectionTime))),{' '},string(month(min(clusterData.DetectionTime),'shortname')),{' '},string(year(min(clusterData.DetectionTime))),{' ('},num2str(length(unique(clusterData.Cluster))),{' distinct clusters identified)'}),'FontSize',size1,'Interpreter','latex');
    end
else
    if length(unique(clusterData.Cluster)) == 1
        title(strcat(num2str(day(min(clusterData.DetectionTime))),{' '},string(month(min(clusterData.DetectionTime),'shortname')),{' '},string(year(min(clusterData.DetectionTime))),{' ('},num2str(length(unique(clusterData.Cluster))),{' cluster identified)'}),'FontSize',size1,'Interpreter','latex');        
    else
        title(strcat(num2str(day(min(clusterData.DetectionTime))),{' '},string(month(min(clusterData.DetectionTime),'shortname')),{' '},string(year(min(clusterData.DetectionTime))),{' ('},num2str(length(unique(clusterData.Cluster))),{' distinct clusters identified)'}),'FontSize',size1,'Interpreter','latex');        
    end
end

ax.LatitudeAxis.TickLabelInterpreter  = 'Latex';
ax.LongitudeAxis.TickLabelInterpreter = 'Latex';
ax.LatitudeAxis.TickLabels = strrep(ax.LatitudeAxis.TickLabels, '°', '$^{\circ}$');
ax.LongitudeAxis.TickLabels = strrep(ax.LongitudeAxis.TickLabels, '°', '$^{\circ}$');
ax.LatitudeLabel.Interpreter = 'latex';
ax.LongitudeLabel.Interpreter = 'latex';

exportgraphics(gca,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month19\Fire\Figure_SpatioTemporalClusters_examples_',num2str(i),'.pdf'),"ContentType","image")

% exportgraphics(gcf,'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month19\Fire\Figure_SpatioTemporalClusters_examples.pdf',"ContentType","image")
% close("all")

%% Compute total wildfire area and fire duration - V2 
clc
% Spatial resolution
pixelArea_km2 = 0.25; % Area in square kilometers for a 500m resolution

% Initialize variables for cluster statistics
uniqueClusters = unique(fireData.Cluster);
clusterAreas_km2 = zeros(length(uniqueClusters), 1);
clusterBrightness = zeros(length(uniqueClusters), 1);
clusterBright_t31 = zeros(length(uniqueClusters), 1);
clusterFRP = zeros(length(uniqueClusters), 1);
clusterConfidence = zeros(length(uniqueClusters), 1);
clusterDate_first = NaT(length(uniqueClusters), 1);
clusterDate_last = NaT(length(uniqueClusters), 1);
clusterDate_first.Format = "uuuu-dd-MM";
clusterDate_last.Format = "uuuu-dd-MM";
centroidLats = zeros(length(uniqueClusters), 1);
centroidLons = zeros(length(uniqueClusters), 1);

% Create a WGS84 spheroid object
wgs84 = wgs84Ellipsoid('kilometer');

% Iterate over each cluster to calculate the area
for i = 1:length(uniqueClusters)
    clusterNum = uniqueClusters(i);

    % Extract the latitude and longitude of the current cluster
    clusterIndices = fireData.Cluster == clusterNum;
    lats = fireData.latitude(clusterIndices);
    lons = fireData.longitude(clusterIndices);
    t = fireData.DetectionTime(clusterIndices);
    % Extract mean remote sensing and other features
    clusterBrightness(i) = max(fireData.brightness(clusterIndices));
    clusterBright_t31(i) = max(fireData.bright_t31(clusterIndices));
    clusterFRP(i) = max(fireData.frp(clusterIndices));
    clusterConfidence(i) = mean(fireData.confidence(clusterIndices));

    % Get first and last detection date
    clusterDate_first(i) = t(1);
    clusterDate_last(i) = t(end);

    % Calculate the mean latitude and longitude
    centroidLats(i) = mean(lats);
    centroidLons(i) = mean(lons);

    % Check the number of points in the cluster
    numPoints = length(lats);
    if numPoints == 1 | (numPoints == 2 & lats(1)-lats(2) ==0)
        % Single point represents one pixel area
        clusterAreas_km2(i) = pixelArea_km2;
    elseif numPoints == 2 & lats(1)-lats(2) ~=0
        clusterAreas_km2(i) = estimateAreaFromTwoPoints(lats, lons);
    else
        % Convert lat/lon to Cartesian coordinates
        [x, y, z] = geodetic2ecef(wgs84, lats, lons, zeros(size(lats)), 'degrees');
        try
            DT = delaunayTriangulation(x,y);
            C = convexHull(DT);
            clusterAreas_km2(i) = polyarea(DT.Points(C,1), DT.Points(C,2));
            % [k, area_km2] = convhull(x, y);
        catch
            % Handle cases where convex hull computation fails
            clusterAreas_km2(i) = NaN;
        end
    end
end
dur = between(clusterDate_first, clusterDate_last,'Days');
% Create a table with cluster statistics
clusterStats = table(uniqueClusters, clusterDate_first, clusterDate_last, double(string(regexp(cellstr(dur),'\d*','match')))+1, centroidLats, centroidLons, clusterAreas_km2, clusterBrightness, clusterBright_t31, clusterFRP, clusterConfidence,...
    'VariableNames', {'ClusterID', 'DetectionDate-first', 'DetectionDate-last','Duration','avgLatitude','avgLongitude','totalArea_km2', 'maxBrightness', 'maxBright_t31', 'maxFRP', 'avgConfidence'});
clear ans C centroidLats centroidLons clusterAreas_km2 clusterBright_t31 clusterBrightness clusterConfidence clusterDate_first clusterDate_last clusterFRP clusterIndices clusterNum DT dur dur1 i lats lons numPoints opts path1 pixelArea_km2 randClusters s size1 t uniqueClusters wgs84 x y z

writetable(clusterStats, 'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\STClusteredFires_v2.csv');

%% Plot size and duration of wildfires in geobubble
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');


path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\STClusteredFires_v2.csv');
opts = detectImportOptions(path1);
data = readtable(path1,opts);
data.Duration = categorical(data.Duration);

f = figure('outerposition',[672          -6        1114         929]);

gb = geobubble(data,'avgLatitude','avgLongitude', ...
    'SizeVariable','totalArea_km2','Basemap','topographic');
gb.SourceTable.Durationlevels = discretize(double(data.Duration),[0 5 15 30 60 73],...
                                 'categorical', {'[0-5] days', '[6-15] days', '[16-30] days', '[31-60] days', '> 60 days'});
gb.ColorVariable = 'Durationlevels';
gb.SizeLegendTitle = 'Burnt Area (km²)';
gb.ColorLegendTitle = 'Duration';

[latitudeLimits,~] = geolimits(gb);
longitudeLimits = [-180 -130];

geolimits(gb,latitudeLimits,longitudeLimits)

set(gca,'FontSize',size1-10)

% gb.LatitudeAxis.TickLabelInterpreter  = 'Latex';
% gb.LongitudeAxis.TickLabelInterpreter = 'Latex';
% gb.LatitudeAxis.TickLabels = strrep(gb.LatitudeAxis.TickLabels, '°', '$^{\circ}$');
% gb.LongitudeAxis.TickLabels = strrep(gb.LongitudeAxis.TickLabels, '°', '$^{\circ}$');
% gb.LatitudeLabel.Interpreter = 'latex';
% gb.LongitudeLabel.Interpreter = 'latex';


exportgraphics(gb,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month19\Fire\Figure_FireGeoBubble_dur_burntarea.pdf'),"ContentType","image")

%% Plot distribution of tagets (burnt area and duration)
clear
clc
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\STClusteredFires_v2.csv');
opts = detectImportOptions(path1);
data = readtable(path1,opts);


%duration
data.DurationCat = discretize(double(data.Duration),[0 5 15 30 60 max(data.Duration)],...
                                 'categorical', {'[0-5[ days', '[5-15[ days', '[15-30[ days', '[30-60[ days', '≥60 days'});
data_dur = groupsummary(data,"DurationCat");

%burntarea
data.BurntAreaCat = discretize(double(data.totalArea_km2),[0 100 500 1000 1500 max(data.totalArea_km2)],...
                                 'categorical', {'[0-100[ km^2', '[100-500[ km^2', '[500-1000[ km^2', '[1000-1500[ km^2', '≥1500 km^2'});
data_burnt = groupsummary(data,"BurntAreaCat");

close('all')
figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout(1,2,'TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);
nexttile,
d = donutchart(data_burnt,"GroupCount","BurntAreaCat",'LegendVisible','on',ExplodedWedges=1);
d.LabelStyle = "percent";
d.CenterLabel = ["Wildfire" "Burnt Area"];
d.EdgeColor = "white";
d.LineWidth = 3;
set(d,'FontName','Times','FontSize',15)

nexttile,
d = donutchart(data_dur,"GroupCount","DurationCat",'LegendVisible','on',ExplodedWedges=1);
d.LabelStyle = "percent";
d.CenterLabel = ["Wildfire" "duration"];
d.EdgeColor = "white";
d.LineWidth = 3;
set(d,'FontName','Times','FontSize',15)

% exportgraphics(gcf,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month22\fire\target_distribution.jpg'),"Resolution",300)

%% Plot correlation between wildfire behaviors (burnt area and duration)
clear
clc

size1 = 30;
cmap = get(0, 'defaultaxescolororder');

path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\STClusteredFires_v2.csv');
opts = detectImportOptions(path1);
data = readtable(path1,opts);
data.BurntAreaCat = discretize(double(data.totalArea_km2),[0 50 100 200 500 1000 inf],...
'categorical', {'$[0-50[km^2$', '$[50-100[km^2$', '$[100-200[km^2$', '$[200-500[km^2$', '$[500-1000[km^2$', '$\ge1000km^2$'});
data.DurationCat = discretize(double(data.Duration),[0 5 15 30 60 max(data.Duration)],...
                                 'categorical', {'[0-5[ days', '[5-15[ days', '[15-30[ days', '[30-60[ days', '$\ge60$ days'});

corr_all = corr(data.Duration,data.totalArea_km2);

figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout(1,2,'TileSpacing','tight','Padding','compact','units','normalized','outerposition',[0 0 1 1]);

nexttile

cats = unique(data.DurationCat);
corr_sep = zeros(size(cats));
for i = 1:length(cats)
    corr_sep(i) = corr(data.Duration(data.DurationCat==cats(i)),data.totalArea_km2(data.DurationCat==cats(i)));
end

b = bar(cats, corr_sep,'FaceColor','flat',FaceAlpha = 0.5);%,'stacked'
hold on,
yline(corr_all,'--','{All categories}','HandleVisibility','off','Interpreter','latex',...
    'LabelHorizontalAlignment','center','LabelVerticalAlignment','top',...
    'FontSize',size1-10,'LineWidth',2, 'Color',cmap(2,:));
set(gca,'Box','on','FontSize',size1-5,'TickLabelInterpreter','latex','LineWidth',0.5,'YMinorGrid','on','YMinorTick','on','TickDir','in','YGrid','on');
ylabel('Correlation coefficient','FontSize',size1,'Interpreter','latex')

set(gca,'ylim',[-0.2 1],'XTickLabelRotation',90)
%add text
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(round(b(1).YData,2));
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
'VerticalAlignment','bottom','FontSize',size1-8,Interpreter='latex',Color=[0 0 0])

nexttile

cats = unique(data.BurntAreaCat);
corr_sep = zeros(size(cats));
for i = 1:length(cats)
    corr_sep(i) = corr(data.Duration(data.BurntAreaCat==cats(i)),data.totalArea_km2(data.BurntAreaCat==cats(i)));
end

b = bar(cats, corr_sep,'FaceColor','flat',FaceAlpha = 0.5);%,'stacked'
hold on,
yline(corr_all,'--','{All categories}','HandleVisibility','off','Interpreter','latex',...
    'LabelHorizontalAlignment','center','LabelVerticalAlignment','top',...
    'FontSize',size1-10,'LineWidth',2, 'Color',cmap(2,:));
set(gca,'Box','on','FontSize',size1-5,'TickLabelInterpreter','latex','LineWidth',0.5,'YMinorGrid','on','YMinorTick','on','TickDir','in','YGrid','on');
% ylabel('Correlation coefficient','FontSize',size1,'Interpreter','latex')

set(gca,'ylim',[-0.2 1],'XTickLabelRotation',90)

%add text
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(round(b(1).YData,2));
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
'VerticalAlignment','bottom','FontSize',size1-8,Interpreter='latex',Color=[0 0 0])

exportgraphics(t,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Corr_BurntDuration.pdf'),"ContentType","image")
%% add negative class to the spatio-temporal clustered dataset v1
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork',"Fire_archive_AlaskaOnly.csv");
opts = detectImportOptions(path1);
fireData_neg = readtimetable(path1,opts);
fireData_neg = fireData_neg(fireData_neg.confidence<=5,:);
indz = randperm(height(fireData_neg),height(clusterStats));

NegClass = fireData_neg(indz,:);
NegClass = fireData_neg(indz,:);
NegClass.Properties.VariableNames(1) = "avgLatitude";
NegClass.Properties.VariableNames(2) = "avgLongitude";
NegClass = removevars(NegClass, ["scan","track","acq_time"]);
NegClass = removevars(NegClass, "instrument");
NegClass = removevars(NegClass, "version");
NegClass.Properties.VariableNames(3) = "maxBrightness";
NegClass.Properties.VariableNames(6) = "maxBright_t31";
NegClass = removevars(NegClass, "satellite");
NegClass.Properties.VariableNames(6) = "maxFRP";
NegClass.Properties.VariableNames(4) = "avgConfidence";
NegClass = removevars(NegClass, "type");
NegClass = removevars(NegClass, "daynight");
NegClass.Properties.DimensionNames{1} = 'DetectionDate-first';
NegClass.("DetectionDate-last") = NegClass.Properties.RowTimes;
NegClass.Duration = zeros(height(clusterStats),1);
NegClass.totalArea_km2 = zeros(height(clusterStats),1);
NegClass.ClusterID = (height(clusterStats)+1:2*height(clusterStats))';
NegClass = timetable2table(NegClass);
Task_FireDetection_noRS = [clusterStats table(ones(height(clusterStats),1),'VariableNames',{'isFire'});NegClass table(zeros(height(clusterStats),1),'VariableNames',{'isFire'})];
writetable(Task_FireDetection_noRS, 'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\Task_FireDetection_noRS.csv');


%% add negative class to the spatio-temporal clustered dataset v2
clc
clear
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\STClusteredFires_v2.csv');
opts = detectImportOptions(path1);
% opts.VariableNames(2:3) = {"DetectionDatefirst", "DetectionDatelast"};
opts.VariableTypes(2:3) = {'datetime','datetime'};
opts = setvaropts(opts, "DetectionDate_first", "InputFormat", "yyyy-dd-MM");
opts = setvaropts(opts, "DetectionDate_last", "InputFormat", "yyyy-dd-MM");

fireData_neg = readtable(path1,opts);

fireData_neg.('Duration') = zeros(height(fireData_neg),1);
fireData_neg.('totalArea_km2') = zeros(height(fireData_neg),1);
fireData_neg.('maxBrightness') = zeros(height(fireData_neg),1);
fireData_neg.('maxBright_t31') = zeros(height(fireData_neg),1);
fireData_neg.('maxFRP') = zeros(height(fireData_neg),1);
fireData_neg.('avgConfidence') = zeros(height(fireData_neg),1);
fireData_neg.('ClusterID') = (height(fireData_neg)+1:2*height(fireData_neg))';

%Loop through each row in the table
for i = 1:height(fireData_neg)
    % Generate a random integer between -120 and -90
    shiftDays = randi([-120, -90]);
    
    % Shift the date by the randomly generated number of days
    fireData_neg.DetectionDate_first(i) = fireData_neg.DetectionDate_first(i) + days(shiftDays);
    fireData_neg.DetectionDate_last(i) = fireData_neg.DetectionDate_first(i) + days(shiftDays);
end

writetable(fireData_neg, 'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\Task_FireDetection_noRS_negClass_datesShifted_v2.csv');

%% add negative class to the spatio-temporal clustered dataset v4
clc
clear
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\STClusteredFires_v2.csv');
opts = detectImportOptions(path1);
% opts.VariableNames(2:3) = {"DetectionDatefirst", "DetectionDatelast"};
opts.VariableTypes(2:3) = {'datetime','datetime'};
opts = setvaropts(opts, "DetectionDate_first", "InputFormat", "yyyy-dd-MM");
opts = setvaropts(opts, "DetectionDate_last", "InputFormat", "yyyy-dd-MM");

fireData_neg = readtable(path1,opts);

fireData_neg.('Duration') = zeros(height(fireData_neg),1);
fireData_neg.('maxBrightness') = zeros(height(fireData_neg),1);
fireData_neg.('maxBright_t31') = zeros(height(fireData_neg),1);
fireData_neg.('maxFRP') = zeros(height(fireData_neg),1);
fireData_neg.('avgConfidence') = zeros(height(fireData_neg),1);
fireData_neg.('ClusterID') = (height(fireData_neg)+1:2*height(fireData_neg))';

%Loop through each row in the table
for i = 1:height(fireData_neg)
    shiftLoc_lat = randchangeLat(sqrt(fireData_neg.totalArea_km2(i))) + randchangeLat(15);
    % shiftLoc_lon = randchangeLon(15, fireData_neg.avgLatitude(i));
    
    % Shift the date by the randomly generated number of days
    fireData_neg.avgLatitude(i) = fireData_neg.avgLatitude(i) + shiftLoc_lat;
    % fireData_neg.avgLongitude(i) = fireData_neg.avgLongitude(i) + shiftLoc_lon;
end
% fireData_neg.('totalArea_km2') = zeros(height(fireData_neg),1);

writetable(fireData_neg, 'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\Task_FireDetection_noRS_negClass_V4.csv');

%% remote sensing preprocessing - posClass
clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  regularize RS features  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
features = {'LC_Type1','LC_Type2','LC_Type3','LC_Type4','LC_Type5','Elevation', 'NDVI', 'EVI', 'sur_refl_b01', 'sur_refl_b02', 'sur_refl_b03', 'sur_refl_b07', 'Clear_day_cov','Clear_night_cov','Albedo_BSA_nir','Albedo_BSA_shortwave','Albedo_WSA_nir','Albedo_WSA_shortwave','prcp', 'srad', 'dayl', 'tmax', 'tmin','LST_Day_1km','LST_Night_1km'};

fprintf('Regularizing: ')
for k=1:length(features)
    fprintf('%s, ',features{k})
    if isequal(features{k},'LC_Type1') || isequal(features{k},'LC_Type2') || isequal(features{k},'LC_Type3') || isequal(features{k},'LC_Type4') || isequal(features{k},'LC_Type5') || isequal(features{k},'Elevation')
        path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',strcat(features{k},'_values.csv'));
        opts = detectImportOptions(path1);
        tmpp = readtable(path1,opts);
        eval(sprintf('%s_values = (repmat(tmpp.%s,[1,91]));',features{k},features{k}))
        eval(sprintf('%s_values = fillmissing(%s_values,''nearest'');',features{k},features{k}));
    else
        %load values of an RS feature
        path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',strcat(features{k},'_values.csv'));
        opts = detectImportOptions(path1);
        values = readtimetable(path1,opts);
        values(1,:) = []; %make sure this is needed

        path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',strcat(features{k},'_dates.csv'));
        opts = detectImportOptions(path1);
        opts.VariableTypes = repmat({'datetime'},[1, length(opts.VariableNames)]);%convert all content to datetime
        dates = readtimetable(path1,opts);
        dates(1,:) = []; %make sure this is needed

        %make a regular matrix of time sequences
        times = NaT(height(dates),91);


        %synchronize
        eval(sprintf('%s_values = zeros(height(dates),91);',features{k}))
        eval(sprintf('%s_dates = NaT(height(dates),91);',features{k}))
        for i = 1:height(dates)
            dt = dates.Properties.RowTimes(i);
            eval(sprintf('tmp = rmmissing(timetable(dates{i,:}'',values{i,:}'',''VariableNames'',{''%s''}));',features{k}))
            tmp_regularized = retime(tmp,dt - days(90):days(1):dt,"nearest");
            eval(sprintf('%s_values(i,:) = tmp_regularized.%s;',features{k},features{k}))
        end
    end
end
clear values tmp tmpp tmp_regularized opts i dt dates times path1 k

%%% Export all features for model learning in python
for i = 1:length(features)
    eval(sprintf('writematrix(%s_values,''%s'');',features{i},strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Fires_task1_',features{i},'.csv')))
end
%% remote sensing preprocessing - negClass
clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  regularize RS features  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
features = {'LC_Type1','LC_Type2','LC_Type3','LC_Type4','LC_Type5','Elevation', 'NDVI', 'EVI', 'sur_refl_b01', 'sur_refl_b02', 'sur_refl_b03', 'sur_refl_b07', 'Clear_day_cov','Clear_night_cov','Albedo_BSA_nir','Albedo_BSA_shortwave','Albedo_WSA_nir','Albedo_WSA_shortwave','prcp', 'srad', 'dayl', 'tmax', 'tmin','LST_Day_1km','LST_Night_1km'};

fprintf('Regularizing: ')
for k=1:length(features)
    fprintf('%s, ',features{k})
    if isequal(features{k},'LC_Type1') || isequal(features{k},'LC_Type2') || isequal(features{k},'LC_Type3') || isequal(features{k},'LC_Type4') || isequal(features{k},'LC_Type5') || isequal(features{k},'Elevation')
        path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',strcat(features{k},'_negClass_values.csv'));
        opts = detectImportOptions(path1);
        tmpp = readtable(path1,opts);
        eval(sprintf('%s_values = (repmat(tmpp.%s,[1,91]));',features{k},features{k}))
        eval(sprintf('%s_values = fillmissing(%s_values,''nearest'');',features{k},features{k}));
    else
        %load values of an RS feature
        path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',strcat(features{k},'_negClass_values.csv'));
        opts = detectImportOptions(path1);
        values = readtimetable(path1,opts);
        values(1,:) = []; %make sure this is needed

        path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',strcat(features{k},'_negClass_dates.csv'));
        opts = detectImportOptions(path1);
        opts.VariableTypes = repmat({'datetime'},[1, length(opts.VariableNames)]);%convert all content to datetime
        dates = readtimetable(path1,opts);
        dates(1,:) = []; %make sure this is needed

        %make a regular matrix of time sequences
        times = NaT(height(dates),91);


        %synchronize
        eval(sprintf('%s_values = zeros(height(dates),91);',features{k}))
        eval(sprintf('%s_dates = NaT(height(dates),91);',features{k}))
        for i = 1:height(dates)
            dt = dates.Properties.RowTimes(i);
            eval(sprintf('tmp = rmmissing(timetable(dates{i,:}'',values{i,:}'',''VariableNames'',{''%s''}));',features{k}))
            tmp_regularized = retime(tmp,dt - days(90):days(1):dt,"nearest");
            eval(sprintf('%s_values(i,:) = tmp_regularized.%s;',features{k},features{k}))
        end
    end
end
clear values tmp tmpp tmp_regularized opts i dt dates times path1 k

%%% Export all features for model learning in python
for i = 1:length(features)
    eval(sprintf('writematrix(%s_values,''%s'');',features{i},strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Fires_task1_negClass_',features{i},'.csv')))
end
%% remote sensing preprocessing - Both Classes (Elevation, Slope, Aspect)
clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  regularize RS features  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
features = {'Elevation', 'Slope', 'Aspect'};

fprintf('Regularizing: ')
for k=1:length(features)
    fprintf('%s, ',features{k})
    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',strcat(features{k},'_values.csv'));
    opts = detectImportOptions(path1);
    tmpp = readtable(path1,opts);
    eval(sprintf('%s_values = (repmat(tmpp.%s,[1,91]));',features{k},features{k}))
    eval(sprintf('%s_values = fillmissing(%s_values,''nearest'');',features{k},features{k}));
end
clear values tmp tmpp tmp_regularized opts i dt dates times path1 k

%%% Export all features for model learning in python
for i = 1:length(features)
    eval(sprintf('writematrix(%s_values,''%s'');',features{i},strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Fires_task1_',features{i},'.csv')))
end
%% Location features preprocessing - Both Classes (lat & lon)
clear
clc
%%% load posclass data
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork',"Task_FireDetection_noRS.csv");
opts = detectImportOptions(path1);
opts.VariableTypes(2:3) = repmat({'datetime'},[1, length(opts.VariableNames(2:3))]);%convert all content to datetime
opts = setvaropts(opts, "DetectionDate_first", "InputFormat", "yyyy-dd-MM");
opts = setvaropts(opts, "DetectionDate_last", "InputFormat", "yyyy-dd-MM");
Fires_posClass = readtable(path1,opts);
Fires_posClass(Fires_posClass.isFire==0,:) = [];

%load negclass data
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork',"Task_FireDetection_noRS_negClass_datesShifted_v2.csv");
opts = detectImportOptions(path1);
Fires_negClass = readtable(path1,opts);
Fires_negClass.isFire = zeros(height(Fires_negClass),1);

%Fuze targets
Fires = [Fires_posClass;Fires_negClass];
clear Fires_posClass Fires_negClass opts path1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  regularize Lat Lon coord. features  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
features = {'avgLatitude', 'avgLongitude'};

fprintf('Regularizing: ')
for k=1:length(features)
    fprintf('%s, ',features{k})
    eval(sprintf('%s_values = (repmat(Fires.%s,[1,91]));',features{k},features{k}))
    eval(sprintf('%s_values = fillmissing(%s_values,''nearest'');',features{k},features{k}));
end
clear values tmp tmpp tmp_regularized opts i dt dates times path1 k Fires

%%% Export all features for model learning in python
for i = 1:length(features)
    eval(sprintf('writematrix(%s_values,''%s'');',features{i},strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Fires_task1_',features{i},'.csv')))
end
%% create time-related features
clear
clc
%read TS 
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2','NDVI_dates.csv');
opts = detectImportOptions(path1);
opts.VariableTypes = repmat({'datetime'},[1, length(opts.VariableNames)]);%convert all content to datetime
dates = readtimetable(path1,opts);
dates(1,:) = []; %make sure this is needed

%make a regular matrix of time sequences
times = NaT(height(dates),91);

%synchronize
TS = NaT(height(dates),91);
for i = 1:height(dates)
    dt = dates.Properties.RowTimes(i);
    TS(i,:) = dt - days(90):days(1):dt;
end
clear times opts path1 dt dates

features = {'Y','M', 'D', 'W', 'H', 'S'};
for i = 1:length(features)
    eval(sprintf('%s  = zeros(size(TS));',features{i}));
end

minTS = min(TS,[],'all');
maxTS = max(TS,[],'all');
hldy_list = holidays(minTS,maxTS);
[Y,M,D] = get_YMD(TS);
W = isweekend(TS); %weekdays
[H,~] = ismember(TS,hldy_list);    %holidays    
S = get_astronomical_season(TS); %Season    

features = {'TS','Y','M', 'D', 'W', 'H', 'S'};
for i = 1:length(features)
    eval(sprintf('writematrix(%s,''%s'');',features{i},strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Fires_task1_',features{i},'.csv')))
end
%% create time-related features - NegClass
clear
clc
%read TS 
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2','NDVI_NegClass_dates.csv');
opts = detectImportOptions(path1);
opts.VariableTypes = repmat({'datetime'},[1, length(opts.VariableNames)]);%convert all content to datetime
dates = readtimetable(path1,opts);
dates(1,:) = []; %make sure this is needed

%make a regular matrix of time sequences
times = NaT(height(dates),91);

%synchronize
TS = NaT(height(dates),91);
for i = 1:height(dates)
    dt = dates.Properties.RowTimes(i);
    TS(i,:) = dt - days(90):days(1):dt;
end
clear times opts path1 dt dates

features = {'Y','M', 'D', 'W', 'H', 'S'};
for i = 1:length(features)
    eval(sprintf('%s  = zeros(size(TS));',features{i}));
end

minTS = min(TS,[],'all');
maxTS = max(TS,[],'all');
hldy_list = holidays(minTS,maxTS);
[Y,M,D] = get_YMD(TS);
W = isweekend(TS); %weekdays
[H,~] = ismember(TS,hldy_list);    %holidays    
S = get_astronomical_season(TS); %Season    

features = {'TS','Y','M', 'D', 'W', 'H', 'S'};
for i = 1:length(features)
    eval(sprintf('writematrix(%s,''%s'');',features{i},strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Fires_task1_NegClass_',features{i},'.csv')))
end
%% Reanalysis features preprocessing - both classes
clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  regularize RS features  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
features = {'v_component_of_wind_10m','u_component_of_wind_10m','surface_pressure','runoff_sum','total_evaporation_sum','snowfall_sum','snowmelt_sum','snow_depth','snow_density','snow_cover','dewpoint_temperature_2m','temperature_2m'};

fprintf('Regularizing: ')
for k=1:length(features)
    fprintf('%s, ',features{k})
    %load values of an RS feature
    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',strcat(features{k},'_values.csv'));
    opts = detectImportOptions(path1);
    values = readtimetable(path1,opts);
    values(1,:) = []; %make sure this is needed

    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',strcat(features{k},'_dates.csv'));
    opts = detectImportOptions(path1);
    opts.VariableTypes = repmat({'datetime'},[1, length(opts.VariableNames)]);%convert all content to datetime
    dates = readtimetable(path1,opts);
    dates(1,:) = []; %make sure this is needed

    %make a regular matrix of time sequences
    times = NaT(height(dates),91);


    %synchronize
    eval(sprintf('%s_values = zeros(height(dates),91);',features{k}))
    eval(sprintf('%s_dates = NaT(height(dates),91);',features{k}))
    for i = 1:height(dates)
        dt = dates.Properties.RowTimes(i);
        eval(sprintf('tmp = rmmissing(timetable(dates{i,:}'',values{i,:}'',''VariableNames'',{''%s''}));',features{k}))
        tmp_regularized = retime(tmp,dt - days(90):days(1):dt,"nearest");
        eval(sprintf('%s_values(i,:) = tmp_regularized.%s;',features{k},features{k}))
    end
end
clear values tmp tmpp tmp_regularized opts i dt dates times path1 k

%%% Export all features for model learning in python
for i = 1:length(features)
    eval(sprintf('writematrix(%s_values,''%s'');',features{i},strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Fires_task1_',features{i},'.csv')))
end

%% input feature preprocessing - negClassV3
clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  regularize RS features  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
features = {'LC_Type1','LC_Type2','LC_Type3','LC_Type4','LC_Type5', 'NDVI', 'EVI', 'sur_refl_b01', 'sur_refl_b02', 'sur_refl_b03', 'sur_refl_b07', 'Clear_day_cov','Clear_night_cov','Albedo_BSA_nir','Albedo_BSA_shortwave','Albedo_WSA_nir','Albedo_WSA_shortwave','prcp', 'srad', 'dayl', 'tmax', 'tmin','LST_Day_1km','LST_Night_1km'};

fprintf('Regularizing: ')
for k=1:length(features)
    fprintf('%s, ',features{k})
    if isequal(features{k},'LC_Type1') || isequal(features{k},'LC_Type2') || isequal(features{k},'LC_Type3') || isequal(features{k},'LC_Type4') || isequal(features{k},'LC_Type5') || isequal(features{k},'Elevation')
        path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',strcat(features{k},'_negClassV3_values.csv'));
        opts = detectImportOptions(path1);
        tmpp = readtable(path1,opts);
        eval(sprintf('%s_values = (repmat(tmpp.%s,[1,91]));',features{k},features{k}))
        eval(sprintf('%s_values = fillmissing(%s_values,''nearest'');',features{k},features{k}));
    else
        %load values of an RS feature
        path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',strcat(features{k},'_negClassV3_values.csv'));
        opts = detectImportOptions(path1);
        values = readtimetable(path1,opts);
        values(1,:) = []; %make sure this is needed

        path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',strcat(features{k},'_negClassV3_dates.csv'));
        opts = detectImportOptions(path1);
        opts.VariableTypes = repmat({'datetime'},[1, length(opts.VariableNames)]);%convert all content to datetime
        dates = readtimetable(path1,opts);
        dates(1,:) = []; %make sure this is needed

        %make a regular matrix of time sequences
        times = NaT(height(dates),91);


        %synchronize
        eval(sprintf('%s_values = zeros(height(dates),91);',features{k}))
        eval(sprintf('%s_dates = NaT(height(dates),91);',features{k}))
        for i = 1:height(dates)
            dt = dates.Properties.RowTimes(i);
            eval(sprintf('tmp = rmmissing(timetable(dates{i,:}'',values{i,:}'',''VariableNames'',{''%s''}));',features{k}))
            tmp_regularized = retime(tmp,dt - days(90):days(1):dt,"nearest");
            eval(sprintf('%s_values(i,:) = tmp_regularized.%s;',features{k},features{k}))
        end
    end
end
clear values tmp tmpp tmp_regularized opts i dt dates times path1 k

%%% Export all features for model learning in python
for i = 1:length(features)
    eval(sprintf('writematrix(%s_values,''%s'');',features{i},strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Fires_task1_negClassV3_',features{i},'.csv')))
end

% remote sensing preprocessing - negClassV3 (Elevation, Slope, Aspect)
clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  regularize RS features  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
features = {'Elevation', 'Slope', 'Aspect'};

fprintf('Regularizing: ')
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2','Elevation_Slope_Aspect_negClassV3_values.csv');
opts = detectImportOptions(path1);
tmpp = readtable(path1,opts);
for k=1:length(features)
    fprintf('%s, ',features{k})
    eval(sprintf('%s_values = (repmat(tmpp.%s,[1,91]));',features{k},features{k}))
    eval(sprintf('%s_values = fillmissing(%s_values,''nearest'');',features{k},features{k}));
end
clear values tmp tmpp tmp_regularized opts i dt dates times path1 k

%%% Export all features for model learning in python
for i = 1:length(features)
    eval(sprintf('writematrix(%s_values,''%s'');',features{i},strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Fires_task1_negClassV3_',features{i},'.csv')))
end

% Location features preprocessing - negClassV3 (lat & lon)
clear
clc
% %%% load posclass data
% path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork',"Task_FireDetection_noRS.csv");
% opts = detectImportOptions(path1);
% opts.VariableTypes(2:3) = repmat({'datetime'},[1, length(opts.VariableNames(2:3))]);%convert all content to datetime
% opts = setvaropts(opts, "DetectionDate_first", "InputFormat", "yyyy-dd-MM");
% opts = setvaropts(opts, "DetectionDate_last", "InputFormat", "yyyy-dd-MM");
% Fires_posClass = readtable(path1,opts);
% Fires_posClass(Fires_posClass.isFire==0,:) = [];

%load negclass data
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork',"Task_FireDetection_noRS_negClass_V3.csv");
opts = detectImportOptions(path1);
Fires_negClass = readtable(path1,opts);
Fires_negClass.isFire = zeros(height(Fires_negClass),1);

%Fuze targets
Fires = [Fires_negClass];%Fires_posClass;
clear Fires_posClass Fires_negClass opts path1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  regularize Lat Lon coord. features  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
features = {'avgLatitude', 'avgLongitude'};

fprintf('Regularizing: ')
for k=1:length(features)
    fprintf('%s, ',features{k})
    eval(sprintf('%s_values = (repmat(Fires.%s,[1,91]));',features{k},features{k}))
    eval(sprintf('%s_values = fillmissing(%s_values,''nearest'');',features{k},features{k}));
end
clear values tmp tmpp tmp_regularized opts i dt dates times path1 k Fires

%%% Export all features for model learning in python
for i = 1:length(features)
    eval(sprintf('writematrix(%s_values,''%s'');',features{i},strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Fires_task1_negClassV3_',features{i},'.csv')))
end

% create time-related features - NegClassV3
clear
clc
%read TS 
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2','NDVI_negClassV3_dates.csv');
opts = detectImportOptions(path1);
opts.VariableTypes = repmat({'datetime'},[1, length(opts.VariableNames)]);%convert all content to datetime
dates = readtimetable(path1,opts);
dates(1,:) = []; %make sure this is needed

%make a regular matrix of time sequences
times = NaT(height(dates),91);

%synchronize
TS = NaT(height(dates),91);
for i = 1:height(dates)
    dt = dates.Properties.RowTimes(i);
    TS(i,:) = dt - days(90):days(1):dt;
end
clear times opts path1 dt dates

features = {'Y','M', 'D', 'W', 'H', 'S'};
for i = 1:length(features)
    eval(sprintf('%s  = zeros(size(TS));',features{i}));
end

minTS = min(TS,[],'all');
maxTS = max(TS,[],'all');
hldy_list = holidays(minTS,maxTS);
[Y,M,D] = get_YMD(TS);
W = isweekend(TS); %weekdays
[H,~] = ismember(TS,hldy_list);    %holidays    
S = get_astronomical_season(TS); %Season    

features = {'TS','Y','M', 'D', 'W', 'H', 'S'};
for i = 1:length(features)
    eval(sprintf('writematrix(%s,''%s'');',features{i},strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Fires_task1_negClassV3_',features{i},'.csv')))
end

%Reanalysis features preprocessing - negClass V3
clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  regularize RS features  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
features = {'v_component_of_wind_10m','u_component_of_wind_10m','surface_pressure','runoff_sum','total_evaporation_sum','snowfall_sum','snowmelt_sum','snow_depth','snow_density','snow_cover','dewpoint_temperature_2m','temperature_2m'};

fprintf('Regularizing: ')
for k=1:length(features)
    fprintf('%s, ',features{k})
    %load values of an RS feature
    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',strcat(features{k},'_negClassV3_values.csv'));
    opts = detectImportOptions(path1);
    values = readtimetable(path1,opts);
    values(1,:) = []; %make sure this is needed

    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',strcat(features{k},'_negClassV3_dates.csv'));
    opts = detectImportOptions(path1);
    opts.VariableTypes = repmat({'datetime'},[1, length(opts.VariableNames)]);%convert all content to datetime
    dates = readtimetable(path1,opts);
    dates(1,:) = []; %make sure this is needed

    %make a regular matrix of time sequences
    times = NaT(height(dates),91);


    %synchronize
    eval(sprintf('%s_values = zeros(height(dates),91);',features{k}))
    eval(sprintf('%s_dates = NaT(height(dates),91);',features{k}))
    for i = 1:height(dates)
        dt = dates.Properties.RowTimes(i);
        eval(sprintf('tmp = rmmissing(timetable(dates{i,:}'',values{i,:}'',''VariableNames'',{''%s''}));',features{k}))
        tmp_regularized = retime(tmp,dt - days(90):days(1):dt,"nearest");
        eval(sprintf('%s_values(i,:) = tmp_regularized.%s;',features{k},features{k}))
    end
end
clear values tmp tmpp tmp_regularized opts i dt dates times path1 k

%%% Export all features for model learning in python
for i = 1:length(features)
    eval(sprintf('writematrix(%s_values,''%s'');',features{i},strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Fires_task1_negClassV3_',features{i},'.csv')))
end
%% input feature preprocessing - negClassV4
clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  regularize RS features  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
features = {'LC_Type1','LC_Type2','LC_Type3','LC_Type4','LC_Type5', 'NDVI', 'EVI', 'sur_refl_b01', 'sur_refl_b02', 'sur_refl_b03', 'sur_refl_b07', 'Clear_day_cov','Clear_night_cov','Albedo_BSA_nir','Albedo_BSA_shortwave','Albedo_WSA_nir','Albedo_WSA_shortwave','prcp', 'srad', 'dayl', 'tmax', 'tmin','LST_Day_1km','LST_Night_1km'};

fprintf('Regularizing: ')
for k=1:length(features)
    fprintf('%s, ',features{k})
    if isequal(features{k},'LC_Type1') || isequal(features{k},'LC_Type2') || isequal(features{k},'LC_Type3') || isequal(features{k},'LC_Type4') || isequal(features{k},'LC_Type5') || isequal(features{k},'Elevation')
        path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',strcat(features{k},'_negClassV4_values.csv'));
        opts = detectImportOptions(path1);
        tmpp = readtable(path1,opts);
        eval(sprintf('%s_values = (repmat(tmpp.%s,[1,91]));',features{k},features{k}))
        eval(sprintf('%s_values = fillmissing(%s_values,''nearest'');',features{k},features{k}));
    else
        %load values of an RS feature
        path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',strcat(features{k},'_negClassV4_values.csv'));
        opts = detectImportOptions(path1);
        values = readtimetable(path1,opts);
        values(1,:) = []; %make sure this is needed

        path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',strcat(features{k},'_negClassV4_dates.csv'));
        opts = detectImportOptions(path1);
        opts.VariableTypes = repmat({'datetime'},[1, length(opts.VariableNames)]);%convert all content to datetime
        dates = readtimetable(path1,opts);
        dates(1,:) = []; %make sure this is needed

        %make a regular matrix of time sequences
        times = NaT(height(dates),91);


        %synchronize
        eval(sprintf('%s_values = zeros(height(dates),91);',features{k}))
        eval(sprintf('%s_dates = NaT(height(dates),91);',features{k}))
        for i = 1:height(dates)
            dt = dates.Properties.RowTimes(i);
            eval(sprintf('tmp = rmmissing(timetable(dates{i,:}'',values{i,:}'',''VariableNames'',{''%s''}));',features{k}))
            tmp_regularized = retime(tmp,dt - days(90):days(1):dt,"nearest");
            eval(sprintf('%s_values(i,:) = tmp_regularized.%s;',features{k},features{k}))
        end
    end
end
clear values tmp tmpp tmp_regularized opts i dt dates times path1 k

%%% Export all features for model learning in python
for i = 1:length(features)
    eval(sprintf('writematrix(%s_values,''%s'');',features{i},strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Fires_task1_negClassV4_',features{i},'.csv')))
end

% remote sensing preprocessing - negClassV4 (Elevation, Slope, Aspect)
clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  regularize RS features  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
features = {'Elevation', 'Slope', 'Aspect'};

fprintf('Regularizing: ')
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2','Elevation_Slope_Aspect_negClassV4_values.csv');
opts = detectImportOptions(path1);
tmpp = readtable(path1,opts);
for k=1:length(features)
    fprintf('%s, ',features{k})
    eval(sprintf('%s_values = (repmat(tmpp.%s,[1,91]));',features{k},features{k}))
    eval(sprintf('%s_values = fillmissing(%s_values,''nearest'');',features{k},features{k}));
end
clear values tmp tmpp tmp_regularized opts i dt dates times path1 k

%%% Export all features for model learning in python
for i = 1:length(features)
    eval(sprintf('writematrix(%s_values,''%s'');',features{i},strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Fires_task1_negClassV4_',features{i},'.csv')))
end

% Location features preprocessing - negClassV4 (lat & lon)
clear
clc
% %%% load posclass data
% path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork',"Task_FireDetection_noRS.csv");
% opts = detectImportOptions(path1);
% opts.VariableTypes(2:3) = repmat({'datetime'},[1, length(opts.VariableNames(2:3))]);%convert all content to datetime
% opts = setvaropts(opts, "DetectionDate_first", "InputFormat", "yyyy-dd-MM");
% opts = setvaropts(opts, "DetectionDate_last", "InputFormat", "yyyy-dd-MM");
% Fires_posClass = readtable(path1,opts);
% Fires_posClass(Fires_posClass.isFire==0,:) = [];

%load negclass data
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork',"Task_FireDetection_noRS_negClass_V4.csv");
opts = detectImportOptions(path1);
Fires_negClass = readtable(path1,opts);
Fires_negClass.isFire = zeros(height(Fires_negClass),1);

%Fuze targets
Fires = [Fires_negClass];%Fires_posClass;
clear Fires_posClass Fires_negClass opts path1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  regularize Lat Lon coord. features  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
features = {'avgLatitude', 'avgLongitude'};

fprintf('Regularizing: ')
for k=1:length(features)
    fprintf('%s, ',features{k})
    eval(sprintf('%s_values = (repmat(Fires.%s,[1,91]));',features{k},features{k}))
    eval(sprintf('%s_values = fillmissing(%s_values,''nearest'');',features{k},features{k}));
end
clear values tmp tmpp tmp_regularized opts i dt dates times path1 k Fires

%%% Export all features for model learning in python
for i = 1:length(features)
    eval(sprintf('writematrix(%s_values,''%s'');',features{i},strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Fires_task1_negClassV4_',features{i},'.csv')))
end

% create time-related features - NegClassV4
clear
clc
%read TS 
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2','NDVI_negClassV4_dates.csv');
opts = detectImportOptions(path1);
opts.VariableTypes = repmat({'datetime'},[1, length(opts.VariableNames)]);%convert all content to datetime
dates = readtimetable(path1,opts);
dates(1,:) = []; %make sure this is needed

%make a regular matrix of time sequences
times = NaT(height(dates),91);

%synchronize
TS = NaT(height(dates),91);
for i = 1:height(dates)
    dt = dates.Properties.RowTimes(i);
    TS(i,:) = dt - days(90):days(1):dt;
end
clear times opts path1 dt dates

features = {'Y','M', 'D', 'W', 'H', 'S'};
for i = 1:length(features)
    eval(sprintf('%s  = zeros(size(TS));',features{i}));
end

minTS = min(TS,[],'all');
maxTS = max(TS,[],'all');
hldy_list = holidays(minTS,maxTS);
[Y,M,D] = get_YMD(TS);
W = isweekend(TS); %weekdays
[H,~] = ismember(TS,hldy_list);    %holidays    
S = get_astronomical_season(TS); %Season    

features = {'TS','Y','M', 'D', 'W', 'H', 'S'};
for i = 1:length(features)
    eval(sprintf('writematrix(%s,''%s'');',features{i},strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Fires_task1_negClassV4_',features{i},'.csv')))
end

%Reanalysis features preprocessing - negClass V4
clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  regularize RS features  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
features = {'v_component_of_wind_10m','u_component_of_wind_10m','surface_pressure','runoff_sum','total_evaporation_sum','snowfall_sum','snowmelt_sum','snow_depth','snow_density','snow_cover','dewpoint_temperature_2m','temperature_2m'};

fprintf('Regularizing: ')
for k=1:length(features)
    fprintf('%s, ',features{k})
    %load values of an RS feature
    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',strcat(features{k},'_negClassV4_values.csv'));
    opts = detectImportOptions(path1);
    values = readtimetable(path1,opts);
    values(1,:) = []; %make sure this is needed

    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',strcat(features{k},'_negClassV4_dates.csv'));
    opts = detectImportOptions(path1);
    opts.VariableTypes = repmat({'datetime'},[1, length(opts.VariableNames)]);%convert all content to datetime
    dates = readtimetable(path1,opts);
    dates(1,:) = []; %make sure this is needed

    %make a regular matrix of time sequences
    times = NaT(height(dates),91);


    %synchronize
    eval(sprintf('%s_values = zeros(height(dates),91);',features{k}))
    eval(sprintf('%s_dates = NaT(height(dates),91);',features{k}))
    for i = 1:height(dates)
        dt = dates.Properties.RowTimes(i);
        eval(sprintf('tmp = rmmissing(timetable(dates{i,:}'',values{i,:}'',''VariableNames'',{''%s''}));',features{k}))
        tmp_regularized = retime(tmp,dt - days(90):days(1):dt,"nearest");
        eval(sprintf('%s_values(i,:) = tmp_regularized.%s;',features{k},features{k}))
    end
end
clear values tmp tmpp tmp_regularized opts i dt dates times path1 k

%%% Export all features for model learning in python
for i = 1:length(features)
    eval(sprintf('writematrix(%s_values,''%s'');',features{i},strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Fires_task1_negClassV4_',features{i},'.csv')))
end

%% 1- Read all inputs and targets

clear
clc

fprintf('Loading and processing the historical fires dataset...\n')

%%% load posclass data
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork',"Task_FireDetection_noRS.csv");
opts = detectImportOptions(path1);
opts.VariableTypes(2:3) = repmat({'datetime'},[1, length(opts.VariableNames(2:3))]);%convert all content to datetime
opts = setvaropts(opts, "DetectionDate_first", "InputFormat", "yyyy-dd-MM");
opts = setvaropts(opts, "DetectionDate_last", "InputFormat", "yyyy-dd-MM");
Fires_posClass = readtable(path1,opts);
Fires_posClass(Fires_posClass.isFire==0,:) = [];

%load negclass data
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork',"Task_FireDetection_noRS_negClass_datesShifted_v2.csv");
opts = detectImportOptions(path1);
Fires_negClass = readtable(path1,opts);
Fires_negClass.isFire = zeros(height(Fires_negClass),1);

%Fuze targets
Fires = [Fires_posClass;Fires_negClass];
clear Fires_posClass Fires_negClass

%%% load input RS features
fprintf('Loading and processing the remote sensing features...\n')

features = {'LC_Type1','LC_Type2','LC_Type3','LC_Type4','LC_Type5', 'NDVI', 'EVI', 'sur_refl_b01', 'sur_refl_b02', 'sur_refl_b03', 'sur_refl_b07', 'Clear_day_cov','Clear_night_cov','Albedo_BSA_nir','Albedo_BSA_shortwave','Albedo_WSA_nir','Albedo_WSA_shortwave','prcp', 'srad', 'dayl', 'tmax', 'tmin','LST_Day_1km','LST_Night_1km','Y','M', 'D', 'W', 'H', 'S'};
for i = 1:length(features)
    %posClass
    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\',strcat('Fires_task1_',features{i}));
    opts = detectImportOptions(path1);
    eval(sprintf('%s  = readtable(path1,opts);',features{i}));
    %negClass
    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\',strcat('Fires_task1_negClass_',features{i}));
    opts = detectImportOptions(path1);
    eval(sprintf('%s_negClass  = readtable(path1,opts);',features{i}));
    %fuze data
    eval(sprintf('%s  = [%s(1:height(%s)/2,:);%s_negClass];',features{i},features{i},features{i},features{i}));
    eval(sprintf('clear %s_negClass',features{i}))
end
clear opts path1 i 

%%% load input Reanalysis features and topographic remote sensing features and coords.
fprintf('Loading and processing the Reanalysis features...\n')

features = {'avgLongitude', 'avgLatitude', 'Elevation', 'Slope', 'Aspect','v_component_of_wind_10m','u_component_of_wind_10m','surface_pressure','runoff_sum','total_evaporation_sum','snowfall_sum','snowmelt_sum','snow_depth','snow_density','snow_cover','dewpoint_temperature_2m','temperature_2m'};
for i = 1:length(features)
    %posClass
    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\',strcat('Fires_task1_',features{i}));
    opts = detectImportOptions(path1);
    eval(sprintf('%s  = readtable(path1,opts);',features{i}));
end
clear opts path1 i 

fprintf('Fuzing all data...\n')

% Create boolean masks for rows with all NaNs in each feature table (all RS and Reanalysis features)
features = {'LC_Type1','LC_Type2','LC_Type3','LC_Type4','LC_Type5', 'avgLongitude', 'avgLatitude', 'Elevation', 'Slope', 'Aspect', 'NDVI', 'EVI', 'sur_refl_b01', 'sur_refl_b02', 'sur_refl_b03', 'sur_refl_b07', 'Clear_day_cov','Clear_night_cov','Albedo_BSA_nir','Albedo_BSA_shortwave','Albedo_WSA_nir','Albedo_WSA_shortwave','prcp', 'srad', 'dayl', 'tmax', 'tmin','LST_Day_1km','LST_Night_1km','Y','M', 'D', 'W', 'H', 'S','v_component_of_wind_10m','u_component_of_wind_10m','surface_pressure','runoff_sum','total_evaporation_sum','snowfall_sum','snowmelt_sum','snow_depth','snow_density','snow_cover','dewpoint_temperature_2m','temperature_2m'};

tmp_str = '';
for i = 1:length(features)
    eval(sprintf('mask_%s  = all(ismissing(%s), 2);',features{i},features{i}));
    tmp_str = strcat(tmp_str,['| mask_' (features{i})]);
end
% Combine masks using OR operation
eval(sprintf('combined_mask = %s;',string(tmp_str(3:end))))

% Get indices of rows to drop
indices_to_drop = find(combined_mask);

% Define the features again just in case
features = {'LC_Type1','LC_Type2','LC_Type3','LC_Type4','LC_Type5', 'avgLongitude', 'avgLatitude', 'Elevation', 'Slope', 'Aspect', 'NDVI', 'EVI', 'sur_refl_b01', 'sur_refl_b02', 'sur_refl_b03', 'sur_refl_b07', 'Clear_day_cov','Clear_night_cov','Albedo_BSA_nir','Albedo_BSA_shortwave','Albedo_WSA_nir','Albedo_WSA_shortwave','prcp', 'srad', 'dayl', 'tmax', 'tmin','LST_Day_1km','LST_Night_1km','Y','M', 'D', 'W', 'H', 'S','v_component_of_wind_10m','u_component_of_wind_10m','surface_pressure','runoff_sum','total_evaporation_sum','snowfall_sum','snowmelt_sum','snow_depth','snow_density','snow_cover','dewpoint_temperature_2m','temperature_2m'};

% Drop rows from each table and interpolate
for i = 1:length(features)
    feature = features{i};
    eval([feature ' = ' feature '(~combined_mask, :);']);  % Drop rows
    eval([feature ' = fillmissing(' feature ', ''linear'', ''EndValues'', ''nearest'');']);  % Linear interpolation
end

% Drop rows from Fire table
Fires = Fires(~combined_mask, :);
clear combined_mask mask_dewpoint_temperature_2m mask_runoff_sum mask_snow_cover mask_snow_density mask_snow_depth mask_snowfall_sum mask_snowmelt_sum mask_surface_pressure mask_temperature_2m mask_total_evaporation_sum mask_u_component_of_wind_10m mask_v_component_of_wind_10m mask_Albedo_BSA_nir mask_Albedo_BSA_shortwave mask_Albedo_WSA_nir mask_Albedo_WSA_shortwave mask_Clear_day_cov mask_Clear_night_cov mask_D mask_dayl mask_Elevation mask_EVI mask_H mask_LC_Type1 mask_LC_Type2 mask_LC_Type3 mask_LC_Type4 mask_LC_Type5 mask_LST_Day_1km mask_LST_Night_1km mask_M mask_NDVI mask_prcp mask_S mask_srad mask_sur_refl_b01 mask_sur_refl_b02 mask_sur_refl_b03 mask_sur_refl_b07 mask_tmax mask_tmin mask_W mask_Y

%%% Export all features for model learning in python
for i = 1:length(features)
    eval(sprintf('writetable(%s,''%s'');',features{i},strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Fires_task1_preprocessed_',features{i},'.csv')))
end
writetable(Fires,'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Fires_task1_preprocessed_targets.csv');
%% *1- Read all inputs and targets -V3

clear
clc

fprintf('Loading and processing the historical fires dataset (V3)...\n')

%%% load posclass data
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork',"Task_FireDetection_noRS.csv");
opts = detectImportOptions(path1);
opts.VariableTypes(2:3) = repmat({'datetime'},[1, length(opts.VariableNames(2:3))]);%convert all content to datetime
opts = setvaropts(opts, "DetectionDate_first", "InputFormat", "yyyy-dd-MM");
opts = setvaropts(opts, "DetectionDate_last", "InputFormat", "yyyy-dd-MM");
Fires_posClass = readtable(path1,opts);
Fires_posClass(Fires_posClass.isFire==0,:) = [];

%load negclass V3 data
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork',"Task_FireDetection_noRS_negClass_V3.csv");
opts = detectImportOptions(path1);
Fires_negClass = readtable(path1,opts);
Fires_negClass.isFire = zeros(height(Fires_negClass),1);

%Fuze targets
Fires = [Fires_posClass;Fires_negClass];
clear Fires_posClass Fires_negClass

%%% load input RS features
fprintf('Loading and processing the remote sensing features (V3)...\n')

features = {'LC_Type1','LC_Type2','LC_Type3','LC_Type4','LC_Type5', 'NDVI', 'EVI', 'sur_refl_b01', 'sur_refl_b02', 'sur_refl_b03', 'sur_refl_b07', 'Clear_day_cov','Clear_night_cov','Albedo_BSA_nir','Albedo_BSA_shortwave','Albedo_WSA_nir','Albedo_WSA_shortwave','prcp', 'srad', 'dayl', 'tmax', 'tmin','LST_Day_1km','LST_Night_1km','Y','M', 'D', 'W', 'H', 'S'};
for i = 1:length(features)
    %posClass
    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\',strcat('Fires_task1_',features{i}));
    opts = detectImportOptions(path1);
    eval(sprintf('%s  = readtable(path1,opts);',features{i}));
    %negClass
    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\',strcat('Fires_task1_negClassV3_',features{i}));
    opts = detectImportOptions(path1);
    eval(sprintf('%s_negClass  = readtable(path1,opts);',features{i}));
    %fuze data
    eval(sprintf('%s  = [%s(1:height(%s)/2,:);%s_negClass];',features{i},features{i},features{i},features{i}));
    eval(sprintf('clear %s_negClass',features{i}))
end
clear opts path1 i 

%%% load input Reanalysis features and topographic remote sensing features and coords.
fprintf('Loading and processing the Reanalysis features (V3)...\n')


features = {'avgLongitude', 'avgLatitude', 'Elevation', 'Slope', 'Aspect','v_component_of_wind_10m','u_component_of_wind_10m','surface_pressure','runoff_sum','total_evaporation_sum','snowfall_sum','snowmelt_sum','snow_depth','snow_density','snow_cover','dewpoint_temperature_2m','temperature_2m'};
for i = 1:length(features)
    %posClass
    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\',strcat('Fires_task1_',features{i}));
    opts = detectImportOptions(path1);
    eval(sprintf('%s  = readtable(path1,opts);',features{i}));
    %negClass
    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\',strcat('Fires_task1_negClassV3_',features{i}));
    opts = detectImportOptions(path1);
    eval(sprintf('%s_negClass  = readtable(path1,opts);',features{i}));
    %fuze data
    eval(sprintf('%s  = [%s(1:height(%s)/2,:);%s_negClass];',features{i},features{i},features{i},features{i}));
    eval(sprintf('clear %s_negClass',features{i}))
end
clear opts path1 i 

fprintf('Fuzing all data (V3)...\n')

% Create boolean masks for rows with all NaNs in each feature table (all RS and Reanalysis features)
features = {'LC_Type1','LC_Type2','LC_Type3','LC_Type4','LC_Type5', 'avgLongitude', 'avgLatitude', 'Elevation', 'Slope', 'Aspect', 'NDVI', 'EVI', 'sur_refl_b01', 'sur_refl_b02', 'sur_refl_b03', 'sur_refl_b07', 'Clear_day_cov','Clear_night_cov','Albedo_BSA_nir','Albedo_BSA_shortwave','Albedo_WSA_nir','Albedo_WSA_shortwave','prcp', 'srad', 'dayl', 'tmax', 'tmin','LST_Day_1km','LST_Night_1km','Y','M', 'D', 'W', 'H', 'S','v_component_of_wind_10m','u_component_of_wind_10m','surface_pressure','runoff_sum','total_evaporation_sum','snowfall_sum','snowmelt_sum','snow_depth','snow_density','snow_cover','dewpoint_temperature_2m','temperature_2m'};

tmp_str = '';
for i = 1:length(features)
    eval(sprintf('mask_%s  = all(ismissing(%s), 2);',features{i},features{i}));
    tmp_str = strcat(tmp_str,['| mask_' (features{i})]);
end
% Combine masks using OR operation
eval(sprintf('combined_mask = %s;',string(tmp_str(3:end))))

% Get indices of rows to drop
indices_to_drop = find(combined_mask);

% Define the features again just in case
features = {'LC_Type1','LC_Type2','LC_Type3','LC_Type4','LC_Type5', 'avgLongitude', 'avgLatitude', 'Elevation', 'Slope', 'Aspect', 'NDVI', 'EVI', 'sur_refl_b01', 'sur_refl_b02', 'sur_refl_b03', 'sur_refl_b07', 'Clear_day_cov','Clear_night_cov','Albedo_BSA_nir','Albedo_BSA_shortwave','Albedo_WSA_nir','Albedo_WSA_shortwave','prcp', 'srad', 'dayl', 'tmax', 'tmin','LST_Day_1km','LST_Night_1km','Y','M', 'D', 'W', 'H', 'S','v_component_of_wind_10m','u_component_of_wind_10m','surface_pressure','runoff_sum','total_evaporation_sum','snowfall_sum','snowmelt_sum','snow_depth','snow_density','snow_cover','dewpoint_temperature_2m','temperature_2m'};

% Drop rows from each table and interpolate
for i = 1:length(features)
    feature = features{i};
    eval([feature ' = ' feature '(~combined_mask, :);']);  % Drop rows
    eval([feature ' = fillmissing(' feature ', ''linear'', ''EndValues'', ''nearest'');']);  % Linear interpolation
end

% Drop rows from Fire table
Fires = Fires(~combined_mask, :);
clear combined_mask mask_dewpoint_temperature_2m mask_runoff_sum mask_snow_cover mask_snow_density mask_snow_depth mask_snowfall_sum mask_snowmelt_sum mask_surface_pressure mask_temperature_2m mask_total_evaporation_sum mask_u_component_of_wind_10m mask_v_component_of_wind_10m mask_Albedo_BSA_nir mask_Albedo_BSA_shortwave mask_Albedo_WSA_nir mask_Albedo_WSA_shortwave mask_Clear_day_cov mask_Clear_night_cov mask_D mask_dayl mask_Elevation mask_EVI mask_H mask_LC_Type1 mask_LC_Type2 mask_LC_Type3 mask_LC_Type4 mask_LC_Type5 mask_LST_Day_1km mask_LST_Night_1km mask_M mask_NDVI mask_prcp mask_S mask_srad mask_sur_refl_b01 mask_sur_refl_b02 mask_sur_refl_b03 mask_sur_refl_b07 mask_tmax mask_tmin mask_W mask_Y

%%% Export all features for model learning in python
for i = 1:length(features)
    eval(sprintf('writetable(%s,''%s'');',features{i},strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Fires_task1_preprocessed_V3_',features{i},'.csv')))
end
writetable(Fires,'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Fires_task1_preprocessed_V3_targets.csv');

%% *1- Read all inputs and targets -V4

clear
clc

fprintf('Loading and processing the historical fires dataset (V4)...\n')

%%% load posclass data
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork',"Task_FireDetection_noRS.csv");
opts = detectImportOptions(path1);
opts.VariableTypes(2:3) = repmat({'datetime'},[1, length(opts.VariableNames(2:3))]);%convert all content to datetime
opts = setvaropts(opts, "DetectionDate_first", "InputFormat", "yyyy-dd-MM");
opts = setvaropts(opts, "DetectionDate_last", "InputFormat", "yyyy-dd-MM");
Fires_posClass = readtable(path1,opts);
Fires_posClass(Fires_posClass.isFire==0,:) = [];

%load negclass V4 data
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork',"Task_FireDetection_noRS_negClass_V4.csv");
opts = detectImportOptions(path1);
Fires_negClass = readtable(path1,opts);
Fires_negClass.isFire = zeros(height(Fires_negClass),1);

%Fuze targets
Fires = [Fires_posClass;Fires_negClass];
clear Fires_posClass Fires_negClass

%%% load input RS features
fprintf('Loading and processing the remote sensing features (V4)...\n')

features = {'LC_Type1','LC_Type2','LC_Type3','LC_Type4','LC_Type5', 'NDVI', 'EVI', 'sur_refl_b01', 'sur_refl_b02', 'sur_refl_b03', 'sur_refl_b07', 'Clear_day_cov','Clear_night_cov','Albedo_BSA_nir','Albedo_BSA_shortwave','Albedo_WSA_nir','Albedo_WSA_shortwave','prcp', 'srad', 'dayl', 'tmax', 'tmin','LST_Day_1km','LST_Night_1km','Y','M', 'D', 'W', 'H', 'S'};
for i = 1:length(features)
    %posClass
    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\',strcat('Fires_task1_',features{i}));
    opts = detectImportOptions(path1);
    eval(sprintf('%s  = readtable(path1,opts);',features{i}));
    %negClass
    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\',strcat('Fires_task1_negClassV4_',features{i}));
    opts = detectImportOptions(path1);
    eval(sprintf('%s_negClass  = readtable(path1,opts);',features{i}));
    %fuze data
    eval(sprintf('%s  = [%s(1:height(%s)/2,:);%s_negClass];',features{i},features{i},features{i},features{i}));
    eval(sprintf('clear %s_negClass',features{i}))
end
clear opts path1 i 

%%% load input Reanalysis features and topographic remote sensing features and coords.
fprintf('Loading and processing the Reanalysis features (V4)...\n')


features = {'avgLongitude', 'avgLatitude', 'Elevation', 'Slope', 'Aspect','v_component_of_wind_10m','u_component_of_wind_10m','surface_pressure','runoff_sum','total_evaporation_sum','snowfall_sum','snowmelt_sum','snow_depth','snow_density','snow_cover','dewpoint_temperature_2m','temperature_2m'};
for i = 1:length(features)
    %posClass
    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\',strcat('Fires_task1_',features{i}));
    opts = detectImportOptions(path1);
    eval(sprintf('%s  = readtable(path1,opts);',features{i}));
    %negClass
    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\',strcat('Fires_task1_negClassV4_',features{i}));
    opts = detectImportOptions(path1);
    eval(sprintf('%s_negClass  = readtable(path1,opts);',features{i}));
    %fuze data
    eval(sprintf('%s  = [%s(1:height(%s)/2,:);%s_negClass];',features{i},features{i},features{i},features{i}));
    eval(sprintf('clear %s_negClass',features{i}))
end
clear opts path1 i 

fprintf('Fuzing all data (V4)...\n')

% Create boolean masks for rows with all NaNs in each feature table (all RS and Reanalysis features)
features = {'LC_Type1','LC_Type2','LC_Type3','LC_Type4','LC_Type5', 'avgLongitude', 'avgLatitude', 'Elevation', 'Slope', 'Aspect', 'NDVI', 'EVI', 'sur_refl_b01', 'sur_refl_b02', 'sur_refl_b03', 'sur_refl_b07', 'Clear_day_cov','Clear_night_cov','Albedo_BSA_nir','Albedo_BSA_shortwave','Albedo_WSA_nir','Albedo_WSA_shortwave','prcp', 'srad', 'dayl', 'tmax', 'tmin','LST_Day_1km','LST_Night_1km','Y','M', 'D', 'W', 'H', 'S','v_component_of_wind_10m','u_component_of_wind_10m','surface_pressure','runoff_sum','total_evaporation_sum','snowfall_sum','snowmelt_sum','snow_depth','snow_density','snow_cover','dewpoint_temperature_2m','temperature_2m'};

tmp_str = '';
for i = 1:length(features)
    eval(sprintf('mask_%s  = all(ismissing(%s), 2);',features{i},features{i}));
    tmp_str = strcat(tmp_str,['| mask_' (features{i})]);
end
% Combine masks using OR operation
eval(sprintf('combined_mask = %s;',string(tmp_str(3:end))))

% Get indices of rows to drop
indices_to_drop = find(combined_mask);

% Define the features again just in case
features = {'LC_Type1','LC_Type2','LC_Type3','LC_Type4','LC_Type5', 'avgLongitude', 'avgLatitude', 'Elevation', 'Slope', 'Aspect', 'NDVI', 'EVI', 'sur_refl_b01', 'sur_refl_b02', 'sur_refl_b03', 'sur_refl_b07', 'Clear_day_cov','Clear_night_cov','Albedo_BSA_nir','Albedo_BSA_shortwave','Albedo_WSA_nir','Albedo_WSA_shortwave','prcp', 'srad', 'dayl', 'tmax', 'tmin','LST_Day_1km','LST_Night_1km','Y','M', 'D', 'W', 'H', 'S','v_component_of_wind_10m','u_component_of_wind_10m','surface_pressure','runoff_sum','total_evaporation_sum','snowfall_sum','snowmelt_sum','snow_depth','snow_density','snow_cover','dewpoint_temperature_2m','temperature_2m'};

% Drop rows from each table and interpolate
for i = 1:length(features)
    feature = features{i};
    eval([feature ' = ' feature '(~combined_mask, :);']);  % Drop rows
    eval([feature ' = fillmissing(' feature ', ''linear'', ''EndValues'', ''nearest'');']);  % Linear interpolation
end

% Drop rows from Fire table
Fires = Fires(~combined_mask, :);
clear combined_mask mask_dewpoint_temperature_2m mask_runoff_sum mask_snow_cover mask_snow_density mask_snow_depth mask_snowfall_sum mask_snowmelt_sum mask_surface_pressure mask_temperature_2m mask_total_evaporation_sum mask_u_component_of_wind_10m mask_v_component_of_wind_10m mask_Albedo_BSA_nir mask_Albedo_BSA_shortwave mask_Albedo_WSA_nir mask_Albedo_WSA_shortwave mask_Clear_day_cov mask_Clear_night_cov mask_D mask_dayl mask_Elevation mask_EVI mask_H mask_LC_Type1 mask_LC_Type2 mask_LC_Type3 mask_LC_Type4 mask_LC_Type5 mask_LST_Day_1km mask_LST_Night_1km mask_M mask_NDVI mask_prcp mask_S mask_srad mask_sur_refl_b01 mask_sur_refl_b02 mask_sur_refl_b03 mask_sur_refl_b07 mask_tmax mask_tmin mask_W mask_Y

%%% Export all features for model learning in python
for i = 1:length(features)
    eval(sprintf('writetable(%s,''%s'');',features{i},strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Fires_task1_preprocessed_V4_',features{i},'.csv')))
end
writetable(Fires,'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\Fires_task1_preprocessed_V4_targets.csv');

%% 2-1 Create task1 dataset - hrz=+days  (best test perf: 98%)
clc
for i = [0 7 30 90]
    eval(sprintf('clear Fire%dDay',i))
    % length = 1
    for j = 1:length(features)
        eval(sprintf('Fire%dDay.%s = %s{:,end-i};',i,features{j},features{j}));
    end
    eval(sprintf('Fire%dDay = struct2table(Fire%dDay(:));',i,i))
    eval(sprintf('Fire%dDay.isFire = Fires{:,end};',i))

    eval(sprintf('writetable(Fire%dDay,''%s'');',i,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\FireDetection_datasetV4_',num2str(i),'Days.csv')))
end

% Make time series inputs for deep learning time series benchamrk tasks
clc
selectedrows = zeros(size(LC_Type1));
for k = 1:height(training)
    eval(sprintf('selectedrows(k,:) = LC_Type1{training.Var1(k)+1,:};',features{1}));
end

%% 3-1 Create task2 dataset - hrz=+days  (best test perf: ?%)
clc
for i = [0 7 30 90]
    eval(sprintf('clear FireBurnt%dDay',i))
    % length = 1
    for j = 1:length(features)
        eval(sprintf('FireBurnt%dDay.%s = %s{:,end-i};',i,features{j},features{j}));
    end
    eval(sprintf('FireBurnt%dDay = struct2table(FireBurnt%dDay(:));',i,i))
    eval(sprintf('FireBurnt%dDay.totalArea_km2 = Fires{:,"totalArea_km2"};',i))
    [indices,~] = find(Fires.isFire==1);
    eval(sprintf('FireBurnt%dDay = FireBurnt%dDay(indices,:);',i,i))
    
    eval(sprintf('writetable(FireBurnt%dDay,''%s'');',i,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\FireBurntArea_dataset_',num2str(i),'Days.csv')))
end

%% 4-1 Create task3 dataset - hrz=+days  (best test perf: test RMSE = ?)

for i = [0 7 30 90]
    eval(sprintf('clear FireDuration%dDay',i))
    for j = 1:length(features)
        eval(sprintf('FireDuration%dDay.%s = %s{:,end-i};',i,features{j},features{j}));
    end
    eval(sprintf('FireDuration%dDay = struct2table(FireDuration%dDay(:));',i,i))
    eval(sprintf('FireDuration%dDay.Duration = Fires{:,"Duration"};',i))
    [indices,~] = find(Fires.isFire==1);
    eval(sprintf('FireDuration%dDay = FireDuration%dDay(indices,:);',i,i))

    eval(sprintf('writetable(FireDuration%dDay,''%s'');',i,strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2\FireDuration_dataset_',num2str(i),'Days.csv')))
end

%% Create task3 dataset - length=30  (best test perf: test RMSE = 8.06)

clear FireDuration
% length = 1
for i = 1:length(features)
    eval(sprintf('FireDuration.%s = %s{:,end-31:end}',features{i},features{i}));
end
FireDuration = struct2table(FireDuration(:));
FireDuration.isFire = Fires{:,end};    
FireDuration.Duration = Fires{:,"Duration"};    

FireDuration = FireDuration(FireDuration.isFire==1,:);
FireDuration = removevars(FireDuration, "isFire");
























%% Compute total wildfire area 
clc
% Spatial resolution
pixelArea_km2 = 0.25; % Area in square kilometers for a 500m resolution

% Initialize variables for cluster statistics
uniqueClusters = unique(fireData.Cluster);
clusterAreas_km2 = zeros(length(uniqueClusters), 1);
clusterBrightness = zeros(length(uniqueClusters), 1);
clusterBright_t31 = zeros(length(uniqueClusters), 1);
clusterFRP = zeros(length(uniqueClusters), 1);
clusterConfidence = zeros(length(uniqueClusters), 1);
clusterDate_first = NaT(length(uniqueClusters), 1);
clusterDate_first.Format = "uuuu-dd-MM";
centroidLats = zeros(length(uniqueClusters), 1);
centroidLons = zeros(length(uniqueClusters), 1);

clusterPreNDVI = zeros(length(uniqueClusters), 1);
clusterPostNDVI = zeros(length(uniqueClusters), 1);

% Create a WGS84 spheroid object
wgs84 = wgs84Ellipsoid('kilometer');

% Iterate over each cluster to calculate the area
for i = 1:length(uniqueClusters)
    clusterNum = uniqueClusters(i);

    % Extract the latitude and longitude of the current cluster
    clusterIndices = fireData.Cluster == clusterNum;
    lats = fireData.latitude(clusterIndices);
    lons = fireData.longitude(clusterIndices);
    
    % Extract mean remote sensing and other features
    clusterBrightness(i) = max(fireData.brightness(clusterIndices));
    clusterBright_t31(i) = max(fireData.bright_t31(clusterIndices));
    clusterFRP(i) = max(fireData.frp(clusterIndices));
    clusterConfidence(i) = mean(fireData.confidence(clusterIndices));

    clusterPreNDVI(i) = mean(fireData.preNDVI(clusterIndices));
    clusterPostNDVI(i) = mean(fireData.postNDVI(clusterIndices));

    % Get date
    t = mean(datetime(year(fireData.DetectionTime(clusterIndices)), 1, 1, 'Format', 'yyyy-dd-MM'));  % First day of the year
    clusterDate_first(i) = t + days(mean(week(fireData.DetectionTime(clusterIndices))-1).*7);

    % Calculate the mean latitude and longitude
    centroidLats(i) = mean(lats);
    centroidLons(i) = mean(lons);

    % Check the number of points in the cluster
    numPoints = length(lats);
    if numPoints == 1 | (numPoints == 2 & lats(1)-lats(2) ==0)
        % Single point represents one pixel area
        clusterAreas_km2(i) = pixelArea_km2;
    elseif numPoints == 2 & lats(1)-lats(2) ~=0
        clusterAreas_km2(i) = estimateAreaFromTwoPoints(lats, lons);
    else
        % Convert lat/lon to Cartesian coordinates
        [x, y, z] = geodetic2ecef(wgs84, lats, lons, zeros(size(lats)), 'degrees');
        try
            DT = delaunayTriangulation(x,y);
            C = convexHull(DT);
            clusterAreas_km2(i) = polyarea(DT.Points(C,1), DT.Points(C,2));
            % [k, area_km2] = convhull(x, y);
        catch
            % Handle cases where convex hull computation fails
            clusterAreas_km2(i) = NaN;
        end
    end
end

% Create a table with cluster statistics
clusterStats = table(uniqueClusters, centroidLats, centroidLons, clusterAreas_km2, clusterBrightness, clusterBright_t31, clusterFRP, clusterConfidence, clusterDate_first, clusterPreNDVI, clusterPostNDVI, clusterPostNDVI-clusterPreNDVI,...
    'VariableNames', {'ClusterID', 'avgLatitude','avgLongitude','totalArea_km2', 'maxBrightness', 'maxBright_t31', 'maxFRP', 'avgConfidence', 'Date', 'preNDVI', 'postNDVI', 'dNDVI'});

% Display the cluster statistics
% clusterStats(1:52,:) = [];

%% Plot total area of wildfires
clc
size1 = 30;

figure('units','normalized','outerposition',[0 0 1 1])
% t = tiledlayout(1,1,'TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);
% nexttile,
s = geobubble(clusterStats,"avgLatitude","avgLongitude","SizeVariable","totalArea_km2","Basemap","landcover");
set(gca,'FontSize',20)
s.SizeLegendTitle = 'Total Area (km2)';
%% plot clusters of wildfires
clc
size1 = 30;

figure('units','normalized','outerposition',[0 0 1 1])
% t = tiledlayout(1,1,'TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);
% nexttile,
geoaxes("Basemap","landcover")
s = geoscatter(clusterStats,"avgLatitude","avgLongitude","filled");
s.SizeData = 30;
s.ColorVariable = "ClusterID";
c = colorbar;
c.Label.String = "Temporal-spatio clusters (centroids)";

set(gca,'FontSize',20)
colormap turbo
% exportgraphics(gcf,'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month17\Fire\Figure_densityPlot.jpg',"ContentType","image")

%% plot density of wildfires
size1 = 30;

figure('units','normalized','outerposition',[0 0 1 1])
% t = tiledlayout(1,1,'TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);
% nexttile,
geoaxes("Basemap","landcover")
dp = geodensityplot(clusterStats.avgLatitude,clusterStats.avgLongitude);
dp.FaceColor = 'interp';
colormap turbo
set(gca,'FontSize',size1)
% exportgraphics(gcf,'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month17\Fire\Figure_densityPlot.jpg',"ContentType","image")

%% Compute index 1: wildfire intensity
size1 = 30;

% Find the minimum and maximum values of the array
minValue = min(clusterStats.dNDVI);
maxValue = max(clusterStats.dNDVI);
% Perform min-max normalization
scaledNDVI = (clusterStats.dNDVI - minValue) / (maxValue - minValue);

clusterStats.intensity_index = scaledNDVI;

figure('units','normalized','outerposition',[0 0 1 1])
% t = tiledlayout(1,1,'TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);
% nexttile,
geoaxes("Basemap","landcover")

s = geoscatter(clusterStats,"avgLatitude","avgLongitude","filled");
s.SizeData = 36;
s.ColorVariable = "intensity_index";
c = colorbar();
colormap turbo
c.Label.String = "Intensity Index";
c.Label.Interpreter = 'latex';
set(gca,'FontSize',size1)
% exportgraphics(gcf,'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Month17\Fire\Figure_IntensityIndex.jpg',"ContentType","image")


%% Functions
function area_km2 = estimateAreaFromTwoPoints(lats, lons)
    % Earth's radius in kilometers
    R = 6371; 

    % Convert degrees to radians
    lat1 = deg2rad(lats(1,1));
    lon1 = deg2rad(lons(1,1));
    lat2 = deg2rad(lats(2,1));
    lon2 = deg2rad(lons(2,1));

    % Haversine distance in latitude
    deltaLat = abs(lat2 - lat1);
    a = sin(deltaLat / 2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1 - a));
    distLat_km = R * c;

    % Haversine distance in longitude
    deltaLon = abs(lon2 - lon1);
    a = sin(deltaLon / 2)^2 + cos(lat1) * cos(lat2) * sin(deltaLon / 2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1 - a));
    distLon_km = R * c;

    % Estimate area as product of sides (approximation)
    area_km2 = distLat_km * distLon_km;
end

function [Y,M,D] = get_YMD(x)
    [Y,M,D] = ymd(x);
end

function seasons = get_astronomical_season(dates)
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
end

function degs = distance2lonChange(distKM,latitude)
    degs = distKM/(111.320*cos(latitude));
end

function change = randchangeLat(distKM)
    switch randi([1, 2])
        case 1
            change = round(distKM/110.574,4);
        case 2
            change = round(distKM/110.574,4);
    end
end


function change = randchangeLon(distKM, latitude)
    switch randi([1, 2])
        case 1
            change = round(distance2lonChange(distKM,latitude),3);
        case 2
            change = -round(distance2lonChange(distKM,latitude),3);
    end
end