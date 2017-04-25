function PlotBaselines( filename, master_yyyymmdd, name, switchxy )
%% Plot the baseline of a .cst file
% This program is intendet to plot the baseline configuration
% usage: PLOTBASELINES( Baselines.cst, 2016-02-19 )
% 

switch nargin
    case 1
        error('Specify master date as second input argument ''yyyymmdd''')
    case 2
        warning('No name specified for plot')
        switchxy = 0;
        name = ' ';
    case 3
        switchxy = 0;
    case 4
        if switchxy == 1
            fprintf('Switching temporal and perpendicular baselines\n')
        elseif switchxy == 0
            fprintf('Not switching temporal and perpendicular baselines\n')
        end
    otherwise
        warning('wrong number of input arguments')
end


%% Initialize variables.
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

fclose(fileID);

%% Allocate imported array to column variable names
Date = dataArray{:, 1};
%Bpar = dataArray{:, 2};
Bper = dataArray{:, 3};
%Base = dataArray{:, 4};
%Bhor = dataArray{:, 5};
%Bver = dataArray{:, 6};
%Alpha = dataArray{:, 7};
clearvars fileID delimiter startRow formatSpec fileID dataArray ans;

%%
dates = num2str(Date);
dates = datenum(dates,'yyyymmdd');
mdate = datenum(master_yyyymmdd, 'yyyymmdd');
dates(:,2) = mdate;
Bper(:,2)  = 0;
%%
% figure(1)
% dates(:,2) = datenum(2016,02,19);
% Bper(:,2)  = 0;
% m1 = plot(dates', Bper', 'O', 'markersize', 10);
% hold on
% l1 = plot(dates', Bper', '-', 'linewidth', 2, 'Color', [1 0.5 0] );
% datetick('x',12)

figure()
hold on
dates = dates - mdate;
switch switchxy
    case 1
        plot(dates', Bper', 'O', 'markersize', 10,'linewidth',1.5);
        plot(dates', Bper', '-', 'linewidth', 2, 'Color', [1 0.5 0] );
        plot(0,0,'O','markersize',15,'markerfacecolor',[0 0 0]);
        text(20,0,'Master')
        xlabel('Temporal Baseline [days]')
        ylabel('Perpendicular Baseline [m]')
    case 0
        plot(Bper', dates', 'O', 'markersize', 10,'linewidth',1.5);
        plot(Bper', dates', '-', 'linewidth', 2, 'Color', [1 0.5 0] );
        text(10,0,'Master')
        ylabel('Temporal Baseline [days]')
        xlabel('Perpendicular Baseline [m]')
        plot(0,0,'O','markersize',15,'markerfacecolor',[0 0 0]);
        xlim([round(min(Bper(:,1))-10,-1), round(max(Bper(:,1))+10,-1)])
    otherwise
        error('swicthxy flag can be set to 0 or 1')
end

title(['Baseline configuration: ',name])
set(gcf, 'Units','Normalized','OuterPosition',[0.2 0.2 0.30 0.7])
set(gca,'DataAspectRatio',[1,2.5,1])


print(filename(1:end-4),'-dpng','-r300')

end