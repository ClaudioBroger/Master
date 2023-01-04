function [] = plot_OutV(OutV,paramSpecs,fileName)

%% plot cost

if ~isfield(OutV,'colnames')
    OutV.colnames = paramSpecs.names;
end
NotProjectedIdx = find(ismember(paramSpecs.names,OutV.colnames));
% paraoptestimated = paraopt(NotProjectedIdx); % to add fitted parameter point


paramSpecsestimated = paramSpecs(NotProjectedIdx,:);
data = OutV;
if isfield(OutV,'V'),
data.rowmat = data.V;
end
paramNames = paramSpecsestimated.names;
bmin = paramSpecsestimated.bmin;
bmax = paramSpecsestimated.bmax;
bmin(paramSpecsestimated.islog==1) = log10(bmin(paramSpecsestimated.islog==1));
bmax(paramSpecsestimated.islog==1) = log10(bmax(paramSpecsestimated.islog==1));
paramtoplot = paramNames;
n = length(paramtoplot);

indices = zeros(1,n);
x = zeros(length(data.rowmat(:,1)),n);
xmin = zeros(1,n);
xmax = zeros(1,n);
for k =1:n
indices(k) = find(ismember(paramNames,paramtoplot{k}));
x(:,k) =  data.rowmat(:,k);
xmin(k) = bmin(indices(k));
xmax(k) = bmax(indices(k));
end


% f = figure();

% clf;


factor = 1.25;
map = copper();
for k=1:(n-1),
    subplot(n,n,(k-1)*n+1)
    [~,density,X,Y] = kde2d([x(:,1),x(:,k)],10);
    contour(X,Y,density,'LineWidth',1);
    hold on;
%         plot(log10(paraoptestimated(1)),log10(paraoptestimated(k)),'or','LineWidth',2); % to add fitted parameter point
    pos = get(gca, 'Position');
    pos(3) = 0.7*factor*pos(3);
    pos(4) = factor*pos(4);
    set(gca, 'Position', pos, 'XTick',[]);       
    axis([xmin(1) xmax(1) xmin(k) xmax(k)]);
    for l =2:(k-1),
        subplot(n,n,(k-1)*n+l)
        [~,density,X,Y] = kde2d([x(:,l),x(:,k)],10);
        contour(X,Y,density,'LineWidth',1);
        hold on;
%             plot(log10(paraoptestimated(l)),log10(paraoptestimated(k)),'or','LineWidth',2); % to add fitted parameter point
        pos = get(gca, 'Position');
        pos(3) = 0.7*factor*pos(3);
        pos(4) = factor*pos(4);
        set(gca, 'Position', pos, 'XTick',[],'YTick',[]);       
        axis([xmin(l) xmax(l) xmin(k) xmax(k)]);   
    end
    subplot(n,n,(k-1)*n+k)
    pd = fitdist(x(:,k),'Kernel');
    X = xmin(k):0.1:xmax(k);
    Y = pdf(pd,X);
    plot(X,Y,'Color','black','LineWidth',1);
    title(paramtoplot{k});
    pos = get(gca, 'Position');
    pos(3) = 0.7*factor*pos(3);
    pos(4) = factor*pos(4);
    set(gca, 'Position', pos, 'XTick',[],'YTick',[]);       
    axis([[xmin(k) xmax(k)] 0 max(Y)]);       
end

k = n;
subplot(n,n,(k-1)*n+1)
[~,density,X,Y] = kde2d([x(:,1),x(:,k)],10);
contour(X,Y,density,'LineWidth',1);
hold on;
%     plot(log10(paraoptestimated(1)),log10(paraoptestimated(k)),'or','LineWidth',2); % to add fitted parameter point
pos = get(gca, 'Position');
pos(3) = 0.7*factor*pos(3);
pos(4) = factor*pos(4);
set(gca, 'Position', pos); 
axis([xmin(1) xmax(1) xmin(k) xmax(k)]);

for l =2:(k-1),
    subplot(n,n,(k-1)*n+l)
    [~,density,X,Y] = kde2d([x(:,l),x(:,k)],10);
    contour(X,Y,density,'LineWidth',1);
    hold on;
%         plot(log10(paraoptestimated(l)),log10(paraoptestimated(k)),'or','LineWidth',2); % to add fitted parameter point
    pos = get(gca, 'Position');
    pos(3) = 0.7*factor*pos(3);
    pos(4) = factor*pos(4);
    set(gca, 'Position', pos, 'YTick',[]);      
    axis([xmin(l) xmax(l) xmin(k) xmax(k)]);    
end
l = k;
subplot(n,n,(k-1)*n+l)
pd = fitdist(x(:,k),'Kernel');
X = xmin(k):0.1:xmax(k);
Y = pdf(pd,X);
plot(X,Y,'Color','black','LineWidth',1);
title(paramtoplot{k});
pos = get(gca, 'Position');
pos(3) = 0.7*factor*pos(3);
pos(4) = factor*pos(4);
set(gca, 'Position', pos, 'YTick',[]);      
axis([[xmin(k) xmax(k)] 0 max(Y)]); 


colormap(map);
%     c = colorbar();

%     pos = get(c, 'Position');
%     set (c, 'Position', [pos(1) pos(2)+0.5 pos(3)*3 pos(4)*3]);
%     c.Label.String = 'Probability';

if nargin == 3,    
    set(f,'PaperOrientation','landscape','PaperUnits','Normalized','PaperPosition',[0 0 1 1]);
    print(fileName, '-dpdf');
end
end