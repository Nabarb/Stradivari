%% stradivari - plots a combination of half-violins, and raw datapoints (1d or 2d scatter).
% [h,u] = stradivari(X,Name,Value)
% displays the violin plots of X. X is a data matrix with columns 
% corresponding to data points and rows to different cases. 
% Each case gets its own violin.

% [h,u] = stradivari(AX,...)
% plots into the axes with handle AX. If no axes are given it creates a new
% figure.
%
% Seee below for optional inputs.
%
%
% heavilly based on raincloud_plot,
% https://github.com/RainCloudPlots/RainCloudPlots.
%
% (C) Federico Barban 2019, modification of RaincloudPlot to add 
% functionalities and improve flexibility. 
% Released under GNU General Public License. 

%%
% ---GNU General Public License Copyright---
% This file is part of DataHigh.
% 
% raincloud_plot is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, version 2.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details in COPYING.txt found
% in the main DataHigh directory.
% 
% You should have received a copy of the GNU General Public License
% along with raincloud_plot.  If not, see <http://www.gnu.org/licenses/>.



function [h, u] = stradivari(varargin)

% ---------------------------- INPUT ----------------------------
%
% X                  - vector of data to be plotted, required.
%
% --------------------- OPTIONAL ARGUMENTS ----------------------
%
% color              - color vector for rainclouds (default color gradient
%                      defined below in MYColors). Can either be a single
%                      vector or a matrix M x 3, where M is the first
%                      dimension of X.
% band_width         - band_width of smoothing kernel (default = 1)
% density_type       - choice of density algo ('ks' or 'rath'). Default = 'ks'
% box_on             - logical to turn box plots on/off (default = 0)
% alpha              - scalar positive value to increase cloud alpha (defalut = 1)
% dot_dodge_amount   - scalar value to increase dot dodge amounts (defalut =0.6)
% line_width         - scalar value to set global line width of clouds (default = 2)
% scatter_size       - Matrix the same size as X that defines the size of
%                      each marker in the scatterplot
% scatter_width      - scalar value that defines the percentual amount of
%                      cloud area occupied by the scatter
% vertical           - 1 or zero. Plot clouds vertical or horizzontal
% coupled            - Matrix. determines which rows of X to couple in a 
%                      full violin. Row numbers along its first dimension 
%                      get coupled  Ex:
%                      X = rand(5,20); 
%                      ind = [ 1 3 nan;...
%                              2 4 5]
%                      raincloud_plot(X,'coupled',ind);
%
%                     This would plot rows 1,2 and 3,4 together and 5
%                     alone.
% jitter              - Same size as X. Used to encode a second dimension
%                       in the scatter.
% grid_on             - 1 or zero. Plots small markers on the side
%                       corresponding to the median value
%
% ---------------------------- OUTPUT ----------------------------
% h - figure handle to change more stuff
% u - parameter from kernel density estimate
%


%% check all the inputs and if they do not exist then revert to default settings
% input parsing settings
p = inputParser;
p.CaseSensitive = true;
p.Parameters;
p.Results;
p.KeepUnmatched = true;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);

% set the desired and optional input arguments


if  isgraphics(varargin{1})    
    ax = varargin{1};
    varargin(1) = [];
end

if isnumeric(varargin{1})
    X = varargin{1};
    varargin(1) = [];
end

if ~exist('ax','var')
    ax = axes(figure);
end

if ~exist('X','var')
    error('No input data provided!');
end

addOptional(p, 'color', MyColors(1:size(X,1)), @isnumeric)
addOptional(p, 'band_width', [])
addOptional(p, 'density_type', 'ks', @ischar)
addOptional(p, 'box_on', 0, @isnumeric)
addOptional(p, 'box_dodge', 0, @isnumeric)
addOptional(p, 'box_dodge_amount', 0, @isnumeric)
addOptional(p, 'alpha', 0.6, validScalarPosNum)
addOptional(p, 'dot_dodge_amount', 0.5, @isnumeric)
addOptional(p, 'box_col_match', 0, @isnumeric)
addOptional(p, 'line_width', 0.01, validScalarPosNum)
addOptional(p, 'lwr_bnd', 1, @isnumeric)
addOptional(p, 'bxcl', [0 0 0], @isnumeric)
addOptional(p, 'bxfacecl', [1 1 1], @isnumeric)
addOptional(p, 'cloud_edge_col', [0 0 0], @isnumeric)
addOptional(p, 'vertical', 0, @isnumeric);
addOptional(p, 'scatter_size', 10, @isnumeric);
addOptional(p, 'scatter_width', 0.25, @isnumeric);
addOptional(p, 'grid_on', 0, @isnumeric)

[M,N]=size(X);
addOptional(p, 'coupled', (1:M), @isnumeric);
addOptional(p, 'jitter',rand(M,N), @isnumeric)


parse(p,varargin{:});
% then set/get all the inputs out of this structure
color               = p.Results.color;
if isvector( color), color = (color(:)*ones(1,M))';end
density_type        = p.Results.density_type;
box_on              = p.Results.box_on;
% box_dodge           = p.Results.box_dodge;
% box_dodge_amount    = p.Results.box_dodge_amount;
alpha               = p.Results.alpha;
dot_dodge_amount    = p.Results.dot_dodge_amount;
line_width          = p.Results.line_width;
% lwr_bnd             = p.Results.lwr_bnd;
% bxcl                = p.Results.bxcl;
% bxfacecl            = p.Results.bxfacecl;
% cloud_edge_col      = p.Results.cloud_edge_col;
band_width          = p.Results.band_width;
scatter_size        = p.Results.scatter_size;
vertical            = p.Results.vertical;
coupled             = p.Results.coupled;
jitter              = p.Results.jitter;
scatter_width       = p.Results.scatter_width;
grid_on             = p.Results.grid_on;

% opts = p.Results;
if M>N
 warning('Looks like there are more cases than observations!');
end

if isscalar(scatter_size)
    scatter_size = ones(M,N) * scatter_size;
elseif any(size(scatter_size) ~= [M,N])
    error('scatter_size size does not match input size.');    
end

if any(size(jitter) ~= [M,N])
    error('jitter size does not match input size.');    
elseif min(jitter(:))<=0 || max(jitter(:))>1 
    error('jitter must be (0 1] range');
end

if vertical
    AxDir = 'XLim';
    
else
    AxDir = 'YLim';
    
end
ind = coupled;
ax = gca;

if exist('rst_RASH', 'file') ~= 2
    RSTpresent = false;
elseif exist('rst_RASH', 'file') == 2
    RSTpresent = true;
end

%% Compute distribution and jittering
for ii=1:size(ind,1)
    for jj=1:size(ind,2)
        if(isnan(ind(ii,jj))),Ydist{ii,jj}=0;continue;end
        %% calculate kernel density
        
        if RSTpresent,outliers{ii,jj} = rst_outlier(X(ind(ii,jj),:),3);
        else,         outliers{ii,jj} = isoutlier(X(ind(ii,jj),:),'quartiles');end
        
        switch density_type
            case 'ks'
                [Ydist{ii,jj}, Xdist{ii,jj}, u] = ksdensity(X(ind(ii,jj),~outliers{ii,jj}), 'bandwidth', band_width);
                
            case 'rash'
                if  RSTpresent
                    % must have https://github.com/CPernet/Robust_Statistical_Toolbox
                    % for this to work
                    
                    % check for rst_RASH function (from Robust stats toolbox) in path, fail if not found
                    error('Could not compute density using RASH method. \nDo you have the Robust Stats toolbox on your path?');
                end
                [Xdist{ii,jj}, Ydist{ii,jj}] = rst_RASH(sort(X(ind(ii,jj),:)),150);
                u = NaN; % not sure how to handle this with RASH yet
%                 Xdist{ii,jj}(Ydist{ii,jj}==0)=[];
%                 Ydist{ii,jj}(Ydist{ii,jj}==0)=[];
        end
        Ydist{ii,jj}(1)=0;
        Ydist{ii,jj}(end)=0;
        Ydist{ii,jj}=(Ydist{ii,jj})./trapz(Xdist{ii,jj},Ydist{ii,jj})  * (-1)^(ii-1);
        
        %% Calculate jittering and spacing
        % make some space under the density plot for the boxplot and raindrops
        yl = [0 max(Ydist{ii,jj}* (-1)^(ii-1))];
        
        % width of boxplot
        wdth = yl(2) * scatter_width;
        
        % jitter for raindrops
        jit = (jitter(ind(ii,jj),:) - 0.5) * wdth;         
        % info for making boxplot
        if RSTpresent
            [quartiles(1),CIQ] = rst_hd(X(ind(ii,jj),:),0.25);
            [quartiles(2),CIQ] = rst_hd(X(ind(ii,jj),:),0.75);
            [quartiles(3),CIQ] = rst_hd(X(ind(ii,jj),:),0.5);
        else
            [quartiles(1)] = quantile(X(ind(ii,jj),:),0.25);
            [quartiles(2)] = quantile(X(ind(ii,jj),:),0.75);
            [quartiles(3)] = quantile(X(ind(ii,jj),:),0.5);
        end
        
        iqr         = quartiles(2) - quartiles(1);
        Xs          = sort(X(ind(ii,jj),:));
        
        %    whiskers(1) = min(Xs(Xs > (quartiles(1) - (1.5 * iqr))));
        %    whiskers(2) = max(Xs(Xs < (quartiles(2) + (1.5 * iqr))));
        whiskers(1) = quartiles(1) - (1.5 * iqr);
        whiskers(2) = quartiles(2) + (1.5 * iqr);
        WiskX{ii,jj}      =  whiskers;
        WiskY{ii,jj}      =  zeros(size(whiskers)) * (-1)^(ii-1);
        
        [~,tmp] = max( (ones(size(quartiles,2),1)*Xdist{ii,jj}) > (ones(size(Xdist{ii,jj},2),1) *quartiles)',[],2);
%         [~,tmp]     = max(Xdist{ii,jj}>quartiles',[],2);
        Xqrtls{ii,jj}  = [Xdist{ii,jj}(tmp);Xdist{ii,jj}(tmp)];
        Yqrtls{ii,jj}  = [zeros(1,3);Ydist{ii,jj}(tmp)];
        % raindrops
        drops_posY{ii,jj} = jit  - yl(2) * dot_dodge_amount;
        drops_posY{ii,jj} = drops_posY{ii,jj} * (-1)^(ii);
        drops_posX{ii,jj} = X(ind(ii,jj),:);
        
        % grid points 
        XGrid{ii,jj} = Xqrtls{ii,jj}(1,3);
        YGrid{ii,jj} = Yqrtls{ii,jj}(1,3);
    end
end

%% compute offsets
distW = cellfun(@(x) max(abs(x)),Ydist);
ldists = max(distW(2:2:end,2:end),[],1);
rdists = max(distW(1:2:end,1:end-1),[],1);
tmp =  rdists;
if ~isempty(ldists),tmp = tmp + ldists;end
distances = cumsum([0 tmp]).*1.05;
Yoffs = num2cell(distances);Yoffs = repmat(Yoffs,size(ind,1),1);
Xoffs = num2cell(zeros(size(Yoffs)));

%% plot all
tick = 'ytick';
nViolin =1;
for jj=1:size(ind,2) % number of violins   
    for ii=1:size(ind,1)  % cicle throungh "lobes" of each violin 
        if(isnan(ind(ii,jj))),continue;end
        if vertical
            [Xoffs{ii,jj}        , Yoffs{ii,jj}  ] = deal( Yoffs{ii,jj}     , Xoffs{ii,jj});
            [Xdist{ii,jj}     , Ydist{ii,jj}     ] = deal( Ydist{ii,jj}     , Xdist{ii,jj});
            [drops_posX{ii,jj}, drops_posY{ii,jj}] = deal( drops_posY{ii,jj}, drops_posX{ii,jj});
            [WiskX{ii,jj}     , WiskY{ii,jj}     ] = deal( WiskY{ii,jj}     , WiskX{ii,jj});
            [Xqrtls{ii,jj}    , Yqrtls{ii,jj}    ] = deal( Yqrtls{ii,jj}    , Xqrtls{ii,jj});
            tick = 'xtick';
        end
        
        %% density plot
        h{1} = fill(Xdist{ii,jj} + Xoffs{ii,jj}, Ydist{ii,jj} + Yoffs{ii,jj}, color(nViolin,:));
        hold(ax,'on');
        set(h{1}, 'EdgeColor', 'none');
        set(h{1}, 'LineWidth', line_width);
        set(h{1}, 'FaceAlpha', alpha);
        set(h{1}, 'EdgeAlpha', alpha);
        
        %% scatter
        % robust data
        h{2} = scatter(drops_posX{ii,jj}(~outliers{ii,jj}) + Xoffs{ii,jj}, drops_posY{ii,jj}(~outliers{ii,jj}) + Yoffs{ii,jj});
        h{2}.SizeData = scatter_size(ind(ii,jj),~outliers{ii,jj});
        h{2}.MarkerFaceColor = color(nViolin,:);
        h{2}.MarkerEdgeColor = 'none';
        
        % outliers
        h{3} = scatter(drops_posX{ii,jj}(outliers{ii,jj})  + Xoffs{ii,jj}, drops_posY{ii,jj}(outliers{ii,jj}) + Yoffs{ii,jj},'x');
        h{3}.SizeData = scatter_size(ind(ii,jj),outliers{ii,jj});
        h{3}.MarkerEdgeColor = color(nViolin,:);
        
        % update counter
        nViolin = nViolin + 1;
    end
end

%% Box 
nViolin =1;
for ii=1:size(ind,1)
    for jj=1:size(ind,2)
        if(isnan(ind(ii,jj))),continue;end
        if box_on            
            
            x1 = Xqrtls{ii,jj}(1,1) + Xoffs{ii,jj};
            x2 = Xqrtls{ii,jj}(1,2) + Xoffs{ii,jj};
            xm = Xqrtls{ii,jj}(1,3) + Xoffs{ii,jj};
            ym = Yqrtls{ii,jj}(1,3) + Yoffs{ii,jj};
            hh = (xlim)*0.1; hh=hh(2);
            [x,y] = roundedRect(mean([x1 x2]),ym,mean([x1 x2]),hh,hh./2);
            MpointOffset = hh/3*1.5;
            
            h{4} = line(whiskers,...
                WiskY{ii,jj} + Yoffs{ii,jj}+ 0.8*MpointOffset,...
                'col', [73,73,73]./255, 'LineWidth', line_width);
            
            
            h{5} = patch(x,y,[73,73,73]./255,'EdgeColor','none');
            h{6} = scatter(xm,ym+MpointOffset,MpointOffset*400,[1 1 1],'filled');
            
           
            
            % 'box' of 'boxplot'
            %          h{3} = rectangle('Position', box_pos,'Curvature',0.2);
            %          set(h{3}, 'EdgeColor', bxcl)
            %          set(h{3}, 'LineWidth', line_width);
            %set(h{3}, 'FaceColor', bxfacecl);
            % could also set 'FaceColor' here, etc
            nViolin = nViolin + 1;
        end
    end
end

%% Grid
if grid_on
    yl = get(gca,'ylim');
    xl = get(gca,'xlim');
    nViolin =1;
    if vertical
        [xl, yl] = deal( yl, xl);
    end
    for jj=1:size(ind,2) 
        for ii=1:size(ind,1)
            if(isnan(ind(ii,jj))),continue;end
            YGrid{ii,jj} = yl(1);
            Marker = '^';
            if vertical
                Marker = '>';
                [XGrid{ii,jj}    , YGrid{ii,jj}    ] = deal( YGrid{ii,jj}    , XGrid{ii,jj});
            end
            %% grid
            h{7} = scatter(XGrid{ii,jj},YGrid{ii,jj},Marker);
            h{7}.MarkerFaceColor = color(nViolin,:);
            h{7}.MarkerEdgeColor = 'none';
%             uistack(h{7},'bottom')       
            nViolin = nViolin + 1;
        end
    end
    
    % restore axes limits
    if vertical
        [xl, yl] = deal( yl, xl);
    end
    set(gca,'xlim',xl);
    set(gca,'ylim',yl);
end

%% Ticks
if vertical
    ax.XTick = unique([Xoffs{:}]);
    ax.XTickLabel = arrayfun(@(x) {num2str(x)},1:size(ind,2));
else
    ax.YTick = unique([Yoffs{:}]);
    ax.YTickLabel = arrayfun(@(x) {num2str(x)},1:size(ind,2));
end
hold(ax,'off');
box(ax,'off');
end

%% Colors for my thesis
function Clrs=MyColors(ind)
%%%% RGB

Clrs=[
      191 24  24 ;
      55  123 191;
      106 196 100;
      232 143 34 ;
      
      224 82  82 ;
      106 182 193;
      232 171 27 ;
      
      122 122 122;
      179 179 179;
      55  69  72 ;
      
      126 47  142;
      119 172 48 ;
      77  190 238;
      ]./255;

%   Clrs = cubehelix(16,2.8,1.9,1.45,0.6,[0.05 0.95],[0 .7]);
   
  if nargin<1
      ind=1:size(Clrs,1);
      f=figure('Color',[1 1 1],'Units','normalized');
      f.Position=[.01 .15 .07*size(Clrs,1) .7];
      ax=axes(f, 'Position',[0.01 0.2 0.98 0.7]);
      
      hBar=bar(ax,ind,ones(size(ind)),'EdgeColor','flat','FaceColor','flat','CData',Clrs);
      set(ax.YAxis,'Visible','off');
      set(ax,'XLim',[ind(1)-.5 ind(end)+.5],'Box','off');
      ax.XTickLabel=cellstr(num2str(Clrs));
      
      if exist('AnnotationToDataPoint')==2
          for ii=ind
              ann=AnnotationToDataPoint(ax,'textbox',ii,.8,'String',num2str(ii));
              ann.FontSize=20;
          end
      else
          fprintf('Consider downloading the function AnnotationToDataPoint!')
      end
  elseif strcmp(ind,'all')
          ind=1:size(Clrs,1);
      end
  
Clrs=Clrs(mod(ind-1,size(Clrs,1))+1,:);
end

function [x,y] = roundedRect(x0,y0,a,b,r)
% n = 1000;
% k = 1:n;
% x1 = x0-a;
% x2 = x0+a;
% x = (x1+x2)/2+(x2-x1)/2.*cos((2.*k-1)./(2*n).*pi);
x = linspace(-a,a,1000);
expb = 2*b/r;
expa = 2*a/r;
X = abs( x )./a;
y = [abs( b.*(1 - X.^expa ).^(1/expb) )];
x = x + x0;
y = y + y0;
end