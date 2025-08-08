function [y, x, varargout] = binAveraging(y0, x0, varargin)
% binAveraging computes non-overlapping bin-averaged quantities from vectors.
%
% Syntax:
%   [y, x] = binAveraging(y0, x0)
%   [y, x, dispersion] = binAveraging(y0, x0, 'Name', Value, ...)
%
% Inputs:
%   y0: vector [N x 1] of ordinate values (e.g., a time series)
%   x0: vector [N x 1] of abscissa values (e.g., a time vector)
%
% Name-Value Pair Arguments:
%   'newX'      : vector specifying new bin edges [M x 1]
%   'binMethod' : method for automatic binning (default: 'auto')
%   'BMIN'      : scalar specifying minimum bin limit (default: min(x0))
%   'BMAX'      : scalar specifying maximum bin limit (default: max(x0))
%   'binWidth'  : scalar specifying bin width
%   'Nbin'      : scalar specifying number of bins
%   'averaging' : 'mean' or 'median' (default: 'mean')
%   'dispersion': 'std', 'decile_10_90', 'decile_25_75', or 'decile_1_99' (default: 'std')
%
% Outputs:
%   y        : vector [M x 1] of binned values of y0
%   x        : vector [M x 1] of corresponding abscissa values
%   varargout: optional output containing dispersion measure [M x 1] or [M x 2]
%
% Example:
%   [y_avg, x_bins, y_std] = binAveraging(y_data, x_data, 'binWidth', 0.5, 'averaging', 'mean', 'dispersion', 'std');
%
% Author: E. Cheynet  - UiB - Norway
% Last modified: [Your Modification Date]
%
% See also: histcounts, accumarray

%% Input Validation
narginchk(2, inf);
nargoutchk(0, 3);

% Ensure that x0 and y0 are column vectors
x0 = x0(:);
y0 = y0(:);

% Check that x0 and y0 are vectors of the same length
if length(x0) ~= length(y0)
    error('x0 and y0 must be vectors of the same length.');
end

% Input parser
p = inputParser();
p.CaseSensitive = false;
p.addParameter('newX', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('binMethod', 'auto', @(x) ischar(x) || isstring(x));
p.addParameter('BMIN', nanmin(x0), @isnumeric);
p.addParameter('BMAX', nanmax(x0), @isnumeric);
p.addParameter('binWidth', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('Nbin', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('averaging', 'mean', @(x) any(strcmpi(x, {'mean', 'median'})));
p.addParameter('dispersion', 'std', @(x) any(strcmpi(x, {'std', 'decile_10_90', 'decile_25_75', 'decile_1_99'})));
p.parse(varargin{:});

% Assign variables from parsed inputs
newX = p.Results.newX;
binMethod = p.Results.binMethod;
BMIN = p.Results.BMIN;
BMAX = p.Results.BMAX;
binWidth = p.Results.binWidth;
Nbin = p.Results.Nbin;
averaging = lower(p.Results.averaging);
dispersion = lower(p.Results.dispersion);

%% Determine Bin Indices
if ~isempty(newX)
    % Use specified bin edges
    [~, ~, bin] = histcounts(x0, newX);
elseif ~isempty(Nbin)
    % Use specified number of bins
    [~, ~, bin] = histcounts(x0, Nbin);
elseif ~isempty(binWidth)
    % Use specified bin width
    [~, ~, bin] = histcounts(x0, 'BinLimits', [BMIN, BMAX], 'BinWidth', binWidth);
else
    % Use automatic binning method
    [~, ~, bin] = histcounts(x0, 'BinMethod', binMethod, 'BinLimits', [BMIN, BMAX]);
end

% Remove zero bins (values outside bin limits)
validBins = bin > 0;
x0 = x0(validBins);
y0 = y0(validBins);
bin = bin(validBins);

%% Compute Averaged Values
switch averaging
    case 'mean'
        y = accumarray(bin, y0, [], @nanmean, NaN);
        x = accumarray(bin, x0, [], @nanmean, NaN);
    case 'median'
        y = accumarray(bin, y0, [], @nanmedian, NaN);
        x = accumarray(bin, x0, [], @nanmedian, NaN);
    otherwise
        error('Invalid value for ''averaging''. Must be ''mean'' or ''median''.');
end

%% Compute Dispersion
switch dispersion
    case 'std'
        stdY = accumarray(bin, y0, [], @nanstd, NaN);
    case 'decile_10_90'
        lowerDecile = accumarray(bin, y0, [], @(x) quantile(x, 0.10), NaN);
        upperDecile = accumarray(bin, y0, [], @(x) quantile(x, 0.90), NaN);
        stdY = [lowerDecile, upperDecile];
    case 'decile_25_75'
        lowerQuartile = accumarray(bin, y0, [], @(x) quantile(x, 0.25), NaN);
        upperQuartile = accumarray(bin, y0, [], @(x) quantile(x, 0.75), NaN);
        stdY = [lowerQuartile, upperQuartile];
    case 'decile_1_99'
        lowerPercentile = accumarray(bin, y0, [], @(x) quantile(x, 0.01), NaN);
        upperPercentile = accumarray(bin, y0, [], @(x) quantile(x, 0.99), NaN);
        stdY = [lowerPercentile, upperPercentile];
    otherwise
        error('Invalid value for ''dispersion''. Must be one of ''std'', ''decile_10_90'', ''decile_25_75'', or ''decile_1_99''.');
end

%% Remove NaNs
nanMask = isnan(x);
x(nanMask) = [];
y(nanMask) = [];
if exist('stdY', 'var')
    stdY(nanMask, :) = [];
end

%% Output Assignment
if nargout >= 3
    varargout{1} = stdY;
end

end



%%
% function [y,x,varargout] = binAveraging(y0,x0,varargin)
% % function [y,x,varargout] = binAveraging(Y,x0,varargin) computes
% % non-overlapping bin-averaged quantities from vectors.
% % 
% % Input: 
% % x0: vector [1 x N] of absiccsa values, e.g. a time vector
% % y0: vector [1 x N] of ordinate values in each x0, e.g. a time series
% % varargin
% % 
% % 
% % options for Input:
% % newX: vector [1 x M]: new abcissa
% % Nbin: scalar [1 x 1]: number of bins (overwritten by newX, if specified) 
% % binWidth: scalar [1 x 1]: width of bins (overwritten by newX, if specified) 
% % BMAX: scalar: Max value of the bin (x0<=BMAX) 
% % BMIN: scalar: Min value of the bin (x0>=BMIN) 
% % averaging: string: 'mean" or 'median'
% % dispersion: string: 'std', 'decile_10_90' or 'decile_25_75'
% % 
% % Output:
% % y: vector [1 x M]: binned values of y0
% % x: vector [1 x M]: new absiccsa
% % varargout: if specified, provides a [ 1x M] vector for the dispersion
% % 
% % Author: E. Cheynet  - UiB - Norway
% % LAst modified 06/12/2019
% % 
% %  see also histcounts accumarray
% 
% %%  Inputparser
% 
% p = inputParser();
% p.CaseSensitive = false;
% p.addOptional('newX',[]);
% p.addOptional('binMethod','auto');
% p.addOptional('BMIN',nanmin(x0));
% p.addOptional('BMAX',nanmax(x0));
% p.addOptional('binWidth',[]);
% p.addOptional('Nbin',[]);
% p.addOptional('averaging','mean');
% p.addOptional('dispersion','std');
% p.parse(varargin{:});
% % shorthen the variables name
% BMIN = p.Results.BMIN ;
% BMAX = p.Results.BMAX ;
% binMethod = p.Results.binMethod ;
% binWidth = p.Results.binWidth ;
% Nbin = p.Results.Nbin;
% newX = p.Results.newX ;
% averaging = p.Results.averaging ;
% dispersion = p.Results.dispersion ;
% y0 = y0(:);
% x0 = x0(:);
% 
% %%
% if ~isempty(newX)
%     [~,~,bin] = histcounts(x0,newX);
% elseif ~isempty(Nbin)
%     [~,~,bin] = histcounts(x0,Nbin);
% elseif isempty(newX) && isempty(binWidth)
%     [~,~,bin] = histcounts(x0,'BinMethod',binMethod,'BinLimits',[BMIN,BMAX]);
% elseif isempty(newX) && ~isempty(binWidth)
%     [~,~,bin] = histcounts(x0,'BinMethod',binMethod,'BinLimits',[BMIN,BMAX],'binWidth',binWidth);
% else
%     error('unspecfied options')
% end
% % bin(bin==0)=1;
% 
% 
% %% Mean values
% if strcmpi(averaging,'mean')
%     y = accumarray(bin(bin>0), y0(bin>0), [], @nanmean, nan);
%     x = accumarray(bin(bin>0), x0(bin>0), [], @nanmean, nan);
% elseif strcmpi(averaging,'median')
%     y = accumarray(bin(bin>0), y0(bin>0), [], @nanmedian, nan);
%     x = accumarray(bin(bin>0), x0(bin>0), [], @nanmedian, nan);
% else
%     error(' ''averaging'' should be ''mean'' or ''median'' ')
% end
% 
% 
% %% Dispersion
% clear dY
% if strcmpi(dispersion,'std'),
%     dY = accumarray(bin(bin>0), y0(bin>0), [], @nanstd, nan);
% elseif strcmpi(dispersion,'decile_10_90'),
%     dY(1,:) = accumarray(bin(bin>0), y0(bin>0), [], @(x) quantile(x,0.1), nan);
%     dY(2,:) = accumarray(bin(bin>0), y0(bin>0), [], @(x) quantile(x,0.9), nan);
% elseif strcmpi(dispersion,'decile_25_75'),
%     dY(1,:) = accumarray(bin(bin>0), y0(bin>0), [], @(x) quantile(x,0.25), nan);
%     dY(2,:) = accumarray(bin(bin>0), y0(bin>0), [], @(x) quantile(x,0.75), nan);    
% else
%     error(' ''std'' should be ''decile_10_90'' or ''decile_25_75'' ')
% end
% 
% dY = dY';
% if size(dY,1)==2
%     dY(isnan(x),:)=[];
%     y(isnan(x))=[];
% x(isnan(x))=[];
% elseif size(dY,1)==1
%     dY(isnan(x))=[];
%     y(isnan(x))=[];
% x(isnan(x))=[];
% else
%     warning('The size should be 2 or 1')
% end
% 
% 
% 
% 
% if nargout ==3,    varargout = {dY};end
% 
% 
% 
% end