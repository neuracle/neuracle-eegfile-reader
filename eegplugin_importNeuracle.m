% eegplugin_importNeuracle() - EEGLAB plugin for importing Neuracle data files.
%
% Usage:
%   >> eegplugin_importNeuracle(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Notes:
%   This plugins consist of the following Matlab files:
%
% Create a plugin:
%   For more information on how to create an EEGLAB plugin see the
%   help message of eegplugin_besa() or visit http://www.sccn.ucsd.edu/eeglab/contrib.html
%
% Author: Xiaoshan HUANG, hxs@neuracle.cn
%
% Versions:
%    v1.0: 2017-09-27, orignal
%    v1.1: 2019-12-26, update pop_importNeuracle.m

% Copyright (c) 2020 Neuracle, Inc. All Rights Reserved. http://neuracle.cn
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
function eegplugin_importNeuracle(fig, trystrs, catchstrs)
    % add folder to path
    % ------------------
    if exist('pop_importNeuracle', 'file')
        p = which('eegplugin_importNeuracle.m');
        p = p(1:findstr(p,'eegplugin_importNeuracle.m')-1);
        addpath(p);
    end;
    % find import data menu
    % ---------------------
    menu = findobj(fig, 'tag', 'import data');
    
    % menu callbacks
    % --------------
    comcnt = [ trystrs.no_check '[EEG, command] = pop_importNeuracle;' catchstrs.new_and_hist ];
    
    % create menus
    % ------------
    uimenu( menu, 'label', 'From Neuracle EEG data files', 'callback', comcnt);

end

