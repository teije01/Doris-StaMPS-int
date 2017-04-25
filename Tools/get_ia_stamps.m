function get_ia_stamps(insarpath,len,width)
% function that computes the look angle for the stamps results based on the
% look_angle.1.in file and the wdith and length of the multi-looked
% (cropped) SAR image.
%
%     Copyright (C) 2015  Bekaert David - University of Leeds
%     Email: eedpsb@leeds.ac.uk or davidbekaert.com
% 
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 2 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License along
%     with this program; if not, write to the Free Software Foundation, Inc.,
%     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
%
% By Bekaert David - Unviersity of Leeds
% modified from the StaMPS script ps_load_initial
% 
% get_la_stamps(insarpath,len,width)

% the file names
if nargin == 0
    fprintf('Loading variables from current folder\n')
    insarpath = '.';
    fprintf('Insarpath = ''.''\n')
    load('width.txt')
    load('len.txt')
    fprintf('width     = %d\n',width)
    fprintf('len       = %d\n',len)
elseif nargin ~= 3
    error('wrong number of input variables')
end

load([insarpath filesep 'psver.mat']);
ianame = [insarpath filesep 'incidence_angle.1.in'];
ps = load([insarpath filesep  'ps' num2str(psver) '.mat']);
hgt = load([insarpath filesep 'hgt' num2str(psver) '.mat']);
currdir = pwd;

% getting the variables
ij = ps.ij;
hgt= hgt.hgt;


% Get matlab version as function arguments change with the matlab version
matlab_version = version('-release');           % [DB] getting the matlab version
matlab_version = str2num(matlab_version(1:4));  % [DB] the year

% constructing the grid
[gridX,gridY]=meshgrid(linspace(0,width,50),linspace(0,len,50));

% computing the incidence angles
incid_angle=load(ianame);
ia0=incid_angle(1:2:end)*pi/180;
ia0=reshape(ia0,50,50)';
ia1000=incid_angle(2:2:end)*pi/180;
ia1000=reshape(ia1000,50,50)';

ia0_ps=griddata_version_control(gridX,gridY,ia0,ij(:,3),ij(:,2),'linear',matlab_version);           % [DB] fix matlab2012 version and older
ia1000_ps=griddata_version_control(gridX,gridY,ia1000,ij(:,3),ij(:,2),'linear',matlab_version);     % [DB] fix matlab2012 version and older
ia=ia0_ps+(ia1000_ps-ia0_ps).*hgt/1000;
lasavename=[insarpath filesep 'ia',num2str(psver)];
save(lasavename,'ia');

cd(insarpath)
ps_plot(ia)

cd(currdir)
end


