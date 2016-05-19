function [files, subdirs, names] = get_matFiles(basepath, pattern)
%GET_MATFILES search for math files in a given directory with a given
%regexp pattern return full paths, folders and the names of the files
%
%          search for all files recursively, inside the basepath
%          then filter for pattern files contain the complete path while
%          subdirs only the relative path to basepath and names are the
%          filenames
%
%
%Marcell Stippinger, 2016

fprintf('Searching files in %s with pattern %s\n', basepath, pattern);
[f, s, n] = rvisit(basepath,'');

match = regexpi(f,pattern);
% match needs to be a matrix for indexing, and boolean indexing does not
% work...
match = find(~cellfun(@isempty,match));
% {} would give the data contained, () results cells
files = f(match);
subdirs = s(match);
names = n(match);

for ifile = 1 : length(files)
    fprintf('%3d: %s\n', ifile, files{ifile})
end



function [files, subdirs, names] = rvisit(basepath, subpath)
if (~isempty(basepath)) && (basepath(end) ~= '/')
    basepath = [basepath '/'];
end
if (~isempty(subpath)) && (subpath(end) ~= '/')
    subpath = [subpath '/'];
end
fullpath = [basepath subpath];
dinfo = dir(fullpath);
filecount = sum(~[dinfo.isdir]);
if filecount
    names = {dinfo(~[dinfo.isdir]).name};
    subdirs(1:filecount) = {subpath};
    files = strcat(fullpath,names);
else
    names = {};
    subdirs = {};
    files = {};
end

for idir = find([dinfo.isdir])
    if ~(strcmp(dinfo(idir).name, '..') || strcmp(dinfo(idir).name, '.'))
        [f, s, n] = rvisit(basepath,fullfile(subpath, dinfo(idir).name));
        l = length(f);
        files(end+1:end+l) = f;
        subdirs(end+1:end+l) = s;
        names(end+1:end+l) = n;
    end
end
        


function [files, names, roots] = get_matFiles_old(varargin)
%GET_MATFILES search for math files in a given directory with a given pattern
%           returns the names of the files, folders, and full paths
%
%          [files, names, roots] = get_matFiles('/','header_','.mat') search for all mat files
%          recursively, inside the folder '/' that are contained in subfolders with prefix
%          "header_"
%
%
%Ruben Pinzon@2015

basepath = varargin{1};
if nargin == 1
   pattern  = '/*BehavElectrData.mat';
   sufix    = 'i01*';
elseif nargin == 2
   sufix    = 'i01*';
   pattern  = varargin{2};
else
   sufix    = varargin{3};
   pattern  = varargin{2};
end

fprintf('Searching files with pattern %s\n', [basepath sufix pattern])
animals = dir([basepath sufix]);
files = {};
names = {};
roots = {};
for ifolder = 1 : length(animals)
    name = animals(ifolder).name;
    D = dir([basepath name pattern]);
    for ifil = 1 : length(D)
        files{end + 1} = [basepath name '/' D(ifil).name];
        names{end + 1} = name;
        roots{end + 1} = [basepath name '/'];
        disp(sprintf('%d: %s', length(files), files{end}));
    end
end


function [files, names, roots] = get_datFiles(varargin)
%GET_MATFILES search for math files in a given directory with a given pattern
%           returns the names of the files, folders, and full paths
%
%          [files, names, roots] = get_matFiles('/','header_','.mat') search for all mat files
%          recursively, inside the folder '/' that are contained in subfolders with prefix
%          "header_"
%
%
%Stippinger Marcell

basepath = varargin{1};
if nargin == 1
   pattern  = '/*BehavElectrData.mat';
   sufix    = 'i01*';
elseif nargin == 2
   sufix    = 'i01*';
   pattern  = varargin{2};
else
   sufix    = varargin{3};
   pattern  = varargin{2};
end

fprintf('Searching files with pattern %s\n', [basepath sufix pattern])

D = dir([basepath sufix pattern]);
files = cell(1,length(D));
names = cell(1,length(D));
roots = cell(1,length(D));
for ifil = 1 : length(D)
    files{ifil} = [basepath sufix '/' D(ifil).name];
    names{ifil} = strrep(D(ifil).name,'.dat','');
    roots{ifil} = [basepath sufix '/' ];
    fprintf('%d: %s\n', ifil, files{ifil});
end
