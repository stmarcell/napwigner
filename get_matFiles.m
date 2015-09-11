function [files, names, roots] = get_matFiles(basepath)

animals = dir([basepath 'i01*']);
for ifolder = 1 : length(animals)
    names{ifolder} = animals(ifolder).name;
    D = dir([basepath names{ifolder} '/*BehavElectrData.mat']);
    files{ifolder} = [basepath names{ifolder} '/' D.name];
    roots{ifolder} = [basepath names{ifolder} '/' animals(ifolder).name];
end 