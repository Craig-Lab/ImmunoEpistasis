function installMyDdesd()

myfilename = 'ddesd_f5.m';

% Where are we
thisfile = mfilename('fullpath');
thispath = fileparts(thisfile);

% Construct source file
myddesdpath = fullfile(thispath,myfilename);

% Where are we going?
installfolder = fullfile(matlabroot,'toolbox','matlab','funfun',myfilename);

% Copy
copyfile(myddesdpath,installfolder,'f')

rehash toolbox
rehash toolboxcache

end