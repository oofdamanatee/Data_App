%% Analyzing MOF Spectra 

% ---- Variables that require user manual input will be surrounded ----
    % like this
% ----

%%
% ---- IMPORTANT!!! Input the correct data right from the start. This will put
% all the data in the right place. If the date is wrong, it could overwrite
% previous analysis.
date_of_experiment = "2025-07-25";
% ----
year_of_experiment = year(datetime(date_of_experiment));

% ---- here put your path to CHEM-SGR (end with a slash). This will be necessary
% throughout the script.
isilon_path = "/Volumes/CHEM-SGR/";
% ----
data_path = isilon_path + "sgr-ftir.chem.pitt.edu/" + year_of_experiment + "/" + date_of_experiment;

try
    cd(isilon_path)
catch
    warning("An error ocurred accessing the data directory. Make sure you are connected to Isilon.")
end

%% Step 1:Load in the spectra
cd ~
% --- the indices of the spectra you wish to use ----
spectra_range = 1; 
% ----

% --- the spectra file prefix or full name (if only using 1 spectrum) ---
file_name = 'PA_20250725_Methanol_';
% ----

cd(data_path)
[data1,freq] = LoadSpectra();
freq = freq(:,1);

if freq(2) - freq(1) > 0
    freq = flip(freq);
end

fprintf("Successfully imported " + size(data1, 2) + " spectra.\n")
%% Plot the Spectra
plot(freq, data1)


%% Functions

% Find the time of each spectrum in seconds after the first spectrum taken
% after time zero
function timeArray = timeAxis(varargin)
    if nargin == 1
        filenames = uigetfile({'*.spa','Thermo Spectrum (*.spa)'}, ...
            'MultiSelect','on','Select Spectra Files...');
    elseif nargin == 4
        pathname = varargin{1};
        fileroot = varargin{2};
        nums = varargin{3};
        timeZero = varargin{4};
        
        cd(pathname);
        
        files = dir(fileroot + "*");
        filenames = {files(nums).name};
    else
        error("Error: Invalid set of arguments. ...If using arguments, enter (filepath,fileroot,nums)")
    end
    times = [];
    for ii = 1:numel(filenames)
        g = dir(filenames{ii});
        times = [times datetime(g.date)]; %appends new datetime to times
    end
    sec = [];
    for ii = 1:numel(times)
        sec = [sec seconds(times(ii) - timeZero)]; %Appends datetime to sec as seconds after first spectrum
    end
    timeArray = sec(sec>0);
    %             axis = (0:(size(obj.data,2)-1)).*obj.timeInterval;
end



