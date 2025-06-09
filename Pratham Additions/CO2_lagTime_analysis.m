%% FULL DATA ANALYSIS of CO2 Lag Time - from start to finish - 

% ---- Variables that require user manual input will be surrounded ----
    % like this
% ----

%%
% ---- IMPORTANT!!! Input the correct data right from the start. This will put
% all the data in the right place. If the date is wrong, it could overwrite
% previous analysis.
date_of_experiment = "2025-06-05";
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
spectra_range = 1:68; 
% ----

% --- the spectra file prefix ---
file_prefix = 'PA_20250605_Trial1_';
% ----

% --- experimental parameters ---
spacer_size = 12;  % in microns
flow_rate = 4; % in mL/min
timeZero = datetime('05-Jun-2025 14:21:01'); % According to 24hr clock as hr:min:sec
sample_name = "4.00 mL/min flow rate";
your_name = "Pratham";
% ---

cd(data_path)
[data1,freq] = LoadSpectra(data_path,file_prefix,spectra_range);
freq = freq(:,1);

if freq(2) - freq(1) > 0
    freq = flip(freq);
end

fprintf("Successfully imported " + size(data1, 2) + " spectra.\n")
%% Data Analysis

%Get time stamps of every spectrum after time zero
times = timeAxis(data_path,file_prefix,spectra_range, timeZero);

% Subtract the initial spectrum (considering spectra after time zero)
ignored_spectra = numel(spectra_range) - numel(times); % Number of spectra taken before time zero (to be ignored)
sub_data = data1(:, (1+ignored_spectra):end) - data1(:,(1+ignored_spectra)); % Contains difference spectra

% Find the maximum absorbance of CO2 for each used spectrum
abs_CO2 = maxAbsorbance(sub_data, freq);

% Total time taken to reach equilibrium CO2 pressure in brass cell
lagTime = times(end); % in minutes

figure;
plot(times, abs_CO2, "r.")
title('Maximum Absorbance of CO2 in Empty Brass Cell Over Time')
xlabel('Time (s)')
ylabel('Relative Absorbance')

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

% Find max CO2 absorbance for spectra
function absorbance = maxAbsorbance(spectra, freq_list)
    % Indices of CO2 gas line frequencies
    CO2Freq = (2200 < freq_list) & (freq_list < 2400);

    % Max absorbance of CO2 for each spectrum
    max_abs_CO2 = [];
    for k = 1:size(spectra, 2)
        spectrum = spectra(:, k);
        max_abs_CO2 = [max_abs_CO2 max(spectrum(CO2Freq))];
    end

    absorbance = max_abs_CO2;
end
