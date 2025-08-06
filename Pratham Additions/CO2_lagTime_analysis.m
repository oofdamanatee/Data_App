%% FULL DATA ANALYSIS of CO2 Lag Time - from start to finish - 

% ---- Variables that require user manual input will be surrounded ----
    % like this
% ----

%%
% ---- IMPORTANT!!! Input the correct data right from the start. This will put
% all the data in the right place. If the date is wrong, it could overwrite
% previous analysis.
date_of_experiment = "2025-08-05";
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
spectra_range = 1:46; 
% ----

% --- the spectra file prefix ---
file_prefix = 'PA_20250805_Trial7N2purge_';
% ----

cd(data_path)
[data1,freq] = LoadSpectra(data_path,file_prefix,spectra_range);
freq = freq(:,1);

if freq(2) - freq(1) > 0
    freq = flip(freq);
end

fprintf("Successfully imported " + size(data1, 2) + " spectra.\n")

%% Data Analysis

% --- experimental parameters ---
spacer_size = 12;  % in microns
flow_rate = 4.0; % in mL/min
timeZero = datetime('05-Aug-2025 11:03:39'); % According to 24hr clock as hr:min:sec
your_name = "Pratham";
% ---

%Get time stamps of every spectrum after time zero
times = timeAxis(data_path,file_prefix,spectra_range, timeZero);

% Subtract the initial spectrum (considering spectra after time zero)
ignored_spectra = numel(spectra_range) - numel(times); % Number of spectra taken before time zero (to be ignored)
sub_data = data1(:, (2+ignored_spectra):end) - data1(:,(2+ignored_spectra)); % Contains difference spectra

%% flipbook of spectra

n_spectra = size(sub_data, 2);
for ii= 1:n_spectra
    figure(10)
    plot(freq, sub_data(:,ii))
    set(gca, "YLim", [-0.12, 0.2])  
    pause  

end

%%
% Find the integrated absorbance of CO2 for each used spectrum
abs_CO2 = IntegratedAbsorbance(sub_data, freq);

% Total time taken to reach equilibrium CO2 pressure in brass cell
[useless, ind_Max] = max(abs_CO2);
lagTime = times(ind_Max) / 60; % in minutes
%% Plot the Integrated Absorbance

if flow_rate > 2
    marker = 15;
else
    marker = 6;
end
figure;
plot(times(2:end), abs_CO2, "r.", MarkerSize= marker)
title('Integrated absorbance of CO2 in Empty Brass Cell Over Time')
xlabel('Time (s)')
ylabel('Integrated Absorbance')

annotation('textbox', [.6,.4, .3, .3], 'String',{'  Date: ' + date_of_experiment, '  CO2 Flow rate (mL/min): '+ string(flow_rate), '  Spacer size (µm): ' + string(spacer_size)}, 'FitBoxToText', 'on')


%% Export figure
ax = gca;
% ------ Figure Name ------
exportgraphics(ax, '/Users/oofdamanatee/Downloads/Figure_8.pdf')
% ------

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

% Find integrated CO2 absorbance for spectra
function absorbance = IntegratedAbsorbance(spectra, freq_list, varargin)
    range = [2295, 2390];

    if nargin == 3
       range = varargin{1};
    end

    % Indices of CO2 gas line frequencies
    CO2Freq = (range(1) < freq_list) & (freq_list < range(2));



    % Integrated absorbance of CO2 for each spectrum
    integrated_abs_CO2 = [];
    for k = 1:size(spectra, 2)
        spectrum = spectra(CO2Freq, k);
        integrated_abs_CO2 = [integrated_abs_CO2 trapz(spectrum)];
    end

    absorbance = integrated_abs_CO2;
end
