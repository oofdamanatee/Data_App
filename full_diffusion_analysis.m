%% FULL DATA ANALYSIS of diffusion experiment - from start to finish - 

% ---- Variables that require user manual input will be surrounded ----
    % like this
% ----

%%
% ---- IMPORTANT!!! Input the correct data right from the start. This will put
% all the data in the right place. If the date is wrong, it could overwrite
% previous analysis.
date_of_experiment = "2025-05-12";
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
%% Step 1: process the microscope image

%% Load the image from Isilon

% ---- path to the image within CHEM-SGR (end with a slash) ----
image_path = "sgr-kiralux-camera.chem.pitt.edu/2025-05-12/";
% ----

% ---- name of the image file (use the .tif one) ----
image_filename = "75EMIMNTF2inPEGDA_20250512_pre-diffusion_roomtemp.tif";
% ----
%% Process the image

% Get the radius from the image
% ------------
[radius_pixels,displacement_pixels,other] = radius_from_image(isilon_path+image_path+image_filename,"image scale",25,"dx definition","from center","radius definition","centroid");
% ------------

% Convert the radius to distance
% ------------

% ------------
%% Get the scaling

units = "um";
path_to_ref = isilon_path + "sgr-kiralux-camera.chem.pitt.edu/2025-05-09/";
um_per_pixel = get_distance_per_pixel(path_to_ref,"scaled_reticle01.tif",[1476 1939;824 828],2000,units,"image scale",10);

%% Calculate the radius
radius_mm = radius_pixels*um_per_pixel/1000;
displacement_mm = displacement_pixels*um_per_pixel/1000;
fprintf("r = %4f mm; dx = %4f mm\n",radius_mm,displacement_mm)
%% Save the data
% make a neat structure
image_processing.path = image_path + image_filename;
image_processing.radius_mm = radius_mm;
image_processing.displacement_mm = displacement_mm;
image_processing.length_per_pixel = um_per_pixel;
image_processing.units = units;
image_processing.radius_definition = "centroid";
image_processing.other_data = other;

cd(data_path)
save(date_of_experiment + "_image_processing_data.mat",'image_processing')
%% ---- END OF STEP 1 ----

%% Step 2: Fit the FTIR peaks to obtain the uptake curve

%% Load in the spectra
cd ~
% --- the indicies of the spectra you wish to use ----
spectra_range = [1:115]; 
% ----

% --- the spectra file prefix ---
file_prefix = '75EMIMNTF2inPEGDA_20250512_room_';
% ----

% --- experimental parameters ---
volume = 0.07;  % in microliters
spacer_size = 12;  % in microns
gel_radius = radius_mm*1000;  % as previously calculated, but can be overridden
time_delay = 60;  % between spectra, in seconds
sample_name = "75% EMIM NTF2 in PEGDA";
your_name = "Matt";
% ---

cd(data_path)
[data1,freq] = LoadSpectra(data_path,file_prefix,spectra_range);
freq = freq(:,1);

if freq(2) - freq(1) > 0
    freq = flip(freq);
end

% Subtract to the initial spectrum
sub_data = data1 - data1(:,1);

% INITIALIZE OBJECT
f = FTIRexperiment(sub_data,freq,volume,spacer_size,gel_radius,...
    time_delay,sample_name,date_of_experiment,your_name);
f = f.timeAxis(data_path,file_prefix,spectra_range);

fprintf("Successfully imported " + size(f.data,2) + " spectra.\n")
%% Guesses for FTIR peak fitting, by eye
% ---- Which spectrum will you match to? Usually the last one is good.
trial_spectrum = 110;
% ----

% set the fit range. Usually doesn't need to be changed
range1 = [2290 2390];

% ---- User-input starting point values ----
sp.center = 2341.25;
sp.wg = 1.7; 
sp.wl = 1.7;
sp.a1 = 2.25;  % main peak height
sp.a2 = 0.07; % expected Boltzmann factor for bend
sp.a3 = 0.0; % gas lines
sp.c0 = 0.02;
sp.c1 = 0; % baseline slope
% ----

%fit function requires fliipped inputs
freq = flip(f.freqAxis);
s = flip(f.data(:,trial_spectrum));


%get x and y data for the fit
ind1 = find(freq>=range1(1) & freq<range1(2));
x = freq(ind1);
ydata = s(ind1);

%plot the fitted function using user parameters
yfit = co2GasLineFitFunction(x,sp.center,sp.wg,sp.wl,sp.a1,sp.a2,sp.a3,sp.c0,sp.c1);
res = ydata-yfit;
sse = sum(res.^2);

figure(1);clf
plot(x,ydata,'o',x,yfit,x,res-0.1,'r-o')
%% Do the FTIR peak fit
T = tic; %time the fitting for later display
f = gasLineFit(f,sp.center,sp.wg,sp.wl,sp.a1,sp.a2,sp.a3,sp.c0,sp.c1);
stop = toc(T);

%selecte 4 evenly placed fits to plot
n_spectra = size(f.data,2);
iis = ceil([1 n_spectra/4 n_spectra*3/4 n_spectra]);
figure(2);clf
for ii = iis
    plot(f.fittedSpectra(ii).x,f.fittedSpectra(ii).ydata,'o',...
        f.fittedSpectra(ii).x,f.fittedSpectra(ii).yfit,...
        f.fittedSpectra(ii).x,f.fittedSpectra(ii).res-0.1,'ro')
    hold on
end
hold off

%let the user know how it went
review = "";
tl = 0;
for ii = 1:n_spectra
    if f.fittedSpectra(ii).O.exitflag < 1
        review = [review;'Spectrum '+ii+' did not converge!!! Results might not be trustworthy.'];
        tl = tl+1;
    end
end
if tl==0
    review = "All fits were successful.";
end
review = [review;"Fitting took "+stop+" seconds."];
review
 %% Plotting the uptake curve for viewing
figure(3);clf

% number of spectra to show
n = size(f.data,2);

%find the indicies for the amount of spectra desired
spectraIndicies = zeros(1,n);
interval = ceil(size(f.data,2)/n);
for ii = 1:n
    spectraIndicies(ii) = (ii*interval);
end

for ii = spectraIndicies
    temp = f.fittedSpectra(ii).fobj;
    pf = co2GasLineFitFunction(f.fittedSpectra(ii).x,temp.center,temp.w_g,temp.w_l,...
        temp.a1,temp.a2,0,0,0);
    plot(subplot(2,1,1),f.fittedSpectra(ii).x,pf)
    hold on
end
title('Fitted Spectra')
xlabel('Wavenumbers (cm^{-1})')
ylabel('Absorbance (AU)')
box off
set(gca,'TickDir','out')
hold off

plot(subplot(2,1,2),f.timePts,concOverTime(f),'o-','color','blue');
hold on
title('Concentration Over Time')
xlabel('Time (s)')
ylabel('Concentration (M)')
box off
set(gca,'TickDir','out')
hold off

set(gcf,'Units','normalized')
set(gcf,'Color','w')
set(gcf,'Position',[0.5 0 0.35 1])
%% Save the data
FTIR_peak_fitting_params.sp = sp;
FTIR_peak_fitting_params.trial_spectrum = trial_spectrum;
FTIR_peak_fitting_params.fit_range = range1;
FTIR_peak_fitting_params.spectra_range = spectra_range;
FTIR_peak_fitting_params.file_prefix = file_prefix;
cd(data_path)
save(date_of_experiment + "_FTIR_peak_fitting.mat",'FTIR_peak_fitting_params')
%% ---- END OF STEP 2 ----

%% Step 3: Fit for single diffusion coefficient

%% Guesses for uptake curve fitting, by eye
t = f.timePts;
% t = t(1:80);
y = f.concOverTime;
% y = y(1:80);
A = f.radius;
nmax = 150;
rres = 50;
rlim = 350;
sigma = 704;
dy = 0;

% ---- User input starting values
dx = displacement_mm*1000;  % from the image analysis. can be overridden
%     D    C   
sp = [160  0.244]; % put guess here
ub = [1e5 1e3];
lb = [0 0];
% ----

figure(728);clf
plot(t,y)
hold on
ymodel = diffusion_moving_beam(t,sp(1),A,sp(2),nmax,sigma,dx,dy,"rlim",rlim,'dx definition','from center');
plot(t,ymodel)
res = y(:) - ymodel(:);
plot(t,res-0.025,'ro')
errs(1) = sum((res).^2);
errs(2) = std(res);
fprintf("SSE: " + errs(1) + "; Std. res: " + errs(2) + "\n")
%% Do the uptake curve fitting

%set up options and type
opts = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',lb,'Upper',ub,'StartPoint',sp,...
    'Display','Iter',...
    'Weights',ones(size(y))*1/errs(2),...
    'TolFun',1e-16,...
    'TolX',1e-16);

ft = fittype(@(D,C,t) diffusion_moving_beam(t,D,A,C,nmax,sigma,dx,dy,"rlim",rlim,'dx definition','from center'),...
    'independent',{'t'},...
    'dependent','absorbance',...
    'coefficients',{'D','C'},...
    'options',opts);

%set up structure for storing output
out = struct('x',[],'ydata',[],'yfit',[],'res',[],...
    'fobj',[],'G',[],'O',[]);

tic

%do the fit
[fobj,G,O] = fit(t,y',ft);

toc

%get results
yfit = fobj(t);
out.x = t;
out.ydata = y;
out.yfit = yfit;
out.res = y - yfit;
out.fobj = fobj;
out.G = G;
out.O = O;

if out.O.exitflag < 1
    warning('Curve fit did not converge!!! Results might not be trustworthy.');
end

f.diffusionFitResult = out;
%% display fit result
figure(144);clf

plot(f.diffusionFitResult.x,f.diffusionFitResult.ydata,...
    'o','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue')
hold on
plot(f.diffusionFitResult.x,f.diffusionFitResult.yfit,...
    'red','LineWidth',1.5)
residuals = f.diffusionFitResult.ydata(:) - f.diffusionFitResult.yfit(:);
plot(f.diffusionFitResult.x,(residuals*1 - 0.02),'o','MarkerEdgeColor','red')
legend('Data points','Fitted curve','Location','northwest')
xlabel('Time (s)')
ylabel('Concentration (M)')
hold off


% get confidence intervals
ci = confint(f.diffusionFitResult.fobj);

readout = [string(f.diffusionFitResult.fobj.D)]
others = ["95% Confidence Interval is "+ci(1)+" to "+ci(2)+".";...
    "R^2 = "+string(f.diffusionFitResult.G.rsquare)]

f.diffusionFitResult.fobj
%% Save the data

single_diffusion_fitting_params.sp = sp;
single_diffusion_fitting_params.ub = ub;
single_diffusion_fitting_params.lb = lb;
single_diffusion_fitting_params.dx = dx;
single_diffusion_fitting_params.nmax = nmax;
single_diffusion_fitting_params.rlim = rlim;
single_diffusion_fitting_params.sigma = sigma;
save(date_of_experiment + "_single_diffusion_fitting_params.mat",'single_diffusion_fitting_params')

f.fitMethod = 'diffusion_moving_beam.m';
cd(data_path);
save(date_of_experiment + "_single_diffusion_fitting_params.mat",'single_diffusion_fitting_params')
save(f.dateString + "_single_diffusion_fitting","f")

%% Update lab notebook with results

% ---- PUT IN THE CORRECT NOTEBOOK PAGE TITLE ---
notebook = 'Matt Lab Notebook';
folder = 'Experiments';
page_title = '2025-05-12 Diffusion of CO2 in 75% EMIM NTF2 in PEGDA';
% ----

obj = labarchivesCallObj('notebook',notebook,...
    'folder',folder,...
    'page',page_title);
% microscope photo
figure(6)
obj = obj.updateFigureAttachment('caption','Kiralux camera photo of the sample annotated with calculated values');
% uptake curve
figure(3)
obj = obj.updateFigureAttachment;
% single diffusion fitting result
figure(144)
caption = "Single diffusion coefficient fitting: ";
coeffs = coeffnames(f.diffusionFitResult.fobj);
units = ["um^2/s" "M"];
if numel(units) ~= numel(coeffs)
    error("Cannot match all fitting parameters with a unit.")
end
for ii = 1:numel(coeffs)
   std_devs{ii} = (ci(2,ii) - ci(1,ii))/4;
   caption = caption + coeffs{ii} + " = " + f.diffusionFitResult.fobj.(coeffs{ii})...
       + " ± " + std_devs{ii} + " " + units(ii) + ", ";
end
obj = obj.updateFigureAttachment('caption',caption);
%% ---- END OF STEP 3 ----

%% Step 4: Fit for double diffusion coefficient

% Proceed to the double_diffusion_fitting.m script for this part. You will
% need to have run this entire script for it to work as it will draw from
% the .mat files you saved here.