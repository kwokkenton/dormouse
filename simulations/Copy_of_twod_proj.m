% To simulate a 2D pressure field
% with an arbitrary transducer response
% frequency, phase, 
%
% This example demonstrates how to use k-Wave to steer a tone burst from a
% linear array transducer in 2D. It builds on the Simulating Transducer
% Field Patterns Example.
%
% author: Kenton Kwok
% date: 7/2/2022, last updated 28/2/2022


clearvars;
addpath('k-Wave/', 'simulations/')

% simulation settings
% set data_cast to single and save memory
DATA_CAST = 'single';

% =========================================================================
% SET SIMULATION CASE
% =========================================================================

beam_type = 'steer';
r = 25e-3;          % radius of the steer [m] 
steering_angle = 0; % angle of steering [deg]

% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set the size of the perfectly matched layer (PML)
PML_X_SIZE = 20;            % [grid points]
PML_Y_SIZE = 10;            % [grid points]

% set total number of grid points not including the PML
Nx = 256 - 2*PML_X_SIZE;    % [grid points]
Ny = 256 - 2*PML_Y_SIZE;    % [grid points]

% set desired grid size in the x-direction not including the PML
x = 100e-3;                  % [m]

% calculate the spacing between the grid points
dx = x/Nx;                  % grid point spacing in the x direction [m]
dy = dx;                    % grid point spacing in the y direction [m]

% create the k-space grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% =========================================================================
% SIMULATION
% =========================================================================


% define the properties of the propagation medium, with sound speed and 
% POWER LAW Absorption
medium.sound_speed = 1540;  % [m/s]
c0 = 1540; 
medium.density = 1000;      % [kg/m^3]
rho0 = 1000;                    % [kg/m^3]
medium.alpha_coeff = 0.75;     % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;

%medium.BonA = 6;

% create the time array
kgrid.makeTime(medium.sound_speed);

% define source mask for a linear transducer  
num_elements = 64;      % [grid points]
x_offset = 1;          % [grid points]
source.p_mask = zeros(Nx, Ny);
start_index = Ny/2 - round(num_elements/2) + 1;
source.p_mask(x_offset, start_index:start_index + num_elements - 1) = 1;

% define the properties of the tone burst used to drive the transducer
sampling_freq = 1/kgrid.dt;     % [Hz]
x_focus = r * sind(steering_angle);      % [m]
z_focus = r * cosd(steering_angle);      % [m]
element_spacing = dx;           % [m]
tone_burst_freq = 1e6;          % [Hz]
tone_burst_cycles = 5;

% create an element index relative to the centre element of the transducer
element_index = -(num_elements - 1)/2:(num_elements - 1)/2;

offset = 70;

switch beam_type
    % this modifies the tone_burst offsets, note that they are
    % element-indexed
    
    case 'steer'
        % use geometric beam forming to calculate the tone burst offsets for 
        % each transducer element based on the element index
        tone_burst_offset = offset + element_spacing * element_index * ...
            sin(steering_angle * pi/180) / (medium.sound_speed * kgrid.dt);
    case 'steer_wrap'
%         apply a phase wrapping, equivalent to modulo operator
        tone_burst_offset = offset + mod(element_spacing * element_index * ...
            sin(steering_angle * pi/180) / (medium.sound_speed* kgrid.dt), ... 
               1/(tone_burst_freq* kgrid.dt)) ;
    case 'focus'
        r = sqrt(z_focus^2 + x_focus^2);
        tone_burst_offset = offset + (r - sqrt((x_focus - element_spacing * ... 
             element_index).^2 + z_focus^2))/(medium.sound_speed * kgrid.dt);
    case 'focus_wrap'
        r = sqrt(z_focus^2 + x_focus^2);
        tone_burst_offset = offset + mod(r - sqrt((x_focus - element_spacing * ... 
             element_index).^2 + z_focus^2)/(medium.sound_speed * ...
             kgrid.dt), 1/(tone_burst_freq* kgrid.dt));
end 


% create the tone burst signals
source.p = toneBurst(sampling_freq, tone_burst_freq, tone_burst_cycles, ...
    'SignalOffset', tone_burst_offset);

% =========================================================================
% DEFINE THE MEDIUM PROPERTIES
% =========================================================================

% define a random distribution of scatterers for the medium
%background_map_mean = 1;
%background_map_std = 0.008;
%background_map = background_map_mean + background_map_std * randn([Nx, Ny]);

% define a random distribution of scatterers for the highly scattering
% region
scattering_map = randn([Nx, Ny]);
%scattering_c0 = c0 + 25 + 75 * scattering_map;
scattering_c0 = c0 + 50 + 0* scattering_map;
scattering_c0(scattering_c0 > 1600) = 1600;
scattering_c0(scattering_c0 < 1400) = 1400;
scattering_rho0 = scattering_c0 / 1.5;

% define properties
sound_speed_map = c0 * ones(Nx, Ny); %.* background_map;
density_map = rho0 * ones(Nx, Ny); %.* background_map;

% define a sphere for a highly scattering region
radius = 10e-3;
x_pos = 0*32e-3;
y_pos = 0*dy * Ny/2;

scattering_region1 = makeDisc(Nx, Ny, round(x_pos / dx), round(y_pos / dx), round(radius/dx));

% assign region
sound_speed_map(scattering_region1 == 1) = scattering_c0(scattering_region1 == 1);
density_map(scattering_region1 == 1) = scattering_rho0(scattering_region1 == 1);

% assign to the medium inputs
medium.sound_speed = sound_speed_map;
medium.density = density_map;


% =========================================================================
% DETECTION
% =========================================================================

% create a binary sensor mask and make measurements at same location of the
% transducer
sensor.mask = zeros(Nx, Ny);
sensor.mask(x_offset, start_index:start_index + num_elements - 1) = 1;

% set the record mode to capture a time series of pressure to mimic an
% ultrasound transducer
sensor.record = {'p'};

% assign the input options

input_args = {'DisplayMask', source.p_mask, 'PMLInside', false, 'PlotPML', false, ...
    'PlotLayout', true, 'PMLSize', [PML_X_SIZE, PML_Y_SIZE], 'DataCast', DATA_CAST};

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% =========================================================================
% VISUALISATION
% =========================================================================

% VISUALISATIONS FOR THE DELAYS

% get the number of time points in the source signal
% num_source_time_points = length(source.p(1,:));
% 
% % get suitable scaling factor for plot axis
% [~, scale, prefix] = scaleSI(kgrid.t_array(num_source_time_points));
% 
% % plot the input time series
% figure;
% stackedPlot(kgrid.t_array(1:num_source_time_points) * scale, source.p);
% xlabel(['Time [' prefix 's]']);
% ylabel('Input Signals');

%%
% figure;
% plot(element_index, tone_burst_offset*kgrid.dt);
% xlabel('Element Index');
% ylabel('Time delay')
%%
% VISUALISATION FOR THE FINAL PATTERN
% plot the simulated sensor data
% figure;
% 
% mx = max(abs(sensor_data.Ix), [], 'all');
% imagesc(max(sensor_data.Ix, [], 3), [-mx, mx]);
% colormap(getColorMap);
% ylabel('Sensor Position');
% xlabel('Time Step');
% colorbar;

%%
%data = sensor_data;
%name = strcat(datestr(datetime('now'),'mmdd'), '_', ...
%        beam_type,'_', int2str(steering_angle), '.mat');
%save(name, 'data');

%%
%test = sqrt((x_focus - element_spacing * ... 
%             element_index).^2 + z_focus^2)/(medium.sound_speed * kgrid.dt);
             
%test = test - min(test);
%ttt = mod(test, 1/(tone_burst_freq * kgrid.dt));
%plot(ttt);
%hold on;
%plot(test);
