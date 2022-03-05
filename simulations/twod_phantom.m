% To simulate a 2D pressure field
% with an arbitrary transducer response
% frequency, phase,
%
% This example demonstrates how to use k-Wave to steer a tone burst from a
% linear array transducer in 2D. It builds on the Simulating Transducer
% Field Patterns Example.
%
% author: Kenton Kwok
% date: 7/2/2022, last updated 5/3/2022
%
% Code is taken from K-Wave examples


clearvars;
addpath('k-Wave/', 'simulations/')

% Simulation settings
% Set data_cast to single and save memory
DATA_CAST = 'single';

% =========================================================================
% SET SIMULATION CASE
% =========================================================================
beam_type = 'steer';
r = 25e-3;          % radius of the steer [m]


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
% SIMULATION ENVIRONMENT
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
element_spacing = dx;           % [m]
tone_burst_freq = 1e6;          % [Hz]
tone_burst_cycles = 5;


% create an element index relative to the centre element of the transducer
element_index = -(num_elements - 1)/2:(num_elements - 1)/2;


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

% Create a binary sensor mask and make measurements at same location of the
% transducer
sensor.mask = zeros(Nx, Ny);
sensor.mask(x_offset, start_index:start_index + num_elements - 1) = 1;

% Set the record mode to capture a time series of pressure to mimic an
% ultrasound transducer
sensor.record = {'p'};

% =========================================================================
% RUN THE SIMULATION
% =========================================================================

% Range of steering angles to test
steering_angles = -2:2:2;

% Preallocate the storage for the scanned image
number_scan_lines = length(steering_angles);
scan_lines = zeros(number_scan_lines, kgrid.Nt);

% Assign the input options
input_args = {'DisplayMask', source.p_mask, 'PMLInside', false, 'PlotPML', false, ...
    'PlotLayout', true, 'PMLSize', [PML_X_SIZE, PML_Y_SIZE], 'DataCast', DATA_CAST};

% Loop through the range of angles 
for angle_index = 1:number_scan_lines
    
    % Update the command line status
    disp('');
    disp(['Computing scan line ' num2str(angle_index) ' of ' num2str(number_scan_lines)]);

    % update the current steering angle
    steering_angle = steering_angles(angle_index); 
    
    x_focus = r * sind(steering_angle);      % [m]
    z_focus = r * cosd(steering_angle);      % [m]
    
    offset = 70;

    switch beam_type
        % this modifies the tone_burst offsets, note that they are
        % element-indexed
    
        case 'steer'
            % use geometric beam forming to calculate the tone burst offsets for
            % each transducer element based on the element index
            tone_burst_offset = offset + element_spacing * element_index * ...
                sin(steering_angle * pi/180) / (c0 * kgrid.dt);
        case 'steer_wrap'
            % apply a phase wrapping, equivalent to modulo operator
            tone_burst_offset = offset + mod(element_spacing * element_index * ...
                sin(steering_angle * pi/180) / (c0* kgrid.dt), ...
                1/(tone_burst_freq* kgrid.dt)) ;
        case 'focus'
            r = sqrt(z_focus^2 + x_focus^2);
            tone_burst_offset = offset + (r - sqrt((x_focus - element_spacing * ...
                element_index).^2 + z_focus^2))/(c0 * kgrid.dt);
        case 'focus_wrap'
            r = sqrt(z_focus^2 + x_focus^2);
            tone_burst_offset = offset + mod(r - sqrt((x_focus - element_spacing * ...
                element_index).^2 + z_focus^2)/(c0 * ...
                kgrid.dt), 1/(tone_burst_freq* kgrid.dt));
    end


    % create the tone burst signals
    source.p = toneBurst(sampling_freq, tone_burst_freq, tone_burst_cycles, ...
        'SignalOffset', tone_burst_offset);

    % Run the simulation
    sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

    %% Generate Beamformed lines
    % get the current apodization setting
    apodization = ones(num_elements);
    
    % get the current beamforming weights and reverse
    delays = -round(tone_burst_offset);
    
    % offset the received sensor_data by the beamforming delays and
    % apply receive apodization
    for element_index = 1:num_elements
        if delays(element_index) > 0
    
            % shift element data forwards
            sensor_data.p(element_index, :) = apodization(element_index).*[sensor_data.p(element_index, 1 + delays(element_index):end), zeros(1, delays(element_index))];
    
        elseif delays(element_index) < 0
    
            % shift element data backwards
            sensor_data.p(element_index, :) = apodization(element_index).*[zeros(1, -delays(element_index)), sensor_data.p(element_index, 1:end + delays(element_index))];
    
        end
    end
    
    %%
    % form the a-line summing across the elements
    line = sum(sensor_data);

    % Extract the scan line from the sensor data
    scan_lines(angle_index, :) = line;

end

%%
% =========================================================================
% PROCESSING STEPS
% =========================================================================
figure;

% Time-window the sensor_data to not pick up the input signal
% Get the number of time points in the source signal
num_source_time_points = length(source.p(1,:));
p_data = sensor_data.p(1:end, num_source_time_points:end);

% Receive Beamforming
S = sum(p_data, 1);
plot(S)


% Time Gain Compensation

% create time gain compensation function based on attenuation value and
% round trip distance
tgc = exp(tgc_alpha_np_m * 2 * r);

% apply the time gain compensation to each of the scan lines
scan_lines = bsxfun(@times, tgc, scan_lines);

% Frequency Filtering

% Envelope Detection

% Log Compression

% Scan Conversion
%%
% =========================================================================
% VISUALISATION
% =========================================================================
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


