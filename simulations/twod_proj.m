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

% =========================================================================
% SET SIMULATION CASE
% =========================================================================

beam_type = 'focus';
r = 25e-3;          % radius of the steer [m] 
steering_angle = 30; % angle of steering [deg]

% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 256;           % number of grid points in the x (row) direction
Ny = Nx;            % number of grid points in the y (column) direction
dx = 100e-3/Nx;    	% grid point spacing in the x direction [m]
dy = dx;            % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium, with sound speed and 
% POWER LAW Absorption
medium.sound_speed = 1500;  % [m/s]
medium.alpha_coeff = 0;     % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;
medium.density = 1000;      % [kg/m^3]

% create the time array
kgrid.makeTime(medium.sound_speed);

% define source mask for a linear transducer  
num_elements = 64;      % [grid points]
x_offset = 25;          % [grid points]
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
% DETECTION
% =========================================================================

% from example of Transducer field examples

% create a sensor mask covering the entire computational domain using the
% opposing corners of a rectangle
% TODO change this into a line
% check when the 

% define the first line sensor region by specifying the location of
% opposing corners
rect1_x_start = 1;
rect1_y_start = 1;
rect1_x_end = Nx;
rect1_y_end = Ny;

% assign the list of opposing corners to the sensor mask
sensor.mask = [rect1_x_start, rect1_y_start, rect1_x_end, rect1_y_end].';

% set the record mode to capture a time series of pressure and velocity
sensor.record = {'p_max'};
sensor.record = {'p', 'I'};

% assign the input options
%input_args = { 'PMLInside', false, 'PlotPML', false};
input_args = {'DisplayMask', source.p_mask, 'PlotLayout', false};

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
data = sensor_data;
name = strcat(datestr(datetime('now'),'mmdd'), '_', ...
        beam_type,'_', int2str(steering_angle), '.mat');
save(name, 'data');

%%
%test = sqrt((x_focus - element_spacing * ... 
%             element_index).^2 + z_focus^2)/(medium.sound_speed * kgrid.dt);
             
%test = test - min(test);
%ttt = mod(test, 1/(tone_burst_freq * kgrid.dt));
%plot(ttt);
%hold on;
%plot(test);
