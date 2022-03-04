% Defining An Ultrasound Transducer Example
%
% This file deals with the creation of a B-mode signal
%
%
% The creation of a kWaveTransducer object will only work in versions of
% MATLAB recent enough to support user defined classes. 
%
% author: Bradley Treeby
% date: 20th July 2011
% last update: 11th June 2017
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2011-2017 Bradley Treeby

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>. 

clearvars;

% simulation settings
% set data_cast to single and save memory
DATA_CAST = 'single';

% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set the size of the perfectly matched layer (PML)
PML_X_SIZE = 20;            % [grid points]
PML_Y_SIZE = 10;            % [grid points]

% set total number of grid points not including the PML
Nx = 128 - 2*PML_X_SIZE;    % [grid points]
Ny = 128 - 2*PML_Y_SIZE;    % [grid points]

% set desired grid size in the x-direction not including the PML
x = 40e-3;                  % [m]

% calculate the spacing between the grid points
dx = x/Nx;                  % [m]
dy = dx;                    % [m]

% create the k-space grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% define the properties of the propagation medium
%medium.sound_speed = 1540;      % [m/s]
%medium.density = 1000;          % [kg/m^3]
c0 = 1540;                      % [m/s]
rho0 = 1000;                    % [kg/m^3]
medium.alpha_coeff = 0.75;      % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;
medium.BonA = 6;

% create the time array
t_end = 100e-6;                  % [s]
kgrid.makeTime(c0, [], t_end);

% =========================================================================
% DEFINE THE ULTRASOUND TRANSDUCER
% =========================================================================

% define source mask for a linear transducer  
num_elements = 72;      % total number of transducer elements [grid points]
element_spacing = dx;   % [m]
x_offset = 1;           % [grid points]

% PML is outside the computational grid, so we do not need to offset
source.p_mask = zeros(Nx, Ny);
start_index = Ny/2 - round(num_elements/2) + 1;
source.p_mask(x_offset, start_index:start_index + num_elements - 1) = 1;

% =========================================================================
% DEFINE THE INPUT SIGNAL
% =========================================================================

% define properties of the input signal
source_strength = 1e6;          % [Pa]
tone_burst_freq = 0.5e6;        % [Hz]
tone_burst_cycles = 5;

% define the properties of the tone burst used to drive the transducer
sampling_freq = 1/kgrid.dt;     % [Hz]

% define the parameters used for steering the beam
steering_angle = 10;            % [Deg]
r = 20e-3;                      % [m]
x_focus = r * sind(steering_angle);      % [m]
z_focus = r * cosd(steering_angle);      % [m]

% create an element index relative to the centre element of the transducer
element_index = -(num_elements - 1)/2:(num_elements - 1)/2;

% use geometric beam forming to calculate the tone burst offsets for 
% each transducer element based on the element index

offset = 0; % offset to make sure no negative values are present [time grid points]
tone_burst_offset = offset + element_spacing * element_index * ...
            sin(steering_angle * pi/180) / (c0 * kgrid.dt) * 0;

% create the tone burst signals
source.p = toneBurst(sampling_freq, tone_burst_freq, tone_burst_cycles, ...
    'SignalOffset', tone_burst_offset);

%% 

% create the input signal using toneBurst 
%input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% scale the source magnitude by the source_strength divided by the
% impedance (the source is assigned to the particle velocity)
%input_signal = (source_strength ./ (medium.sound_speed * medium.density)) .* input_signal;

% calculate the width of the transducer in grid points
%transducer_width = transducer.number_elements * transducer.element_width ...
    %+ (transducer.number_elements - 1) * transducer.element_spacing;

% use this to position the transducer in the middle of the computational grid
%transducer.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);

% properties used to derive the beamforming delays
%transducer.focus_distance = 20e-3;              % focus distance [m]
%transducer.elevation_focus_distance = 19e-3;    % focus distance in the elevation plane [m]

% apodization
%transducer.transmit_apodization = 'Rectangular';    
%transducer.receive_apodization = 'Rectangular';

% print out transducer properties
%transducer.properties;

% =========================================================================
% DEFINE SENSOR MASK
% =========================================================================

% create a binary sensor mask 
sensor.mask = zeros(Nx, Ny);
sensor.mask(x_offset, start_index:start_index + num_elements - 1) = 1;

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
scattering_c0 = c0 + 25 + 75 * scattering_map;
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
%
scattering_region1 = makeDisc(Nx, Ny, round(x_pos / dx), round(y_pos / dx), radius);

% assign region
sound_speed_map(scattering_region1 == 1) = scattering_c0(scattering_region1 == 1);
density_map(scattering_region1 == 1) = scattering_rho0(scattering_region1 == 1);

% assign to the medium inputs
medium.sound_speed = sound_speed_map;
medium.density = density_map;

% =========================================================================
% RUN THE SIMULATION
% =========================================================================

% set the input settings

input_args = {'DisplayMask', source.p_mask, 'PMLInside', false, 'PlotLayout', true, ... 
    'PMLSize', [PML_X_SIZE, PML_Y_SIZE], 'DataCast', DATA_CAST };

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% calculate the amplitude spectrum of the input signal and the signal
% recorded each of the sensor positions 
% [f_input, as_input] = spect([input_signal, zeros(1, 2 * length(input_signal))], 1/kgrid.dt);
% [~, as_1] = spect(sensor_data(1, :), 1/kgrid.dt);
% [~, as_2] = spect(sensor_data(2, :), 1/kgrid.dt);
% [f, as_3] = spect(sensor_data(3, :), 1/kgrid.dt);

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the input signal and its frequency spectrum
% figure;
% subplot(2, 1, 1);
% plot((0:kgrid.dt:(length(input_signal) - 1) * kgrid.dt) * 1e6, input_signal, 'k-');
% xlabel('Time [\mus]');
% ylabel('Particle Velocity [m/s]');
% 
% subplot(2, 1, 2);
% plot(f_input .* 1e-6, as_input./max(as_input(:)), 'k-');
% hold on;
% line([tone_burst_freq, tone_burst_freq] .* 1e-6, [0 1], 'Color', 'k', 'LineStyle', '--');
% xlabel('Frequency [MHz]');
% ylabel('Amplitude Spectrum [au]');
% f_max = medium.sound_speed / (2 * dx);
% set(gca, 'XLim', [0, f_max .* 1e-6]);
% 
% % plot the recorded time series
% figure;
% stackedPlot(kgrid.t_array * 1e6, {'Sensor Position 1', 'Sensor Position 2', 'Sensor Position 3'}, sensor_data);
% xlabel('Time [\mus]');
% 
% % plot the corresponding amplitude spectrums
% figure;
% plot(f .* 1e-6, as_1 ./ max(as_1(:)), 'k-', ...
%      f .* 1e-6, as_2 ./ max(as_1(:)), 'b-', ...
%      f .* 1e-6, as_3 ./ max(as_1(:)), 'r-');
% legend('Sensor Position 1', 'Sensor Position 2', 'Sensor Position 3');
% xlabel('Frequency [MHz]');
% ylabel('Normalised Amplitude Spectrum [au]');
% f_max = medium.sound_speed / (2 * dx);
% set(gca, 'XLim', [0, f_max .* 1e-6]);