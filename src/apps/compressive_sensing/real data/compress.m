clc;

% Load the data.
load 'data.mat';
load 'background.mat';

% Calculate the full spectra after subtracting the background.
full_spectra.cm = data.cm;
full_spectra.rho = data.rho ./ background.rho;
full_spectra.theta = unwrap(data.theta - background.theta);

% Settings for compressed sensing.
offset = 6000;
range = 11000 - offset;
n_samples = length(data.interf);
n_cs_samples = 200;
sigma = 1e-4;
reconstruction_method = 'MP';

% Create the measurement matrix.
measurement_matrix = zeros(n_cs_samples, n_samples);
samples_permutation = offset + randperm(range, n_cs_samples);
measurement_cs_sub_ind = [(1:n_cs_samples)' samples_permutation'];
measurement_cs_lin_ind = sub2ind(size(measurement_matrix), measurement_cs_sub_ind(:,1), measurement_cs_sub_ind(:,2));
measurement_matrix(measurement_cs_lin_ind) = 1;

% Obtain measurements.
measurements = measurement_matrix * (data.interf - background.interf) + sigma * randn(n_cs_samples,1);

% Reconstruct.
reconstruction_spectra = fftshift(reconstruct(measurements, measurement_matrix, reconstruction_method));

% Populate cs_spectra.
cs_spectra.cm = full_spectra.cm;
cs_spectra.re = real(reconstruction_spectra);
cs_spectra.im = imag(reconstruction_spectra);
[cs_spectra.theta, cs_spectra.rho] = magphase(cs_spectra.re, cs_spectra.im);

% Not sure why that is done in the original code.
point = find(full_spectra.cm >= 1900,1);
full_spectra.theta = full_spectra.theta - full_spectra.theta(point) - 0.1;
cs_spectra.theta = cs_spectra.theta - full_spectra.theta(point) - 0.1;

% Plot results.
plotrange=[750 1840];
figure;
subplot(2,2,3); plot(full_spectra.cm, full_spectra.rho);
ylabel('magnitude (AU)');
xlim(plotrange);
subplot(2,2,1); plot(full_spectra.cm, -full_spectra.theta * 180 / pi);
ylabel('-phase (deg.)');
xlim(plotrange);
title('Full spectra');
subplot(2,2,4); plot(cs_spectra.cm, cs_spectra.rho);
xlim(plotrange);
subplot(2,2,2); plot(cs_spectra.cm, -cs_spectra.theta * 180 / pi);
xlim(plotrange);
title('CS spectra');