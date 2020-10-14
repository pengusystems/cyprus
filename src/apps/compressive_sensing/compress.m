% Settings.
clc;
n_cs_samples = 50;
sigma = 1e-4;
reconstruction_method = 'MP';
example = 3;

% Synthsize an interferogram(time) for a one second interval of a periodic function.
% Also create a sub-sampled version to compare with.
fs = 1000;
n_samples = fs + 1;
t = (0:1/fs:1)';
f = (-fs/2:fs/2)';
t_subsampled = (0:1/n_cs_samples:1)';
f_subsampled = (-n_cs_samples/2:n_cs_samples/2)';
switch example
    case 1
        interferogram_handle = @(t) cos(2*pi*100*t) + cos(2*pi*30*t) + cos(2*pi*10*t);
        interferogram = interferogram_handle(t);
        interferogram_subsampled = interferogram_handle(t_subsampled);

    case 2
        interferogram_handle = @(t) exp(-5*t) .* cos(2*pi*100*t + pi/2);
        interferogram = interferogram_handle(t);
        interferogram_subsampled = interferogram_handle(t_subsampled);
 
    case 3
        ftir_handle = @(f) cauchypdf(f,30,1) + cauchypdf(f,-30,1);
        interferogram_handle = @(f) real(ifft((ftir_handle(f))));
        interferogram = interferogram_handle(f);
        interferogram_subsampled = interferogram_handle(f_subsampled);
    
    case 4     
        ftir_handle = @(f) 0.5*(f==100) + 0.5*(f==-100);% +  5*(f > 300) + 5*(f < -300);
        %ftir_handle = @(f) cauchypdf(f,30,1) + cauchypdf(f,-30,1) + 5*(f < -100) + 5*(f > 100);
        interferogram_handle = @(f) real(ifft((ftir_handle(f))));
        plot(interferogram_handle(f));
        exit;
        fp = 100;
        d = fdesign.lowpass('Fp,Fst,Ap,Ast',fp,fp+30,0.5,50,fs);
        Hd = design(d);        
        interferogram = filtfilt(Hd.Numerator,1,interferogram_handle(f));
        interferogram_subsampled = interferogram_handle(f_subsampled);        
end

% Find the spectrum for the full and sub-sampled interferograms.
spectra = fft(interferogram);
subsampled_fft = fftshift(fft(interferogram_subsampled));

% Create the measurement matrix.
measurement_matrix = zeros(n_cs_samples, n_samples);
samples_permutation = randperm(n_samples);
measurement_cs_sub_ind_ = [(1:n_cs_samples)' samples_permutation(1:n_cs_samples)'];
measurement_cs_lin_ind = sub2ind(size(measurement_matrix), measurement_cs_sub_ind_(:,1), measurement_cs_sub_ind_(:,2));
measurement_matrix(measurement_cs_lin_ind) = 1;

% Obtain measurements
measurements = measurement_matrix * interferogram + sigma * randn(n_cs_samples,1);

% Reconstruct.
reconstruction = reconstruct(measurements, measurement_matrix, reconstruction_method);

% Calculate MSE between original and reconstructed.
mse = norm(spectra - reconstruction)^2 / length(spectra);

% Compute masked data for plotting.
masked_samples = samples_permutation(n_cs_samples+1:end);
interferogram_masked = interferogram;
interferogram_masked(masked_samples) = 0;

% Plot results.
figure;
subplot(3,3,1); plot(t, interferogram); title('Original');
subplot(3,3,4); plot(f, abs(fftshift(spectra))); axes_f = gca;
subplot(3,3,7); plot(t, interferogram); axes_t = gca;
subplot(3,3,2); plot(t_subsampled, interferogram_subsampled); title('Uniform Sub-sampled'); axis([axes_t.XLim, axes_t.YLim]);
subplot(3,3,5); plot(f_subsampled, abs(fftshift(subsampled_fft))); %axis([axes_f.XLim, axes_f.YLim]);
subplot(3,3,8); plot(t_subsampled, real(ifft(subsampled_fft))); axis([axes_t.XLim, axes_t.YLim]);
subplot(3,3,3); plot(t, interferogram_masked); title('CS'); axis([axes_t.XLim, axes_t.YLim]);
subplot(3,3,6); plot(f, abs(fftshift(reconstruction))); %axis([axes_f.XLim, axes_f.YLim]);
subplot(3,3,9); plot(t, real(ifft(reconstruction))); axis([axes_t.XLim, axes_t.YLim]);
