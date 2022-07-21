
% Gerchberg Saxton Algorithm (Fienup improved)
% ------------------------------------------------------------
% Gerchberg Saxton Algorithm
% an example of the GS algorithm is shown in this program
% ------------------------------------------------------------
% pseudo code of GS algorithm
% Gerchberg Saxton Algorithm(Source, Target, Retrieved_Phase)
%  A = IFT(Target)
%  while error criterion is not satisfied
%    B = Amplitude(Source) * exp(i*Phase(A))
%    C = FT(B)
%    D = Amplitude(Target) * exp(i*Phase(C))
%    A = IFT(D)
%  end while
%  Retrieved_Phase = Phase(A)
% ------------------------------------------------------------

%% start

    clc; clearvars; close all;
    tic;

%% test parameters
    
    test_num = 0;
    border_lines = +5;              % positive for white borders, negative for black borders (in pixels)
    flip_result_X = 1;              % 1 for flipping the output hologram, otherwise the projected image would be flipped
    flip_result_Y = 0;
    include_lens = 1;               % 1 for including an output with lens
    iteration_num = 1000;
    iteration_target = 0.999;
    lcos_res_X = 1920;
    lcos_res_Y = 1080;
    lcos_pixel_pitch = 6.4e-6;
    wavelength = 517.8e-9;
    focus_length = 0.5;

%% generate grid
    
    resolution = max( lcos_res_X, lcos_res_Y );
    x = linspace( -lcos_pixel_pitch * resolution/2, lcos_pixel_pitch * ( resolution/2 - 1 ), resolution );
    y = linspace( -lcos_pixel_pitch * resolution/2, lcos_pixel_pitch * ( resolution/2 - 1 ), resolution );
    [X,Y] = meshgrid(x,y);

%% calculate required fresnel lens profile
    
    k = ( 2 * pi ) / wavelength;
    fresnel_lens = exp( -1i * k * ( X.^2 + Y.^2 ) / ( 2 * focus_length ) );

%% calculate input beam and intensity （gaussian beam was selected）
    
	beam_x0 = 0e-3;     		% beam center
	beam_y0 = 0e-3;     		% beam center
	beam_waist = 3e-3;      	% beam waist
	beam_peak = 1;              % peak of the beam
	input_intensity = beam_peak * exp( - ( ( X-beam_x0 ).^2 + ( Y-beam_y0 ).^2 ) / ( 2 * beam_waist^2 ) );
    % surf(input_intensity);
    % shading interp
    % axis tight
    
%% target image input

    target = rgb2gray( imread( [ 'images\' num2str(test_num) '_target.png' ] ) );
    target = imresize( target,[resolution,resolution],'bilinear' );
	target = im2double(target);
    if border_lines > 0
        target(1:border_lines,:) = 1;
        target(resolution-border_lines+1:resolution,:) = 1;
        target(:,1:border_lines) = 1;
        target(:,resolution-border_lines+1:resolution) = 1;
    elseif border_lines < 0
        border_lines = abs(border_lines);
        target(1:border_lines,:) = 0;
        target(resolution-border_lines+1:resolution,:) = 0;
        target(:,1:border_lines) = 0;
        target(:,resolution-border_lines+1:resolution) = 0;
    end
    % ...and transfer to sweet sweet GPU arrays
    if gpuDeviceCount > 0
        fprintf( 'GPU online.\n' );
        target = gpuArray(target);
    end
    
%% correlation logging
    
	correlation = zeros( 1, iteration_num );

%% modified GS algorithm iterations

	rand_phase = 2 * pi * rand(resolution,resolution);
    A = fftshift( ifft2( ifftshift( target .* exp( 1i*rand_phase ) ) ) );
    for iteration_count = 1 : iteration_num
        B = abs( input_intensity ) .* exp( 1i * angle( A ) );
        C = fftshift( fft2( fftshift( B ) ) );
        result_diff = target - abs(C)/max(max(abs(C)));
        correlation(iteration_count) = corr2( abs(target), abs(C) );
        if correlation(iteration_count) > iteration_target
            fprintf( 'Done! correlation: %0.6f\n', correlation(iteration_count) );
            fprintf( 'Total iterations: %i\n', iteration_count );
            break;
        else
            fprintf( 'iterations: %4i/%4i, correlation: %0.6f\n', iteration_count, iteration_num, correlation(iteration_count) );
            D = abs( target + rand * result_diff ) .* exp( 1i * angle( C ) );
            A = fftshift( ifft2( ifftshift( D ) ) );
        end
    end
    
%% extract results

    % reconstruct
    B = abs( input_intensity ) .* exp( 1i * angle( A ) );
    C = fftshift( fft2( fftshift( B ) ) );
    result = abs(C) / max( max( abs(C) ) );
    % phase info
    hologram_phase = angle(A) + pi;
    if flip_result_X == 1
        hologram_phase = flip(hologram_phase,2);
    end
    if flip_result_Y == 1
        hologram_phase = flip(hologram_phase,1);
    end
    hologram_phase_w_lens = angle( exp( 1i * ( angle(A) + angle(fresnel_lens) ) ) ) + pi;
    % Cut to actual resolution of LCoS
    position_top = ( resolution - lcos_res_Y ) / 2 + 1;
    position_bottom = ( resolution + lcos_res_Y ) / 2;
    position_left = ( resolution - lcos_res_X ) / 2 + 1;
    position_right = ( resolution + lcos_res_X ) / 2;
    hologram_display = hologram_phase( position_top:position_bottom, position_left:position_right );
    hologram_display_w_lens = hologram_phase_w_lens( position_top:position_bottom, position_left:position_right );

%% image output as files

    % result
    result_bitmap = uint8( interp1( [0,max(max(result))], [0,255], result, 'linear' ) );
    result_filename = [ 'images\' num2str(test_num) '_result.bmp' ];
    imwrite( result_bitmap, result_filename );
    % hologram w/o lens
    holo_bitmap = uint8( interp1( [0 2*pi], [0 255], hologram_display, 'linear' ) );
    holo_filename = [ 'holograms\' num2str(test_num) '_hologram_none.bmp' ];
    imwrite( holo_bitmap, holo_filename );
    % hologram w lens
    if include_lens == 1 
        holo_lens_bitmap = uint8( interp1( [0 2*pi], [0 255], hologram_display_w_lens, 'linear' ) );
        holo_lens_filename = [ 'holograms\' num2str(test_num) '_hologram_' num2str(1e9*wavelength) 'nm_' num2str(100*focus_length) 'cm.bmp' ];
        imwrite( holo_lens_bitmap, holo_lens_filename );
    end

%% result display

	figure('Name','Result','NumberTitle','off','WindowState','maximized','Color','#7F7F7F')
    % Fresenl Lens
    subplot(2,3,1);
	imagesc( angle(fresnel_lens)+pi )
    colorbar
	title( 'Fresenl Lens' );
    axis image
    % Target Image
    subplot(2,3,2);
	imshow(target,[])
    colorbar
	title([ 'Target Image ' num2str(test_num) ]);
    axis image
    % Reconstructed Image
	subplot(2,3,3);
	imshow(result,[])
    colorbar
	title([ 'Reconstructed Image ' num2str(test_num) ]);
    axis image
    % Correlation
    subplot(2,3,4);
	i = 0:1:iteration_count;
	plot( i, [0 correlation(1:iteration_count)] );
	title( 'Correlation' );
    axis([0 iteration_count 0 1])
    % Hologram
	subplot(2,3,5);
	imagesc(hologram_phase/pi)
    colorbar
	title([ 'Hologram ' num2str(test_num) ]);
    axis image
    % Hologram with Lens
    subplot(2,3,6);
	imagesc(hologram_phase_w_lens/pi)
    colorbar
	title([ 'Hologram ' num2str(test_num) ' with Lens' ]);
    axis image
	
%% end
    
    toc;
    