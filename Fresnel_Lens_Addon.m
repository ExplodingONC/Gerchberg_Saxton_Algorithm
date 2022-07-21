
% Add fresnel lens to existing holograms
% ------------------------------------------------------------

% %-------------------------------------------------------
% start

    clc; clearvars; close all;
    tic;

% %-------------------------------------------------------
% test parameters
    
    lens_only = 0;          % no base hologram
    test_num = 8;
    lcos_res_X = 1920;      % default value
    lcos_res_Y = 1080;      % default value
    lcos_pixel_pitch = 6.4e-6;
    wavelength = 520e-9;
    focus_length = 2.0;

% %-------------------------------------------------------
% hologram image input

    if lens_only == 0
        input = imread( [ 'holograms\' num2str(test_num) '_hologram_none.bmp' ] );
        input = im2double(input);
        lcos_res_X = size(input,2);
        lcos_res_Y = size(input,1);
    else
        input = zeros( lcos_res_Y, lcos_res_X );
    end
    input = input * 2*pi - pi;
    
% %-------------------------------------------------------
% generate grid
    
    x = linspace( -lcos_pixel_pitch * lcos_res_X / 2, lcos_pixel_pitch * lcos_res_X / 2, lcos_res_X );
    y = linspace( -lcos_pixel_pitch * lcos_res_Y / 2, lcos_pixel_pitch * lcos_res_Y / 2, lcos_res_Y );
    [X,Y] = meshgrid(x,y);

% %-------------------------------------------------------
% calculate required fresnel lens profile
    
    k = ( 2 * pi ) / wavelength;
    fresnel_lens = exp( -1i * k * ( X.^2 + Y.^2 ) / ( 2 * focus_length ) );
    
% %-------------------------------------------------------
% lens add-on

    output = angle( exp( 1i * ( input + angle(fresnel_lens) ) ) ) + pi;
    
% %-------------------------------------------------------
% image output as files

    holo_lens_bitmap = uint8( interp1( [0 2*pi], [0 255], output, 'linear' ) );
    if lens_only == 0
        holo_lens_filename = [ 'holograms\' num2str(test_num) '_hologram_' num2str(1e9*wavelength) 'nm_' num2str(100*focus_length) 'cm.bmp' ];
    else
        holo_lens_filename = [ 'holograms\fresnel_lens_' num2str(1e9*wavelength) 'nm_' num2str(100*focus_length) 'cm.bmp' ];
    end
    imwrite( holo_lens_bitmap, holo_lens_filename );
    
% %-------------------------------------------------------
% end
    
    toc;
    