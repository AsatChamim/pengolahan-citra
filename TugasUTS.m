pkg load image;

disp('--- Memulai Praktikum Pengolahan Citra ---');

% 1. Konversi Citra: RGB, Grayscale, Biner
img_rgb = imread('E:\pengolahan citra\rgb.jpeg');
figure('Name', '1. RGB, Grayscale, Biner');
subplot(1,3,1); imshow(img_rgb); title('Citra Asli (RGB)');
img_gray = rgb2gray(img_rgb);
subplot(1,3,2); imshow(img_gray); title('Grayscale');
img_biner = img_gray > 128;
subplot(1,3,3); imshow(img_biner); title('Biner (Threshold 128)');

% 2. Kuantisasi Grayscale (2, 4, 8 Tingkat)
img_q_orig = imread('E:\pengolahan citra\abu-abu.jpeg');
figure('Name', '2. Kuantisasi Grayscale');
subplot(2, 2, 1); imshow(img_q_orig); title('Citra Grayscale Asli');
img_q2 = uint8(round(img_q_orig / 255) * 255);
subplot(2, 2, 2); imshow(img_q2); title('2 Tingkat Keabuan');
img_q4 = uint8(64 * round(img_q_orig / 64));
subplot(2, 2, 3); imshow(img_q4); title('4 Tingkat Keabuan');
img_q8 = uint8(32 * round(img_q_orig / 32));
subplot(2, 2, 4); imshow(img_q8); title('8 Tingkat Keabuan');

% 3. Meningkatkan Kecerahan Citra
img_gelap = imread('E:\pengolahan citra\gambar gelap.jpeg');
img_gelap_gray = rgb2gray(img_gelap);
figure('Name', '3. Peningkatan Kecerahan Citra');
subplot(2, 2, 1); imshow(img_gelap_gray); title('Citra Asli (Gelap)');
subplot(2, 2, 2); imhist(img_gelap_gray); title('Histogram Citra Asli');
img_bright = img_gelap_gray + 50;
img_bright(img_bright > 255) = 255;
subplot(2, 2, 3); imshow(uint8(img_bright)); title('Kecerahan Ditingkatkan');
subplot(2, 2, 4); imhist(uint8(img_bright)); title('Histogram Ditingkatkan');

% 4. Ekualisasi Histogram
img_eq_orig = imread('E:\pengolahan citra\gambar gelap.jpeg');
img_eq_gray = rgb2gray(img_eq_orig);
figure('Name', '4. Ekualisasi Histogram');
subplot(2, 2, 1); imshow(img_eq_gray); title('Citra Asli');
subplot(2, 2, 2); imhist(img_eq_gray); title('Histogram Citra Asli');
img_eq_result = histeq(img_eq_gray);
subplot(2, 2, 3); imshow(img_eq_result); title('Ekualisasi Histogram');
subplot(2, 2, 4); imhist(uint8(img_eq_result)); title('Histogram Ekualisasi');

% 5. Filter Median dan Filter Rata-Rata (Mean)
img_filter_orig = imread('E:\pengolahan citra\bunga-hibiscus-atau-kembang-sepatu_169.jpeg');
if size(img_filter_orig, 3) == 3
    img_filter_disp = rgb2gray(img_filter_orig);
else
    img_filter_disp = img_filter_orig;
end
figure('Name', '5. Filter Median dan Rata-Rata');
subplot(1, 3, 1); imshow(img_filter_disp); title('Citra Asli (Grayscale)');
img_median = medfilt2(img_filter_disp, [3 3]);
subplot(1, 3, 2); imshow(img_median); title('Filter Median');
h_mean = fspecial('average', [3 3]);
img_mean = imfilter(img_filter_disp, h_mean);
subplot(1, 3, 3); imshow(img_mean); title('Filter Rata-Rata');

% 6. Filter High-Boost (untuk citra warna)
img_daun_color = imread('E:\pengolahan citra\daun.jpeg');
figure('Name', '6. Filter High-Boost (Warna)');
subplot(1, 2, 1); imshow(img_daun_color); title('Citra Asli (Warna)');
sigma_gauss = 1; kernel_size_gauss = 3;
x_gauss = -floor(kernel_size_gauss/2):floor(kernel_size_gauss/2);
g_gauss = exp(-(x_gauss.^2)/(2*sigma_gauss^2)); g_gauss = g_gauss/sum(g_gauss);
h_gaussian = g_gauss'*g_gauss;
A_highboost = 1.5;
img_highboost = zeros(size(img_daun_color));
for i = 1:3
    channel = double(img_daun_color(:,:,i));
    smoothed_channel = conv2(channel, h_gaussian, 'same');
    highboost_channel = channel * A_highboost - (A_highboost - 1) * smoothed_channel;
    highboost_channel(highboost_channel < 0) = 0;
    highboost_channel(highboost_channel > 255) = 255;
    img_highboost(:,:,i) = highboost_channel;
end
subplot(1, 2, 2); imshow(uint8(img_highboost)); title(['High-Boost (A = ' num2str(A_highboost) ')']);

% 7. Rotasi dan Perbesaran Citra (untuk citra warna)
img_slanted = imread('E:\pengolahan citra\miring.jpeg');
if size(img_slanted, 3) == 1
    img_slanted = cat(3, img_slanted, img_slanted, img_slanted);
end
angle = 30;
img_rotated = imrotate(img_slanted, angle, 'bilinear', 'crop');
scale_factor = 2;
img_enlarged = imresize(img_rotated, scale_factor, 'bilinear');
figure('Name', '7. Rotasi dan Perbesaran Citra (Warna)');
subplot(1,3,1); imshow(img_slanted); title('Citra Asli (Warna)');
subplot(1,3,2); imshow(img_rotated); title(['Rotasi ' num2str(angle) '\circ']);
subplot(1,3,3); imshow(img_enlarged); title(['Perbesaran ' num2str(scale_factor) 'x']);
%8. menerapkan effek twirl
function out_twirl = twirl_image_local(img_in_double, strength, radius)
    [rows, cols, channels] = size(img_in_double);
    [x, y] = meshgrid(1:cols, 1:rows);
    center = [cols/2, rows/2];
    xc = x - center(1);
    yc = y - center(2);
    angle = strength * exp(-(xc.^2 + yc.^2)/(2*radius^2));
    xn = xc.*cos(angle) - yc.*sin(angle) + center(1);
    yn = xc.*sin(angle) + yc.*cos(angle) + center(2);
    out_twirl = zeros(rows, cols, channels);
    for k = 1:channels
        out_twirl(:,:,k) = interp2(img_in_double(:,:,k), xn, yn, 'linear', 0);
    end
    out_twirl = uint8(out_twirl);
end
img_wajah = imread('E:\pengolahan citra\foto aku.jpeg');
img_twirl_result = twirl_image_local(double(img_wajah), 100.0, min(size(img_wajah))/200);
figure('Name', '8. Efek Twirl');
subplot(1,2,1), imshow(img_wajah), title('Asli');
subplot(1,2,2), imshow(img_twirl_result), title('Efek Twirl Kuat');

% 9. Filter Gaussian untuk Pengaburan
img_street = imread('E:\pengolahan citra\gambar jalan.jpeg');
img_street_double = double(img_street);
sigma_gaussian = 15;
kernel_size_gaussian = 2*ceil(3*sigma_gaussian)+1;
if mod(kernel_size_gaussian, 2) == 0
    kernel_size_gaussian = kernel_size_gaussian + 1;
end
h_gaussian_filter = fspecial('gaussian', kernel_size_gaussian, sigma_gaussian);
img_blurred = zeros(size(img_street_double));
for k = 1:size(img_street_double, 3)
    img_blurred(:,:,k) = imfilter(img_street_double(:,:,k), h_gaussian_filter, 'replicate');
end
img_blurred = uint8(img_blurred);
figure('Name', '9. Gaussian Blur');
subplot(1,2,1), imshow(img_street), title('Asli');
subplot(1,2,2), imshow(img_blurred), title(['Gaussian Blur (Sigma = ' num2str(sigma_gaussian) ')']);

% 10. Transformasi Affine
img_affine_orig = imread('E:\pengolahan citra\kotak.jpeg');
img_affine_double = double(img_affine_orig);
[rows_orig, cols_orig, channels_orig] = size(img_affine_double);

[x_orig, y_orig] = meshgrid(1:cols_orig, 1:rows_orig);

theta = pi/6;
s_scale = 1.2;
shx = 0.2;
tx = 50;
ty = -30;

T_rotation_scale = [s_scale*cos(theta), -s_scale*sin(theta); s_scale*sin(theta), s_scale*cos(theta)];
T_shear = [1, shx; 0, 1];
T_affine = T_shear * T_rotation_scale;
T_affine_inv = inv(T_affine);

corners_orig = [1 cols_orig cols_orig 1; 1 1 rows_orig rows_orig];
transformed_corners_x = T_affine(1,1)*corners_orig(1,:) + T_affine(1,2)*corners_orig(2,:) + tx;
transformed_corners_y = T_affine(2,1)*corners_orig(1,:) + T_affine(2,2)*corners_orig(2,:) + ty;

min_x_out = floor(min(transformed_corners_x)); max_x_out = ceil(max(transformed_corners_x));
min_y_out = floor(min(transformed_corners_y)); max_y_out = ceil(max(transformed_corners_y));

new_cols = max_x_out - min_x_out + 1;
new_rows = max_y_out - min_y_out + 1;

[x_out, y_out] = meshgrid(min_x_out:max_x_out, min_y_out:max_y_out);

coords_src_temp = T_affine_inv * [x_out(:)' - tx; y_out(:)' - ty];
x_src = reshape(coords_src_temp(1,:), size(x_out));
y_src = reshape(coords_src_temp(2,:), size(y_out));

img_affine_result = zeros(new_rows, new_cols, channels_orig);
for k = 1:channels_orig
    img_affine_result(:,:,k) = interp2(x_orig, y_orig, img_affine_double(:,:,k), x_src, y_src, 'linear', 0);
end
img_affine_result = uint8(img_affine_result);

figure('Name', '10. Transformasi Affine');
subplot(1,2,1), imshow(img_affine_orig), title('Citra Asli');
subplot(1,2,2), imshow(img_affine_result), title('Transformasi Affine');

disp('--- Semua proses pengolahan citra telah selesai! ---');
