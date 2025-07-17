%Menggunakan filter Batas%
pkg load image;

F= imread('E:\pengolahan citra/bunga-hibiscus-atau-kembang-sepatu_169.jpeg');
if size(F, 3) == 3
    F = rgb2gray(F); % Ubah ke grayscale jika RGB
end

figure(1); clf;

% Gambar asli di kolom pertama
subplot(1, 6, 1);
imshow(F);
title('Gambar Asli');

% Tambahkan noise dan tampilkan di kolom 2-6
for i = 1:5
    switch i
        case 1
            F_noise = imnoise(F, 'gaussian');
            noise_type = 'Gaussian';
        case 2
            F_noise = imnoise(F, 'salt & pepper', 0.05);
            noise_type = 'Salt & Pepper';
        case 3
            F_noise = imnoise(F, 'speckle');
            noise_type = 'Speckle';
        case 4
            F_noise = imnoise(F, 'poisson');
            noise_type = 'Poisson';
        case 5
            noise = uint8(20 * randn(size(F)));
            F_noise = im2uint8(mat2gray(double(F) + double(noise)));
            noise_type = 'Uniform';
    end

    subplot(1, 6, i + 1);
    imshow(F_noise);
    title(sprintf('Noise: %s', noise_type));
end


%Menggunkana Filter Median%

pkg load image;
% Membaca gambar asli
F = imread('E:\pengolahan citra\foto aku.jpeg');
if size(F, 3) == 3
    F = rgb2gray(F); % Ubah ke grayscale jika gambar RGB
end

figure(2); clf;

% Tambahkan noise dan filter median untuk setiap jenis noise
for i = 1:5
    switch i
        case 1
            F_noise = imnoise(F, 'gaussian');
            noise_type = 'Gaussian';
        case 2
            F_noise = imnoise(F, 'salt & pepper', 0.05);
            noise_type = 'Salt & Pepper';
        case 3
            F_noise = imnoise(F, 'speckle');
            noise_type = 'Speckle';
        case 4
            F_noise = imnoise(F, 'poisson');
            noise_type = 'Poisson';
        case 5
            noise = uint8(20 * randn(size(F)));
            F_noise = im2uint8(mat2gray(double(F) + double(noise)));
            noise_type = 'Uniform';
    end

    % Terapkan filter median
    [tinggi, lebar] = size(F_noise);
    G = F_noise; % Salin gambar noise untuk hasil filter median
    for baris = 2:tinggi-1
        for kolom = 2:lebar-1
            % Ambil nilai piksel dari jendela 3x3
            data = [F_noise(baris-1, kolom-1) ...
                    F_noise(baris-1, kolom) ...
                    F_noise(baris-1, kolom+1) ...
                    F_noise(baris, kolom-1) ...
                    F_noise(baris, kolom) ...
                    F_noise(baris, kolom+1) ...
                    F_noise(baris+1, kolom-1) ...
                    F_noise(baris+1, kolom) ...
                    F_noise(baris+1, kolom+1)];

            % Urutkan data untuk mendapatkan median
            data_sorted = sort(data);
            G(baris, kolom) = data_sorted(5); % Ambil nilai median (nilai ke-5 setelah diurutkan)
        end
    end

    % Tampilkan gambar noise di kolom pertama
    subplot(2, 5, (i-1)*5 + 1);
    imshow(F_noise);
    title(sprintf('%s Noise', noise_type));

    % Tampilkan hasil filter median di kolom kedua
    subplot(2, 5, (i-1)*5 + 2);
    imshow(G);
    title('Median');
end



