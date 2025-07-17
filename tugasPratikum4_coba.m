pkg load image
img_distorted = imread('E:\pengolahan citra/bunga-hibiscus-atau-kembang-sepatu_169.jpeg');
imshow(img_distorted);
title('Citra Terdistorsi');
[rows, cols, ~] = size(img_distorted);
center_x = cols / 2;
center_y = rows / 2;

k = -0.0000015; % Contoh koefisien distorsi barel (negatif untuk "membuka" distorsi)
                % Nilai ini sangat tergantung pada tingkat distorsi Anda.
img_corrected = uint8(zeros(rows, cols, 3)); % Citra hasil koreksi

for y_corrected = 1:rows
    for x_corrected = 1:cols
        % Hitung koordinat relatif terhadap pusat di citra hasil
        dx = x_corrected - center_x;
        dy = y_corrected - center_y;

        % Hitung jarak radial dari pusat di citra hasil
        r_corrected_sq = dx^2 + dy^2;
        r_corrected = sqrt(r_corrected_sq);

        % Hitung faktor skala distorsi (transformasi balik)
        % Untuk distorsi barel, kita "menekan" kembali ke dalam
        % Persamaan ini adalah penyederhanaan dari model distorsi radial
        factor = 1 + k * r_corrected_sq;

        % Hitung koordinat di citra terdistorsi
        x_distorted = round(center_x + dx * factor);
        y_distorted = round(center_y + dy * factor);

        % Periksa batas citra
        if (x_distorted >= 1 && x_distorted <= cols && y_distorted >= 1 && y_distorted <= rows)
            img_corrected(y_corrected, x_corrected, :) = img_distorted(y_distorted, x_distorted, :);
        end
    end
end

figure;
imshow(img_corrected);
title('Citra Ter skeling');
