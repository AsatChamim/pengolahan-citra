##Tugas Praktikum 2%
##
##Memuat package yang dibutuhkan histogram
 pkg load image;
##
## 1. Menampilkan Histogram
img = imread('wallhaven-jxprop.jpg');
img = im2uint8(img); %gambar dalam format uint8
figure;
imhist(img);

R =img(:, :, 1); % Pisahkan Channel RGB
G =img(:, :, 2);
B =img(:, :, 3);

figure;
subplot(3,1,1);
imhist(R);
title('Histogram Merah');
xlim([0 255]);

subplot(3,1,2);
imhist(G);
title('Histogram Hijau');
xlim([0 255]);

subplot(3,1,3);
imhist(B);
title('Histogram Biru');
xlim([0 255]);

## 2. Meningkatkan Kecerahan
figure;
subplot(2,1,1);
imgGray = rgb2gray(img);
imshow(imgGray);
imgCerah = imgGray + 50;
imshow(imgCerah);

subplot(2,1,2);
title('Histogram Setelah Dicerahkan');
imhist(imgCerah);
##3.Meregangkan Kontras
figure;
subplot(2,1,1);
imshow(imgGray);
title('Citra Asli');


imgKontras = imadjust(imgGray, stretchlim(imgGray), []);
subplot(2,1,2);
imshow(imgKontras);
title('Citra Setelah Kontras Stretching');

subplot(2,1,3);
imhist(imgKontras);
title('Histogram Kontras Stretching');


##4.Kombinasi kecerahan dan kontras
 Kombinasi kecerahan dan kontras
figure;
subplot(2,1,1);
##imgKombinasi = imadjust(imgGray, [], [], 1.2) + 30; % kontras 1.2 dan brightness +30
imshow(imgKombinasi);
title('Kombinasi Kecerahan & Kontras');


subplot(2,1,2);
imhist(imgKombinasi);
title('Histogram Kombinasi');

#5.Membalik citra
figure;
subplot(2,1,1);
imgbalik = 255 - imgGray;
imshow(imgbalik);
title('Gambar Dibalik');

subplot(2,1,2);
imhist(imgbalik);
title('Histogram Dibalik');


#6.Pemetaan nonlinear
imgLog = uint8(100 * log(1 + double(imgGray)));
figure;
subplot(2,1,1);
imshow(imgLog);
title('Gambar Hasil Pemetaan Log');
subplot(2,1,2);
imhist(imgLog);
title('Histogram Hasil Log Transformasi');




#7.Pemotongan aras keabuan
threshold = 100; % nilai threshold
imgThreshold = imgGray > threshold;
figure;
imshow(imgThreshold);
title('Hasil Thresholding');



#8.Ekualisasi Histogram
imgEqual = histeq(imgGray);
figure;
subplot(2,1,1);
imshow(imgEqual);
title('Gambar Setelah Ekualisasi Histogram');

subplot(2,1,2);
imhist(imgEqual);
title('Histogram Setelah Ekualisasi');



