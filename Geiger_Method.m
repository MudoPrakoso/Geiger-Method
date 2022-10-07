clear
clf
% Import data koordinat dari masing - masing stasiun seismik
filename = input('Masukkan file lokasi stasiun seismik ', 's');
import_data = importdata(filename);
% Kode dari masing - masing stasiun seismik
sta_code_x = import_data.textdata;
% Koordinat dari masing - masing stasiun seismik
sta_loc = import_data.data;

% Import data waktu tiba gelombang seismik di masing - masing stasiun
% seismik
filename = input('Masukkan file waktu tiba gelombang seismik ', 's');
import_data = importdata(filename);
% Kode dari masing - masing stasiun seismik yang mencatat gelombang seismik
sta_code_t = import_data.textdata;
% Data waktu tiba masing - masing stasiun seismik yang mencatat gelombang seismik
to = import_data.data;

% Import model kecepatan gelombang seismik
filename = input('Masukkan file model kecepatan gelombang seismik ', 's');
v_mod = importdata(filename);

%--------------------------------------------------------------------------

% Memilih stasiun seismik yang merekam gelombang seismik gempa bumi
nsta_loc = zeros(length(to), size(sta_loc, 2));
for i = 1 : length(to)
    for j = 1 : size(sta_loc, 1)
        if strcmp(sta_code_x(j,1), sta_code_t(i,1))
            for k = 1 : size(sta_loc, 2)
                nsta_loc(i,k) = sta_loc(j,k);
            end
        end
    end
end 

%--------------------------------------------------------------------------

% Memilih salah satu stasiun seismik sebagai stasiun referensi
% Stasiun seismik referensi yang dipilih adalah stasiun seismik SNJI
phi0 = -7.78;
lambda0 = 111.759;
z0 = 1341;

% Mengubah sistem koordinat stasiun seismik menjadi sistem koordinat
% kartesian dengan menggunakan metode Richter
sta_loc = zeros(size(nsta_loc, 1), 3);
for i = 1 : size(nsta_loc, 1)
    phi = 0.5*(nsta_loc(i,1) + phi0);
    sin_phi = sin(phi*pi/180);
    A = 1.8553654 + 0.0062792*(sin_phi)^2 + 0.0000319*(sin_phi)^4;
    B = 1.8428071 + 0.0187098*(sin_phi)^2 + 0.0001583*(sin_phi)^4;
    sta_loc(i,1) = 60*A*(nsta_loc(i,2)-lambda0);
    sta_loc(i,2) = 60*B*(nsta_loc(i,1)-phi0);
    sta_loc(i,3) = -0.001*nsta_loc(i,3);
end

%--------------------------------------------------------------------------

% Konversi waktu tiba gelombang seismik dalam satuan detik
to_s = zeros(size(to,1), 1);
for i = 1 : size(to,1)
    to_s(i,1) = to(i,1)*3600 + to(i,2)*60 + to(i,3);
end

% Menentukan stasiun dengan waktu tiba tersingkat
[~, idx] = min(to_s);

% Menentukan model awal
m = [sta_loc(idx,1); sta_loc(idx,2); 0; 0];

% Menentukan waktu tiba gelombang seismik pada stasiun referensi, SNJI
for i = 1 : size(sta_code_t, 1)
    if strcmp(sta_code_t(i,1), 'SNJI')
        to0 = to_s(i,1);
    end
end

% Menentukan waktu tiba gelombang seismik relatif terhadap waktu tiba
% gelombang seismik pada stasiun referensi, SNJI
to = zeros(size(to_s, 1), 1);
for i = 1 : size(to_s, 1)
    to(i,1) = to_s(i,1) - to0;
end

%--------------------------------------------------------------------------

% Menentukan pada lapisan ke berapa model kedalaman berada, berdasarkan
% model kecepatan yang digunakan
if m(3) <= 20
    j = 1;
elseif m(3) > 20 && m(3) <= 35
    j = 2;
elseif m(3) > 35 && m(3) <= 50
    j = ceil(m(3)/5) - 5;
elseif m(3) > 50 && m(3) <= 800
    j = ceil(m(3)/10);
elseif m(3) > 800 && m(3) <= 1500
    j = ceil(m(3)/100) + 72;
elseif m(3) > 1500 && m(3) <= 2500
    j = ceil(m(3)/500) + 84;
elseif m(3) > 2500 && m(3) <= 2700
    j = 90;
elseif m(3) > 2700 && m(3) <= 2740
    j = 91;
elseif m(3) > 2740 && m(3) <= 2750
    j = 92;
elseif m(3) > 2750 && m(3) <= 2800
    j = 93;
elseif m(3) > 2800 && m(3) <= 2850
    j = 94;
elseif m(3) > 2850 && m(3) <= 2889
    j = 95;
elseif m(3) > 2889 && m(3) <= 2900
    j = 96;
elseif m(3) > 2900 && m(3) <= 3500
    j = ceil(m(3)/100) + 67;
elseif m(3) > 3500 && m(3) <= 4500
    j = ceil(m(3)/500) + 95;
elseif m(3) > 4500 && m(3) <= 5153.9
    j = 105;
elseif m(3) > 5153.9 && m(3) <= 5500
    j = 106;
elseif m(3) > 5500 && m(3) <= 6000
    j = 107;
elseif m(3) > 6000 && m(3) <= 6371
    j = 108;
end

% Mendefinisikan variabel yang akan digunakan
tp = zeros(size(sta_loc, 1), 1);
deltas = zeros(size(sta_loc, 1), 1);
ethab = zeros(size(sta_loc, 1), 2);
sin_phib = zeros(size(sta_loc, 1), j);
sin_phibj = zeros(size(sta_loc, 1), j, 2);
phibj = zeros(size(sta_loc, 1), j, 2);
ethabj = zeros(size(sta_loc, 1), j, 2);
deltab = zeros(size(sta_loc, 1), 2);
etha = zeros(size(sta_loc, 1), 1);
sin_phi = zeros(size(sta_loc, 1), 1);
sin_phij = zeros(size(sta_loc, 1), j);
phij = zeros(size(sta_loc, 1), j);
ethaj = zeros(size(sta_loc, 1), j);
delta = zeros(size(sta_loc, 1), 1);
ttj = zeros(size(sta_loc, 1), j);

for i = 1 : size(sta_loc, 1)
    
    % Menghitung waktu tiba gelombang seismik apabila model berada di
    % lapisan pertama
    if j == 1;
        tp(i,1) = m(4) + sqrt((m(1)-sta_loc(i,1))^2 + (m(2)-sta_loc(i,2))^2 + (m(3)-sta_loc(i,3))^2)/v_mod(j,2);
    
    % Menghitung waktu tiba gelombang seismik apabila model berada di 
    % lapisan kedua atau seterusnya    
    else
        
        % Menghitung waktu tiba apabila stasiunnya merupakan stasiun
        % yang pertama kali menerima gelombang gempa
        if i == idx
            for k = 1 : j-1
                if k == 1
                    ttj(i,k) = (v_mod(k,1)-sta_loc(i,3))/v_mod(k,1);
                else
                    ttj(i,k) = v_mod(k,1)/v_mod(k,2);
                end
            end
            tp(idx,1) = m(4) + ((m(3)-v_mod(j,1))/v_mod(j,2)) + sum(ttj(i,:));
        
        % Menghitung waktu tiba apabila stasiunnya merupakan bukan stasiun 
        % yang pertama kali merima gelombang seismik    
        else
            % Menentukan batas bawah dan batas atas dari etha dan delta
            deltas(i,1) = sqrt((m(1)-sta_loc(i,1))^2 + (m(2)-sta_loc(i,2))^2);
            ethab(i,1) = deltas(i,1)*(m(3)-v_mod(j-1,1))/(m(3)-sta_loc(i,3));
            ethab(i,2) = deltas(i,1);
            for k = 1 : 2
                sin_phib(i,k) = ethab(i,k)/sqrt(ethab(i,k)^2 + (m(3)-v_mod(j-1,1))^2);
                for l = 1 : j-1
                    sin_phibj(i,l,k) = v_mod(l,2)/v_mod(j,2)*sin_phib(i,k);
                    phibj(i,l,k) = asin(sin_phibj(i,l,k));
                    if l == 1;
                        ethabj(i,l,k) = (v_mod(l,1)-sta_loc(i,3))*tan(phibj(i,l,k));
                    else
                        ethabj(i,l,k) = (v_mod(l,1)-v_mod(l-1,1))*tan(phibj(i,l,k));
                    end
                end
                deltab(i,k) = ethab(i,k) + sum(ethabj(i,:,k));
            end
        
            % Proses iterasi dalam menentukan jarak sumber gempa bumi dengan
            % masing - masing stasiun seismik
            etha(i,1) = (deltas(i,1)-deltab(i,1))/(deltab(i,2)-deltas(i,1))*(ethab(i,2)-ethab(i,1))+ethab(i,1);
            sin_phi(i,1) = etha(i,1)/sqrt(etha(i,1)^2 + (m(3) - v_mod(j-1,1))^2);
            for k = 1 : j-1
                sin_phij(i,k) = v_mod(k,2)/v_mod(j,2)*sin_phi(i,1);
                phij(i,k) = asin(sin_phij(i,k));
                if k == 1
                    ethaj(i,k) = (v_mod(k,1)-sta_loc(i,3))*tan(phij(i,k));
                else
                    ethaj(i,k) = (v_mod(k,1)-v_mod(k-1,1))*tan(phij(i,k));
                end
            end
            delta(i,1) = etha(i,1) + sum(ethaj(i,:));
        
            % Apabila jarak antara delta dengan delta sumber kurang dari 10
            % meter iterasi akan dihentikan dan delta dianggap sebagai delta
            % sumber
            while abs(deltas(i,1)-delta(i,1)) > 0.01
            
                if delta(i,1) < deltas(i,1)
                    ethab(i,1) = etha(i,1);
                    deltab(i,1) = delta(i,1);
                elseif delta(i,1) > deltas(i,1)
                    ethab(i,2) = etha(i,1);
                    deltab(i,2) = delta(i,1);
                end
                etha(i,1) = ethab(i,1) + 0.5*(ethab(i,2)-ethab(i,1));
                sin_phi(i,1) = etha(i,1)/sqrt(etha(i,1)^2 + (m(3) - v_mod(j-1,1))^2);
                for k = 1 : j-1
                    sin_phij(i,k) = v_mod(k,2)/v_mod(j,2)*sin_phi(i,1);
                    phij(i,k) = asin(sin_phij(i,k));
                    if k == 1
                        ethaj(i,k) = (v_mod(k,1)-sta_loc(i,3))*tan(phij(i,k));
                    else
                        ethaj(i,k) = (v_mod(k,1)-v_mod(k-1,1))*tan(phij(i,k));
                    end
                end
                delta(i,1) = etha(i,1) + sum(ethaj(i,:));
            end
        
            % Menghitung waktu tiba gelombang seismik untuk j lapisan
            for k = 1 : j-1
                ttj(i,k) = ethaj(i,k)/(v_mod(k,2)*sin_phij(i,k));
            end
            tp(i,1) = m(4) + (sqrt(etha(i,1)^2 + (m(3)-v_mod(j-1,1))^2)/v_mod(j,2)) + sum(ttj(i,:));
        end
    end
end

%--------------------------------------------------------------------------

% Melakukan proses inversi dengan menggunakan metode least square
res = to-tp;
sres = sum(res.^2);
rms = sqrt(sres/length(res));
lambda = input('Masukkan nilai faktor redaman ');
lambdaup = input('Masukkan nilai faktor penambahan dari faktor redaman ');
lambdadown = input('Masukkan nilai faktor pengurangan dari faktor redaman ');
alpha = zeros(length(res), 1);
ethax = zeros(length(res), 1);
ethay = zeros(length(res), 1);
rj = zeros(length(res), 1);
G = zeros(length(res), 3);
if j == 1
    for i = 1 : length(res)
        rj(i,1) = sqrt((m(1)-sta_loc(i,1))^2 + (m(2)-sta_loc(i,2))^2 + (m(3)-sta_loc(i,3))^2);
        G(i,1) = (m(1)-sta_loc(i,1))/(v_mod(1,2)*rj(i,1));
        G(i,2) = (m(2)-sta_loc(i,2))/(v_mod(1,2)*rj(i,1));
        G(i,3) = (m(3)-sta_loc(i,3))/(v_mod(1,2)*rj(i,1));
        G(i,4) = 1;
    end
else
    for i = 1 : length(res)
        if i == idx;
            G(i,1) = 0;
            G(i,2) = 0;
            G(i,3) = 1/v_mod(j,2);
            G(i,4) = 1;
        else
            alpha(i,1) = atan((sta_loc(i,2)-m(2))/(sta_loc(i,1)-m(1)));
            ethax(i,1) = etha(i,1)*cos(alpha(i,1));
            ethay(i,1) = etha(i,1)*sin(alpha(i,1));
            rj(i,1) = sqrt(etha(i,1)^2 + (m(3)-v_mod(j-1,1))^2);
            G(i,1) = -ethax(i,1)/(v_mod(j,2)*rj(i,1));
            G(i,2) = -ethay(i,1)/(v_mod(j,2)*rj(i,1));
            G(i,3) = (m(3)-v_mod(j-1,1))/(v_mod(j,2)*rj(i,1));
            G(i,4) = 1;
        end
    end
end
deltam = (G.'*G + lambda.*eye(length(4)))\(G.'*res);
% Perubahan posisi terhadap model awal atau model sebelumnya
deltax = sqrt(sum(deltam(1:3).^2));
% Memperbaruhi model awal atau model sebelumnya
m = m + deltam;

%--------------------------------------------------------------------------

for a = 1:100;
    
    %----------------------------------------------------------------------
    
    % Menentukan pada lapisan ke berapa model kedalaman berada, berdasarkan
    % model kecepatan yang digunakan
    if m(3) <= 20
        j = 1;
    elseif m(3) > 20 && m(3) <= 35
        j = 2;
    elseif m(3) > 35 && m(3) <= 50
        j = ceil(m(3)/5) - 5;
    elseif m(3) > 50 && m(3) <= 800
        j = ceil(m(3)/10);
    elseif m(3) > 800 && m(3) <= 1500
        j = ceil(m(3)/100) + 72;
    elseif m(3) > 1500 && m(3) <= 2500
        j = ceil(m(3)/500) + 84;
    elseif m(3) > 2500 && m(3) <= 2700
        j = 90;
    elseif m(3) > 2700 && m(3) <= 2740
        j = 91;
    elseif m(3) > 2740 && m(3) <= 2750
        j = 92;
    elseif m(3) > 2750 && m(3) <= 2800
        j = 93;
    elseif m(3) > 2800 && m(3) <= 2850
        j = 94;
    elseif m(3) > 2850 && m(3) <= 2889
        j = 95;
    elseif m(3) > 2889 && m(3) <= 2900
        j = 96;
    elseif m(3) > 2900 && m(3) <= 3500
        j = ceil(m(3)/100) + 67;
    elseif m(3) > 3500 && m(3) <= 4500
        j = ceil(m(3)/500) + 95;
    elseif m(3) > 4500 && m(3) <= 5153.9
        j = 105;
    elseif m(3) > 5153.9 && m(3) <= 5500
        j = 106;
    elseif m(3) > 5500 && m(3) <= 6000
        j = 107;
    elseif m(3) > 6000 && m(3) <= 6371
        j = 108;
    end

    % Mendefinisikan variabel yang akan digunakan
    tp = zeros(size(sta_loc, 1), 1);
    deltas = zeros(size(sta_loc, 1), 1);
    ethab = zeros(size(sta_loc, 1), 2);
    sin_phib = zeros(size(sta_loc, 1), j);
    sin_phibj = zeros(size(sta_loc, 1), j, 2);
    phibj = zeros(size(sta_loc, 1), j, 2);
    ethabj = zeros(size(sta_loc, 1), j, 2);
    deltab = zeros(size(sta_loc, 1), 2);
    etha = zeros(size(sta_loc, 1), 1);
    sin_phi = zeros(size(sta_loc, 1), 1);
    sin_phij = zeros(size(sta_loc, 1), j);
    phij = zeros(size(sta_loc, 1), j);
    ethaj = zeros(size(sta_loc, 1), j);
    delta = zeros(size(sta_loc, 1), 1);
    ttj = zeros(size(sta_loc, 1), j);

    for i = 1 : size(sta_loc, 1)

        % Menghitung waktu tiba gelombang seismik untuk lapisan pertama
        if j == 1;
            tp(i,1) = m(4) + sqrt((m(1)-sta_loc(i,1))^2 + (m(2)-sta_loc(i,2))^2 + (m(3)-sta_loc(i,3))^2)/v_mod(j,2);

        % Menentukan batas bawah dan batas atas dari etha dan delta
        else
            deltas(i,1) = sqrt((m(1)-sta_loc(i,1))^2 + (m(2)-sta_loc(i,2))^2);
            ethab(i,1) = deltas(i,1)*(m(3)-v_mod(j-1,1))/(m(3)-sta_loc(i,3));
            ethab(i,2) = deltas(i,1);
            for k = 1 : 2
                sin_phib(i,k) = ethab(i,k)/sqrt(ethab(i,k)^2 + (m(3)-v_mod(j-1,1))^2);
                for l = 1 : j-1
                    sin_phibj(i,l,k) = v_mod(l,2)/v_mod(j,2)*sin_phib(i,k);
                    phibj(i,l,k) = asin(sin_phibj(i,l,k));
                    if l == 1;
                        ethabj(i,l,k) = (v_mod(l,1)-sta_loc(i,3))*tan(phibj(i,l,k));
                    else
                        ethabj(i,l,k) = (v_mod(l,1)-v_mod(l-1,1))*tan(phibj(i,l,k));
                    end
                end
                deltab(i,k) = ethab(i,k) + sum(ethabj(i,:,k));
            end

            % Proses iterasi dalam menentukan jarak sumber gempa bumi dengan
            % masing - masing stasiun seismik
            etha(i,1) = (deltas(i,1)-deltab(i,1))/(deltab(i,2)-deltas(i,1))*(ethab(i,2)-ethab(i,1))+ethab(i,1);
            sin_phi(i,1) = etha(i,1)/sqrt(etha(i,1)^2 + (m(3) - v_mod(j-1,1))^2);
            for k = 1 : j-1
                sin_phij(i,k) = v_mod(k,2)/v_mod(j,2)*sin_phi(i,1);
                phij(i,k) = asin(sin_phij(i,k));
                if k == 1
                    ethaj(i,k) = (v_mod(k,1)-sta_loc(i,3))*tan(phij(i,k));
                else
                    ethaj(i,k) = (v_mod(k,1)-v_mod(k-1,1))*tan(phij(i,k));
                end
            end
            delta(i,1) = etha(i,1) + sum(ethaj(i,:));

            % Apabila jarak antara delta dengan delta sumber kurang dari 10
            % meter iterasi akan dihentikan dan delta dianggap sebagai delta
            % sumber
            while abs(deltas(i,1)-delta(i,1)) > 0.01

                if delta(i,1) < deltas(i,1)
                    ethab(i,1) = etha(i,1);
                    deltab(i,1) = delta(i,1);
                elseif delta(i,1) > deltas(i,1)
                    ethab(i,2) = etha(i,1);
                    deltab(i,2) = delta(i,1);
                end
                etha(i,1) = ethab(i,1) + 0.5*(ethab(i,2)-ethab(i,1));
                sin_phi(i,1) = etha(i,1)/sqrt(etha(i,1)^2 + (m(3) - v_mod(j-1,1))^2);
                for k = 1 : j-1
                    sin_phij(i,k) = v_mod(k,2)/v_mod(j,2)*sin_phi(i,1);
                    phij(i,k) = asin(sin_phij(i,k));
                    if k == 1
                        ethaj(i,k) = (v_mod(k,1)-sta_loc(i,3))*tan(phij(i,k));
                    else
                        ethaj(i,k) = (v_mod(k,1)-v_mod(k-1,1))*tan(phij(i,k));
                    end
                end
                delta(i,1) = etha(i,1) + sum(ethaj(i,:));
            end

            % Menghitung waktu tiba gelombang seismik untuk j lapisan
            for k = 1 : j-1
                ttj(i,k) = ethaj(i,k)/(v_mod(k,2)*sin_phij(i,k));
            end
            tp(i,1) = m(4) + (sqrt(etha(i,1)^2 + (m(3)-v_mod(j-1,1))^2)/v_mod(j,2)) + sum(ttj(i,:));
        end
    end

    %----------------------------------------------------------------------
    
    % Melakukan proses inversi dengan menggunakan metode least square
    res = to-tp;

    nsres = sum(res.^2);
    if nsres - sres < 0
        lambda = lambda / lambdadown;
    else
        lambda = lambda * lambdaup;
    end
    sres = nsres;
    rms = sqrt(sres/length(res));

    alpha = zeros(length(res), 1);
    ethax = zeros(length(res), 1);
    ethay = zeros(length(res), 1);
    rj = zeros(length(res), 1);
    G = zeros(length(res), 4);
    if j == 1
        for i = 1 : length(res)
            rj(i,1) = sqrt((m(1)-sta_loc(i,1))^2 + (m(2)-sta_loc(i,2))^2 + (m(3)-sta_loc(i,3))^2);
            G(i,1) = (m(1)-sta_loc(i,1))/(v_mod(1,2)*rj(i,1));
            G(i,2) = (m(2)-sta_loc(i,2))/(v_mod(1,2)*rj(i,1));
            G(i,3) = (m(3)-sta_loc(i,3))/(v_mod(1,2)*rj(i,1));
            G(i,4) = 1;
        end
    else
        for i = 1 : length(res)
            alpha(i,1) = atan((sta_loc(i,2)-m(2))/(sta_loc(i,1)-m(1)));
            ethax(i,1) = etha(i,1)*cos(alpha(i,1));
            ethay(i,1) = etha(i,1)*sin(alpha(i,1));
            rj(i,1) = sqrt(etha(i,1)^2 + (m(3)-v_mod(j-1,1))^2);
            G(i,1) = -ethax(i,1)/(v_mod(j,2)*rj(i,1));
            G(i,2) = -ethay(i,1)/(v_mod(j,2)*rj(i,1));
            G(i,3) = (m(3)-v_mod(j-1,1))/(v_mod(j,2)*rj(i,1));
            G(i,4) = 1;
        end
    end
    deltam = (G.'*G + lambda.*eye(4))\(G.'*res);
    deltax = sqrt(sum(deltam(1:3).^2));
    m = m + deltam;
    
    %----------------------------------------------------------------------
    
    % Menentukan model variasi dan kovariasi
    s2 = sum(res.^2)/(length(res)-length(size(G, 2)));
    sigma2 = (inv(G.'*G + lambda.*eye(4))*G.'*G*inv(G.'*G + lambda.*eye(4))).*s2;
    sigma = zeros(size(sigma2, 1), 1);
    for i = 1 : size(sigma2, 1)
        for j = 1 : size(sigma2, 2)
            if i == j
                sigma(i,1) = sqrt(sigma2(i,j));
            end
        end
    end
    
    %----------------------------------------------------------------------
    % Membuat rangkuman hasil inversi dari setiap iterasi
    resume(a,1) = a;
    resume(a,2) = rms;
    resume(a,3) = deltax;
    resume(a,4:7) = m(1:4,1);
    resume(a,8:11) = sigma(1:4,1);
end

%--------------------------------------------------------------------------
% Mentranformasi kembali koordinat hasil perhitungan dalam sistem bola
[~, idx] = min(resume(:,2));
y = @(phi)(60*(1.8428071 + 0.0187098*(sind(0.5*(phi+phi0)))^2 + 0.0001583*(sind(0.5*(phi+phi0)))^4)*(phi-phi0) - resume(idx,5));
phia = 0;
phib = -70;
ya = y(phia);
dphi = abs(phia-phib);
while dphi > 0.0001
    phic = (phia+phib)/2;
    yc = y(phic);
    if ya*yc < 0
        phib = phic;
    elseif ya*yc > 0
        phia = phic;
        ya = y(phic);
    end
    dphi = abs(phia-phib);
end
phi = (phia+phib)/2;

%--------------------------------------------------------------------------
% Hasil posisi pusat gempa bumi dan waktu relatif terjadinya gempa bumi
m(2,1) = phi;
m(1,1) = resume(idx,4)/(60*(1.8553654 + 0.0062792*(sind(0.5*(phi+phi0)))^2 + 0.0000319*(sind(0.5*(phi+phi0)))^4)) + lambda0;
m(3:4,1) = resume(idx,6:7);
dm(1:4,1) = resume(idx,8:11);
rms = resume(idx,2);

%--------------------------------------------------------------------------
% Menentukan waktu terjadinya gempa bumi
t0 = to0 + m(4,1);
if t0 >= 3600
hh = floor(t0/3600);
minute = t0 - hh*3600;
    if minute >= 60
        mm = floor(minute/60);
        ss = minute - mm*60;
    else
        mm = 0;
        ss = minute;
    end
else
    hh = 0;
    minute = t0;
    if minute >= 60
        mm = floor(minute/60);
        ss = minute - mm*60;
    else
        mm = 0;
        ss = minute;
    end
end

%--------------------------------------------------------------------------
% Menampilkan grafik rms setiap iterasi
stem(resume(:,1),resume(:,2))
title('Perubahan nilai Root Mean Square Error terhadap iterasi')
xlabel('Iterasi')
ylabel('Root Mean Square Error(s)')
hold on
stem(resume(idx,1),resume(idx,2),'r')

%--------------------------------------------------------------------------
% Menampilkan hasil posisi pusat gempa bumi dan waktu terjadinya gempa bumi
fprintf('Gempa bumi terjadi pukul %d:%d:%.2f +- %.2f s\n', hh, mm, ss, dm(4,1))
disp('Posisi pusat gempa bumi:')
fprintf('Bujur    : %.3f derajat +- %.3f kilometer\n', m(1,1), dm(1,1))
fprintf('Lintang  : %.3f derajat +- %.3f kilometer\n', m(2,1), dm(2,1))
fprintf('Kedalaman: %.3f kilometer +- %.3f kilometer\n', m(3,1), dm(3,1))
fprintf('Root Mean Square Error %.3f s\n', rms)