%% README
%in this code we write code for three aims:
%1) pca of one dimension data(without matlab function)
%2) find eignfaces of a dataset consisting different faces using pca 
%3) using k-mean algorithm in finding the color map of a simple photo
%% problem3(normal pca)
clear
clc
load 'C:\Users\Sana\OneDrive\Desktop\semiterm6\neuroscience_karbalaee\homework\HW2_SanaAminnaji_98104722\1d_pcadata.mat' ;
x1 = X(: , 1) ;
x2 = X(: , 2) ;
sz = 25 ;
scatter(x1 , x2 , sz , 'filled') ;
mean1 = [mean(x1) , mean(x2)] ;
x = [x1 , x2] ;
w = zeros(2 ,2) ;
for i = 1:50
    a = x(i , :) - mean1 ;
    w = w + (transpose(a))*a ;
end
w = w./(50-1) ;
[e1 , v1] = eig(w) ;
e11 = e1(: , 1) ;
e21 = e1(: , 2) ;
d1 = sqrt(diag(v1)) ;
hold on ;
quiver(mean(x1) , mean(x2) , e1(1 , 2) , e1(2 , 2) , d1(2) , 'k' , 'LineWidth' , 1) ;
quiver(mean(x1) , mean(x2) , -e1(1 , 2) , -e1(2 , 2) , d1(2) , 'k' , 'LineWidth' , 1) ;
quiver(mean(x1) , mean(x2) , e1(1 , 1) , e1(2 , 1) , d1(1) , 'k' , 'LineWidth' , 1) ;
quiver(mean(x1) , mean(x2) , -e1(1 , 1) , -e1(2 , 1) , d1(1) , 'k' , 'LineWidth' , 1) ;

%% problem3 part2(eign faces)
clear 
clc
load 'C:\Users\Sana\OneDrive\Desktop\semiterm6\neuroscience_karbalaee\homework\HW2_SanaAminnaji_98104722\faces.mat' ;

% %converting the data matrix to pictures
% image2 = zeros(5000 , 32 , 32) ;
% for k = 1:5000
%     image1 = X(k , :) ;
%     for i = 1:32
%         for j = 1:32
%             image2(k , j , i) = image1(j + 32*(i-1)) ;
%         end
%     end
% end
% finalimage = zeros(25*32 , 20*32) ;
% for o = 1: 25
%     for p = 1:20
%         for i = 1:32
%            for j = 1:32
%               finalimage(32*(o-1) + j , 32*(p-1) +i) = image2(p + 20*(o-1) , i , j) ; 
%            end
%         end
%     end
% end
% finalimage1 = mat2gray(finalimage) ;
% finalimage1 = transpose(finalimage1) ;
% imshow(finalimage1) ;

mean1 = zeros(1 , 1024) ;
wight = zeros(1024 , 1024) ;
for i = 1:1024
    mean1(i) = mean(X(: , i)) ;
end

for i = 1:5000
    a = X(i , :) - mean1 ;
    wight = wight + (transpose(a))*a ;
end
wight = wight./(5000 - 1) ;
[vector , value] = eig(wight) ;
d1 = diag(value) ;
d1 = sort(d1) ;

finalimage1 = zeros(32*32 , 32*32) ;
for i = 1:32
    for j = 1:32
        for r = 1:32
            for s = 1:32
              finalimage1(32*(i-1) + r , 32*(j-1) + s) =  vector(s + 32*(r-1) , i + (j-1)*32) ;
            end
        end
    end
end

finalimage3  = zeros(8*32 , 8*32) ;
for i = 1:8
    for j = 1:8
        for r = 1:32
            for s = 1:32
                finalimage3(32*(i-1) + r , 32*(j-1) + s) = vector(32*(r-1) + s , 1024 - 64 + j + (i-1)*8) ;
            end
        end
    end
end
finalimage2 = mat2gray(finalimage3) ;
finalimage2 = transpose(finalimage2) ;
imshow(finalimage2) ;


%% problem3 part3(eignfaces)
clear 
clc
load 'C:\Users\Sana\OneDrive\Desktop\semiterm6\neuroscience_karbalaee\homework\HW2_SanaAminnaji_98104722\faces.mat' ;
X = X(1:20 , :) ;
mean1 = zeros(1 , 1024) ;
wight = zeros(1024 , 1024) ;
for i = 1:1024
    mean1(i) = mean(X(: , i)) ;
end

for i = 1:20
    a = X(i , :) - mean1 ;
    wight = wight + (transpose(a))*a ;
end
wight = wight./(20 - 1) ;
[vector , value] = eig(wight) ;

finalimage1 = zeros(32*32 , 32*32) ;
for i = 1:32
    for j = 1:32
        for r = 1:32
            for s = 1:32
              finalimage1(32*(i-1) + r , 32*(j-1) + s) =  vector(s + 32*(r-1) , i + (j-1)*32) ;
            end
        end
    end
end

finalimage3  = zeros(8*32 , 8*32) ;
for i = 1:8
    for j = 1:8
        for r = 1:32
            for s = 1:32
                finalimage3(32*(i-1) + r , 32*(j-1) + s) = vector(32*(r-1) + s , 1024 - 64 + j + (i-1)*8) ;
            end
        end
    end
end
finalimage2 = mat2gray(finalimage3) ;
finalimage2 = transpose(finalimage2) ;
imshow(finalimage2) ;
%% problem4 (k-mean)
clear 
clc
file = 'C:\Users\Sana\OneDrive\Desktop\semiterm6\neuroscience_karbalaee\homework\HW2_SanaAminnaji_98104722\ponyo.jpeg' ;
[img_first] = imread(file) ;
imshow(img_first) ;
k = 10 ;
new_picture = zeros(225 , 225 , 3) ;
img = double(img_first) ;
condition = 0 ;
o = 0 ;
cluster_center = zeros(k , 3) ;
cluster = randi(k,225) ;
distance1 = [] ;
d = zeros(k , 225 , 225) ;
while (condition == 0)
    sum_center = zeros(k , 3) ;
    size_cluster = zeros(k , 1) ;
    for i = 1:225
        for j = 1:225
            sum_center(cluster(i , j) , 1) = sum_center(cluster(i , j) , 1) + img(i , j , 1) ;
            sum_center(cluster(i , j) , 2) = sum_center(cluster(i , j) , 2) + img(i , j , 2) ;
            sum_center(cluster(i , j) , 3) = sum_center(cluster(i , j) , 3) + img(i , j , 3) ;
            size_cluster(cluster(i , j)) = size_cluster(cluster(i , j)) + 1 ;
        end
    end
    for i = 1:k
        cluster_center(i , 1) = sum_center(i , 1)./size_cluster(i) ;
        cluster_center(i , 2) = sum_center(i , 2)./size_cluster(i) ;
        cluster_center(i , 3) = sum_center(i , 3)./size_cluster(i) ;
    end
    distance = 0 ;
    distance2 = 0 ;
    cluster_t = zeros(225 , 225) ;
    for i = 1:225
        for j = 1:225
            vary = 100000 ;
            for c = 1:k
                distance = sqrt((img(i , j , 1) - cluster_center(c , 1))^2 + (img(i , j , 2) - cluster_center(c , 2))^2 + (img(i , j , 3) - cluster_center(c , 3))^2) ;
                d(c , i , j) = distance ;
                if (c == cluster(i , j))
                    distance2 = distance2 + distance ;
                end
                if (distance < vary)
                    vary = distance ;
                    cluster_t(i , j) = c ;
                end
            end
        end
    end
    if isequal(cluster_t , cluster)
        condition = 1 ;
    end
    cluster = cluster_t ;
    o = o + 1 ;
    distance1 = [distance1 , distance2] ;
end

figure
plot(distance1) ;

for i = 1:225
    for j = 1:225
       new_picture(i , j , 1) = round(cluster_center(cluster(i , j) , 1)) ; 
       new_picture(i , j , 2) = round(cluster_center(cluster(i , j) , 2)) ; 
       new_picture(i , j , 3) = round(cluster_center(cluster(i , j) , 3)) ; 
    end
end

x = new_picture;
xtype = class(x) ;
y = uint8(x) ;
figure
imshow(y) ;
