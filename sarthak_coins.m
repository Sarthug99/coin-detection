close all;
clear all;
clc;
I = imread(input("Enter input image"));  %Asking the user to input the image
[m,n] = size(I);         %Storing the image size

I_gray = rgb2gray(I);     %Converting to grayscale image
subplot 221
imshow(I_gray);             %Plotting the grayscale image
title('Grayscale Image');
%% Gaussian filtering of the input image
kernel = [1,2,1 ; 2,4,2 ; 1,2,1]/16;      %Defining 3*3 Gaussian kernel
zero_padded = zeros(m+2,(n/3 +2));
for i= 1:m
    for j= 1:n/3
        if(i<m+1 && j<n/3+1)                   
            zero_padded(i+1,j+1) = I_gray(i,j);    %zero_padding the grayscale image
        end
    end
end
I_lpf = zeros(m,n/3);
a = 2;
for k = a:m+a-1                  %loop for convolution with Gaussian kernel
    for l= a:n/3+a-1
        for r = -(a-1):a-1
            for s = -(a-1):a-1 
                I_lpf(k-(a-1),l-(a-1)) = I_lpf(k-(a-1),l-(a-1)) + (zero_padded(k+r, l+s)*kernel(r+a,s+a)) ;
            end
        end
    end
end
subplot 222
imshow(uint8(I_lpf));              %plotting the Gaussian Fitered image
title('After Gaussian Filtering');

%% Image Binarization
threshold = graythresh(uint8(I_lpf));     %calculating threshold through graythresh function
I_binary = zeros(m,n/3);
for i = 1:m
    for j = 1:n/3
        if (I_lpf(i,j) > threshold * 256)     %Obtaining the binary image through the threshold 
            I_binary(i,j) = 0;
        else
            I_binary(i,j) = 1;           
        end
    end
end
subplot 223
imshow(I_binary);                 %plotting the binary image
title('Binarized Image'); 

%% Dilation and Erosion
I_dilated = zeros(m,n/3);
for i = 2 : m-1
    for j = 2 : (n/3)-1     %performing dilation with square as the structuring element
        if(I_binary(i,j) == 1 || I_binary(i+1,j) == 1 || I_binary(i+1,j+1) == 1 || I_binary(i,j+1) == 1)
            I_dilated(i,j) = 1;
        end
    end    
end

I_double = zeros(m,n/3);
for i = 2 : m-1
    for j = 2 : (n/3)-1     %performing dilation one more time
        if(I_dilated(i,j) == 1 || I_dilated(i+1,j) == 1 || I_dilated(i+1,j+1) == 1 || I_dilated(i,j+1) == 1)
            I_double(i,j) = 1;
        end
    end    
end

subplot 224 
imshow(I_double);         %plotting the image after dilation
title('After dilation');
SE = [1,1,1 ; 1,1,1 ; 1,1,1];
I_ero = imerode(I_double, SE);   %Performing erosion on the dilated image

%% Algorithm for bwlabel function
B= [0,1,0; 1,1,1; 0,1,0];    %Defining the structuring element for dilation    
A=I_ero;           
p=find(A==1);                % Find a non-zero element's position.
p=p(1);                      %Storing the first non-zero element.
Label=zeros(m, n/3);         %Initialising the label matrix
N=0;

while(~isempty(p))           
    N=N+1;                   %Label for each new component
    p=p(1);
X= zeros(m, n/3);
X(p)=1;                      %Matrix with the first non-zero element of initial matrix

Y=A&imdilate(X,B);            %Performing dilation of previous matrix with structuring element and taking intersection with original matrix.
while(~isequal(X,Y))
    X=Y;                      %Repeating intersection and dilation till dilated matrix and obtained matrix are not the same.
    Y=A&imdilate(X,B);
end

Pos=find(Y==1);               %Storing the positions of non-zero elements in obtained matrix

A(Pos)=0;                     %Converting the pixel values at those positions to zero in the original image.

Label(Pos)=N;                 %Labeling the components

p=find(A==1);                 %Finding the next non-zero element positions in original image and storing in p.

end                           %Continuing this until all non-zero elements in original image have been labelled.

%% Calculating number of coins

number = Label(1,1);            %Different coins will have different labels
for i = 1:m
    for j = 1: n/3
        if(Label(i,j) > number)        %The maximum value in the label matrix will correspond to the number of coins in the image.
            number = Label(i,j);      
        end
    end
end
fprintf('The number of coins in the input image are: %d \n', number);


%% Calculating pixel values

pixel_count = zeros(1,number);        %Initializing the pixel matrix wrt the number of coins obtained.
for k = 1:number
    for i = 1:m
        for j = 1:n/3                  %Pixel positions with a particular label will add to the pixel count of that particular coin
            if(Label(i,j) == k)
                pixel_count(1,k) = pixel_count(1,k) + 1;     
            end
        end 
    end
end
fprintf('Pixel count is:\n');
disp(pixel_count);


%% Displaying the coins alongwith their pixel values

for k = 1:number
    Im=zeros(m, n/3);
    ele=find(Label==k);      %finding pixel positions with the corresponding label.
    Im(ele)=1;
    figure,imshow(Im);          %displaying the image with pixel count and coin number.
    title(['Size in pixel: ',num2str(pixel_count(k)),' ', 'Coin #',num2str(k)])
end