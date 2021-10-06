# coin-detection
Wrote a MATLAB code for detecting the number of coins along with their pixel size in the given image.

The approach followed in the code is as follows:
1) Convert image to grayscale
2) Convolution with low pass filter
3) Image Binarization
4) Dilation and Erosion
5) bwlabel algorithm
6) Counting the number of coins along with pixel sizes.
