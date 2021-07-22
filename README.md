# Udacity Sensor Fusion Nano Degree - Radar Final Project
## Submitted by Ryan Wicks

### CA-CFAR Processing

For my CA-CFAR processing, I decided to build a kernel that matched the structure discussed in the course notes and perform at 2D convolution with the reshaped signal to generate a noise threshold array. This was done because the implentation was straightforward, and the efficiency is likely to be higher than any naive hand coded implementation. The reason for this likely improvement in efficiency is the conv2 programmers have optimized the code already, and might use different methods to improve the performance depending on the size of the kernel. Furthermore, the use of a continuous kernel rather than a set of conditionals to determine guard or training cells allows to compiler more opportunities to prefetch data and instructions.

To create my CA-CFAR kernel, I created an array of zeroes of size (2*Tr + 2*Gr + 1, 2*Td + 2*Gd + 1) and set each value to 1/Nt, where Nt is the total number of training cells, which is (2*Tr + 2*Gr + 1)*(2*Td + 2*Gd + 1) - (2*Gr + 1)*(2*Gd + 1). I then set the central square of guard cells to 0 to create the guard section. 

Applying the kernel to the data was a simple matter of running the conv2 function with the "same" parameter to produce the same size output as the input data array.

After the appropriate transformations and adding the offset to get the noise threshold, the comparions of the signal to noise was performed as a binary thresholding that was converted back to integers for plotting.

### Choice of Guard and Training Cells

The choice of kernel parameters was made by looking at the size of a signal with a reasonably large velocity. The guard cells were chosen by looking at the width of the peak in the data in both directions, and using half that vale for the guard signal. The training cells were chosen to be twice the size of the guard cells. 

The offset was increased until the noise disappeared.

### Edge Suppression

After generating the binary mask where the signal is greater than the noise threshold, the edges of the data were set to 0. In these regions, the cell was not completely applied, and would tend to be noisier. As a result, I simply set these regions to 0 using array slicing techniques.
