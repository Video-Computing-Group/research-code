Cell segmentation code to use for cell tracking

##Instructions

The images should have 1, 2, 3, names, and there can’t be gaps between names. All the images should be tif format.

In run.m 
1. Change line 7 to the directory containing the images. Note the path should end with ‘/’. Example: dir = ‘./images/’
2. Based on the images folder the images are sorted alphabatically.

3. Then execute the run.m code.

Plotting segmentation results:
1. Use code plotSegment.m as follows
    For image number 7th alphabatically ordered image in the folder specified:
	plotSegment(7)
