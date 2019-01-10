Cell segmentation code to use for cell tracking

##Instructions

The images should have 1, 2, 3, names, and there can’t be gaps between names. All the images should be tif format.

In analyticsGathering.m 
1. Change line 7 to the directory containing the images. Note the path should end with ‘/’. Example: dir = ‘./images/’
2. Based on number of images in your images folder change line 11. If you have 10 image sin your images folder line 11 will be:
		imageInd = [1 10];

3. Then run the  analyticsGathering.m code.

Plotting segmentation results:
1. Use code plotSegment.m as follows
    For image number 7.tif:
	plotSegment(7)
