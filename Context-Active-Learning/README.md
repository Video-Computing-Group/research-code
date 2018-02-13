# Active Learning on Data with Contextual Relationships

## Overview
This package is an implementation of the of the paper [Non-Uniform Subset Selection for Active Learning in Structured Data](http://openaccess.thecvf.com/content_cvpr_2017/papers/Paul_Non-Uniform_Subset_Selection_CVPR_2017_paper.pdf), by [Sujoy Paul](www.ee.ucr.edu/~supaul/
), [Jawad Bappy](www.ee.ucr.edu/~mbappy/
) and [Amit K Roy-Chowdhury](http://www.ee.ucr.edu/~amitrc/) and published at [CVPR 2017](http://cvpr2017.thecvf.com/).

## Dependencies
This package uses or depends on the the following package:
1. [Submodular Function Optimization](https://www.mathworks.com/matlabcentral/fileexchange/20504-submodular-function-optimization) (Used after some modification)
2. [LibSVM](https://www.csie.ntu.edu.tw/~cjlin/libsvm/) (Used unmodified)
3. [Undirected Graphical Model](https://www.cs.ubc.ca/~schmidtm/Software/UGM.html) (Used unmodified)

All these packages are included in this package.

## Data
This package contains a sample dataset cora.mat which is a processed version of the [CORA](https://linqs.soe.ucsc.edu/node/236) dataset

## Running
Please run the main.m code in order to get results on the cora.mat dataset.

## Citation
Please cite the following work if you use this package.
```javascript
@inproceedings{paul2017non,
  title={Non-uniform subset selection for active learning in structured data},
  author={Paul, Sujoy and Bappy, Jawadul H and Roy-Chowdhury, Amit K},
  booktitle={2017 IEEE Conference on Computer Vision and Pattern Recognition (CVPR)},
  pages={830--839},
  year={2017},
  organization={IEEE}
}
```

## Contact 
Please contact the first author of the associated paper - Sujoy Paul (supaul@ece.ucr.edu) for any further queries.

