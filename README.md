LibClustering
=============

This is a C++ port of the MATLAB code by L. Zelnik-Manor and P. Perona http://www.vision.caltech.edu/lihi/Demos/SelfTuningClustering.html. The self-tuning spectral clustering automatically finds the ideal number of clusters of an affinity matrix.

Build
-----

You will need the Eigen2 library (libeigen2-devel) and CMake (only tested under Linux). 

	$ mkdir build
	$ cd build
	$ cmake ../
	$ make


Overview
--------

Spectral clustering derives its name from spectral analysis of a graph, which is how the data is represented. Each object to be clustered can initially be represented as an n-dimensional numeric vector, but there must also be some method for performing a comparison between each object and expressing this comparison as a scalar.

This n by n comparison of all objects with all others forms the affinity matrix, which can be intuitively thought of as a rough representation of an underlying undirected, weighted, and fully-connected graph whose edges express the relative relationships, or affinities, between each pair of objects in the original data.

This affinity matrix is then decomposed into its eigenvectors, forming a matrix where each column is an eigenvector.

Using K-means, each row of this matrix is clustered using a standard K-means algorithm.

Using self-tuning clustering, this matrix is Givens rotated to try to get each row with only one non-zero element. The algorithm then removes the smallest eigenvector and performs the same rotation. The number of eigenvectors where first rotation exists with only one non-zero element is the ideal number of clusters. For details on the algorithm, see http://www.vision.caltech.edu/lihi/Publications/SelfTuningClustering.pdf.


Example
-------

Follow the instructions to build the library, then:

	$ cd ../example
	$ cmake .
	$ make

Implementation Details
----------------------

Due to the algorithm being ported from MATLAB, the evrot.cpp implementation uses a slightly odd style where objects are allocated on the heap, but passed around as references rather than pointers. This is so that the operator overloads provided by the Eigen2 library could be used naturally, but the structure of the program kept the same as in MATLAB. 
