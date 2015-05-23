LibClustering
=============

This is a C++ port of the MATLAB code by L. Zelnik-Manor and P. Perona http://www.vision.caltech.edu/lihi/Demos/SelfTuningClustering.html. The self-tuning spectral clustering automatically finds the ideal number of clusters of an affinity matrix.

Build
-----

You will need the Eigen2 library (libeigen2-devel) and CMake (only tested under Linux). 

	$ mkdir release
	$ cd release
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

	#include <vector>
	#include <algorithm>
	#include <iterator>
	#include <cmath>
	#include <Eigen/Core>
	#include <clustering/SpectralClustering.h>
	
	int main() {
		std::vector<int> items = {1,2,3,4,5,6,7,8,9,10};
	
		// generate similarity matrix
		unsigned int size = items.size();
		Eigen::MatrixXd m = Eigen::MatrixXd::Zero(size,size);
	
		for (unsigned int i=0; i < size; i++) {
			for (unsigned int j=0; j < size; j++) {
				// generate similarity
				int d = items[i] - items[j];
				int similarity = exp(-d*d / 100);
				m(i,j) = similarity;
				m(j,i) = similarity;
			}
		}
	
		// the number of eigenvectors to consider. This should be near (but greater) than the number of clusters you expect. Fewer dimensions will speed up the clustering
		int numDims = size;
	
		// do eigenvalue decomposition
		SpectralClustering* c = new SpectralClustering(m, numDims);
	
		// whether to use auto-tuning spectral clustering or kmeans spectral clustering
		bool autotune = true;
	
		std::vector<std::vector<int> > clusters;
		if (autotune) {
			// auto-tuning clustering
			clusters = c->clusterRotate();
		} else {
			// how many clusters you want
			int numClusters = 5;
			clusters = c->clusterKmeans(numClusters);
		}
	
		// output clustered items
		// items are ordered according to distance from cluster centre
		for (unsigned int i=0; i < clusters.size(); i++) {
			std::cout << "Cluster " << i << ": " << "Item ";
			std::copy(clusters[i].begin(), clusters[i].end(), std::ostream_iterator<int>(std::cout, ", "));
			std::cout << std::endl;
		}
	}


Implementation Details
----------------------

Due to the algorithm being ported from MATLAB, the evrot.cpp implementation uses a slightly odd style where objects are allocated on the heap, but passed around as references rather than pointers. This is so that the operator overloads provided by the Eigen2 library could be used naturally, but the structure of the program kept the same as in MATLAB. 
