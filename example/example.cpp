#include <vector>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <Eigen/Core>
#include "SpectralClustering.h"

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