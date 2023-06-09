/*
 * ClusterRotate.cpp
 *
 *  Created on: 04-Mar-2009
 *      Author: sbutler
 */

#include "ClusterRotate.h"

#include <map>
#include <Eigen/Dense>

ClusterRotate::ClusterRotate(int method):
	mMethod(method),
	mMaxQuality(0)
{

}

std::vector<std::vector<int> > ClusterRotate::cluster(Eigen::MatrixXd& X) {
	mMaxQuality = 0;
	std::vector<std::vector<int> > clusters;
	Eigen::MatrixXd vecRot;
	Eigen::MatrixXd vecIn = X.block(0,0,X.rows(),2);
	Evrot* e = NULL;
	for (int g=2; g <= X.cols(); g++) {
		// make it incremental (used already aligned vectors)
		if( g > 2 ) {
			vecIn.resize(X.rows(),g);
			vecIn.block(0,0,vecIn.rows(),g-1) = e->getRotatedEigenVectors();
			vecIn.block(0,g-1,X.rows(),1) = X.block(0,g-1,X.rows(),1);
			delete e;
		}
		//perform the rotation for the current number of dimensions
		e = new Evrot(vecIn, mMethod);

		//save max quality
		if (e->getQuality() > mMaxQuality) {
			mMaxQuality = e->getQuality();
		}
		//save cluster data for max cluster or if we're near the max cluster (so prefer more clusters)
		if ((e->getQuality() > mMaxQuality) || (mMaxQuality - e->getQuality() <= 0.001)) {
			clusters = e->getClusters();
			vecRot = e->getRotatedEigenVectors();
		}
	}

	Eigen::MatrixXd clusterCentres = Eigen::MatrixXd::Zero(clusters.size(),vecRot.cols());
	for (unsigned int i=0; i < clusters.size(); i++) {
		for (unsigned int j=0; j < clusters[i].size(); j++) {
			//sum points within cluster
			clusterCentres.row(i) += vecRot.row(clusters[i][j]);
		}
	}
	for (unsigned int i=0; i < clusters.size(); i++) {
		//find average point within cluster
		clusterCentres.row(i) = clusterCentres.row(i) / clusters[i].size();
	}

	//order clustered points by (ascending) distance to cluster centre
	for (unsigned int i=0; i < clusters.size(); i++) {
		std::multimap<double,int> clusterDistance;
		for (unsigned int j=0; j < clusters[i].size(); j++) {
			double d2 = (vecRot.row(clusters[i][j]) - clusterCentres.row(i)).squaredNorm();
			clusterDistance.insert(std::make_pair(d2, clusters[i][j]));
		}
		//the map will be sorted based on the key so just loop through it
		//to get set of data indices sorted on the distance to cluster
		clusters[i].clear();
		for (std::multimap<double,int>::iterator it = clusterDistance.begin(); it != clusterDistance.end(); it++) {
			clusters[i].push_back(it->second);
		}
	}

	return clusters;
}

