/*
 * Kmeans.cpp
 *
 *  Created on: 04-Mar-2009
 *      Author: sbutler
 */

#include "Kmeans.h"

#include <random>
#include <iostream>

#include <map>

Kmeans::Kmeans() {
}

Kmeans::~Kmeans() {
}

/**
 * kmeans clustering
 * based on kmeans matlab script by Ian T Nabney (1996-2001)
 *
 * @param data 	affinity matrix
 * @oaram ncentres	number of clusters
 * @return		ordered set of data row indices (the path id) for each cluster
 * 				(order is based on distance to cluster centre)
 */
std::vector<std::vector<int> > Kmeans::cluster(Eigen::MatrixXd& data, int ncentres) {
	int ndims = data.cols();
	int ndata = data.rows();

	//int ncentres = 2;
	Eigen::MatrixXd centres = Eigen::MatrixXd::Zero(ncentres, ndims);
	Eigen::MatrixXd old_centres;

	std::vector<int> rands;
	for (int i=0; i < ncentres; i++) {
		//randomly initialise centers
		bool flag;
		do {
			flag = false;
			int randIndex = rand() % ndata;
			//make sure same row not chosen twice
			for (unsigned int j=0; j < rands.size(); ++j) {
				if (randIndex == rands[j]) {
					flag = true;
					break;
				}
			}
			if (!flag) {
				centres.row(i) = data.row(randIndex);
				rands.push_back(randIndex);
			}
		} while (flag);
	}
	Eigen::MatrixXd id = Eigen::MatrixXd::Identity(ncentres, ncentres);
	//maps vectors to centres.
	Eigen::MatrixXd post(ndata, ncentres);

	Eigen::VectorXd minvals(ndata);

	double old_e = 0;
	int niters = 100;
	for (int n=0; n < niters; n++) {
		//Save old centres to check for termination
		old_centres = centres;

		// Calculate posteriors based on existing centres
		Eigen::MatrixXd d2(ndata, ncentres);
		for (int j = 0; j < ncentres; j++) {
			for (int k=0; k < ndata; k++) {
				d2(k,j) = (data.row(k)-centres.row(j)).squaredNorm();
			}
		}

		int r,c;
		// Assign each point to nearest centre
		for (int k=0; k < ndata; k++) {
			//get centre index (c)
			minvals[k] = d2.row(k).minCoeff(&r, &c);
			//set centre
			post.row(k) = id.row(c);
		}

		Eigen::VectorXd num_points = post.colwise().sum();
		// Adjust the centres based on new posteriors
		for (int j = 0; j < ncentres; j++) {
			if (num_points[j] > 0) {
				Eigen::MatrixXd s = Eigen::MatrixXd::Zero(1,ndims);
				for (int k=0; k<ndata; k++) {
					if (post(k,j) == 1) {
						s += data.row(k);
					}
				}
				centres.row(j) = s / num_points[j];
			}
		}

		// Error value is total squared distance from cluster centres
		double e = minvals.sum();
		double ediff = fabs(old_e - e);
		double cdiff = (centres-old_centres).cwiseAbs().maxCoeff();
		std::cout << "Cycle " << n << " Error " << e << " Movement " << cdiff << ", " << ediff << std::endl;

		if (n > 1) {
			//Test for termination
			if (cdiff < 0.0000000001 && ediff < 0.0000000001) {
				break;
			}
		}
		old_e = e;
	}

	//------- finished kmeans ---------

	//find the item closest to the centre for each cluster
	std::vector<std::vector<int> > clustered_items;
	for (int j=0; j < ncentres; j++) {
		//put data into map (multimap because minvals[k] could be the same for multiple units)
		std::multimap<double,int> cluster;
		for (int k=0; k < ndata; k++) {
			if (post(k,j) == 1) {
				cluster.insert(std::make_pair(minvals[k], k));
			}
		}
		//extract data from map
		std::vector<int> units;
		//the map will be sorted based on the key (the minval) so just loop through it
		//to get set of data indices sorted on the minval
		for (std::multimap<double,int>::iterator it = cluster.begin(); it != cluster.end(); it++) {
			units.push_back(it->second);
		}

		clustered_items.push_back(units);
	}
	//return the ordered set of item indices for each cluster centre
	return clustered_items;
}
