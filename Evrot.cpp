/*
 * evrot.cpp
 *
 *  Created on: 04-Mar-2009
 *      Author: sbutler
 */

#include "Evrot.h"

#include <map>

Evrot::Evrot(Eigen::MatrixXd& X, int method):
	mMethod(method),
	mNumDims(X.cols()),
	mNumData(X.rows()),
	mNumAngles((int)(mNumDims*(mNumDims-1)/2)), // get the number of angles
	ik(Eigen::VectorXi(mNumAngles)),
	jk(Eigen::VectorXi(mNumAngles)),
	mX(X),
	mClusters(std::vector<std::vector<int> >(mNumDims)) //allocate clusters vector
{
	// build index mapping (to index upper triangle)
	int k = 0;
	for( int i=0; i<mNumDims-1; i++ ){
		for( int j=i+1; j<=mNumDims-1; j++ ){
			ik[k] = i;
			jk[k] = j;
			k++;
		}
	}

	evrot();
}

Evrot::~Evrot() {

}

void Evrot::evrot() {

	// definitions
	int max_iter = 200;
	double dQ,Q,Q_new,Q_old1,Q_old2,Q_up,Q_down;
	double alpha;
	int iter,d;

	Eigen::VectorXd theta = Eigen::VectorXd::Zero(mNumAngles);
	Eigen::VectorXd theta_new = Eigen::VectorXd::Zero(mNumAngles);

	Q = evqual(mX); // initial quality
	if( DEBUG )
		std::cout << "Q = " << Q << std::endl;
	Q_old1 = Q;
	Q_old2 = Q;
	iter = 0;
	while( iter < max_iter ){ // iterate to refine quality
		iter++;
		for( d = 0; d < mNumAngles; d++ ){
			if( mMethod == 2 ){ // descend through numerical drivative
				alpha = 0.1;
				{
					// move up
					theta_new[d] = theta[d] + alpha;
					Eigen::MatrixXd& Xrot = rotate_givens(theta_new);
					Q_up = evqual(Xrot);
					delete &Xrot;
				}
				{
					// move down
					theta_new[d] = theta[d] - alpha;
					Eigen::MatrixXd& Xrot = rotate_givens(theta_new);
					Q_down = evqual(Xrot);
					delete &Xrot;
				}

				// update only if at least one of them is better
				if( Q_up > Q || Q_down > Q){
					if( Q_up > Q_down ){
						theta[d] = theta[d] + alpha;
						theta_new[d] = theta[d];
						Q = Q_up;
					} else {
						theta[d] = theta[d] - alpha;
						theta_new[d] = theta[d];
						Q = Q_down;
					}
				}
			} else { // descend through true derivative
				alpha = 1.0;
				dQ = evqualitygrad(theta, d);
				theta_new[d] = theta[d] - alpha * dQ;
				Eigen::MatrixXd& Xrot = rotate_givens(theta_new);
				Q_new = evqual(Xrot);
				delete &Xrot;

				if( Q_new > Q){
					theta[d] = theta_new[d];
					Q = Q_new;
				}
				else{
					theta_new[d] = theta[d];
				}
			}
		}
		// stopping criteria
		if( iter > 2 ){
			if( Q - Q_old2 < 1e-3 ){
				break;
			}
		}
		Q_old2 = Q_old1;
		Q_old1 = Q;
	}

	if (DEBUG)
		std::cout << "Done after " << iter << " iterations, Quality is " << Q << std::endl;

	mXrot = rotate_givens(theta_new);
	cluster_assign();

	//output
	mQuality = Q;
}

void Evrot::cluster_assign() {
	// find max of each row
	Eigen::VectorXi max_index_col(mNumData);
	int i,j;
	for (int i=0; i<mNumData; i++ ) {
		int row, col;
		mXrot.row(i).cwise().abs().maxCoeff(&row, &col);
		max_index_col[i] = col;
	}

	// prepare cluster assignments
	for( j=0; j<mNumDims; j++ ){  // loop over all columns
		for( i=0; i<mNumData; i++ ){ // loop over all rows
			if( max_index_col[i] == j ){
				mClusters[j].push_back(i);
			}
		}
	}
}

double Evrot::evqual(Eigen::MatrixXd& X) {
	// take the square of all entries and find max of each row
	Eigen::MatrixXd X2 = X.cwise().pow(2);
	Eigen::VectorXd max_values = X2.rowwise().maxCoeff();

	// compute cost
	for (int i=0; i<mNumData; i++ ) {
		X2.row(i) = X2.row(i) / max_values[i];
	}
	double J = 1.0 - (X2.sum()/mNumData -1.0)/mNumDims;
	if( DEBUG )
		std::cout << "Computed quality = "<< J << std::endl;

	return J;
}

double Evrot::evqualitygrad(Eigen::VectorXd& theta, int angle_index) {
	// build V,U,A
	Eigen::MatrixXd& V = gradU(theta, angle_index);

	Eigen::MatrixXd& U1 = build_Uab(theta, 0,angle_index-1);
	Eigen::MatrixXd& U2 = build_Uab(theta, angle_index+1,mNumAngles-1);

	Eigen::MatrixXd A = mX*U1*V*U2;

	delete &V;
	delete &U1;
	delete &U2;

	// rotate vecs according to current angles
	Eigen::MatrixXd& Y = rotate_givens(theta);

	// find max of each row
	Eigen::VectorXd max_values(mNumData);
	Eigen::VectorXi max_index_col(mNumData);
	for (int i=0; i<mNumData; i++ ) {
		int row, col;
		Y.row(i).cwise().abs().maxCoeff(&row, &col);
		max_values[i] = Y(i,col);
		max_index_col[i] = col;
	}

	// compute gradient
	double dJ=0, tmp1, tmp2;
	for( int j=0; j<mNumDims; j++ ){  // loop over all columns
		for( int i=0; i<mNumData; i++ ){ // loop over all rows
			tmp1 = A(i,j) * Y(i,j) / (max_values[i]*max_values[i]);
			tmp2 = A(i,max_index_col[i]) * (Y(i,j)*Y(i,j)) / (max_values[i]*max_values[i]*max_values[i]);
			dJ += tmp1-tmp2;
		}
	}
	dJ = 2*dJ/mNumData/mNumDims;
	if( DEBUG )
		std::cout << "Computed gradient = " << dJ << std::endl;

	delete &Y;

	return dJ;
}

Eigen::MatrixXd& Evrot::rotate_givens(Eigen::VectorXd& theta) {
	Eigen::MatrixXd& G = build_Uab(theta, 0, mNumAngles-1);
	Eigen::MatrixXd& Y = *new Eigen::MatrixXd(mX.rows(),mX.cols());
	Y = mX*G;
	delete &G;
	return Y;
}

Eigen::MatrixXd& Evrot::build_Uab(Eigen::VectorXd& theta, int a, int b) {
	int k,i;
	//set Uab to be an identity matrix
	Eigen::MatrixXd& Uab = *new Eigen::MatrixXd(mNumDims,mNumDims);
	Uab.setZero();
	Uab.setIdentity();

	if( b < a ) {
		return Uab;
	}

	double tt,u_ik;
	for( k=a; k<=b; k++ ){
		tt = theta[k];
		for( i=0; i<mNumDims; i++ ) {
			u_ik = 			Uab(i,ik[k]) * cos(tt) - Uab(i,jk[k]) * sin(tt);
			Uab(i,jk[k]) = 	Uab(i,ik[k]) * sin(tt) + Uab(i,jk[k]) * cos(tt);
			Uab(i,ik[k]) = u_ik;
		}
	}
	return Uab;
}

Eigen::MatrixXd& Evrot::gradU(Eigen::VectorXd& theta, int k) {
	Eigen::MatrixXd& V = *new Eigen::MatrixXd(mNumDims,mNumDims);
	V.setZero();

	V(ik[k],ik[k]) = -sin(theta[k]);
	V(ik[k],jk[k]) = cos(theta[k]);
	V(jk[k],ik[k]) = -cos(theta[k]);
	V(jk[k],jk[k]) = -sin(theta[k]);

	return V;
}
