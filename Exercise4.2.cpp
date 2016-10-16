// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <math.h>

#define lamda 0.000001
#define step 0.001
#define max_iteration 1000000

using namespace Rcpp;
typedef Eigen::SparseVector<double> Svd;
typedef Svd::InnerIterator InIterVec;
typedef Eigen::VectorXd Vd;
typedef Eigen::VectorXi Vi;
typedef Eigen::MappedSparseMatrix<double> MSpMat;
typedef MSpMat::InnerIterator InIterMat;


void SGD(const MSpMat& X, const Vd& Y, Vd& beta, int index, Vi& lastUpdate, int iteration){
  double y = Y(index);

  double xb = 0.0;
  double temp = 0.0;

  //compute x*b
  for (InIterMat i_(X, index); i_; ++i_){
    temp = beta(i_.index()) * pow(1 - lamda, iteration - lastUpdate(i_.index()));//optimize by precompute power of 1 - lamda
 
    beta(i_.index()) = temp;
    
    lastUpdate(i_.index()) = iteration;
    xb += i_.value() * temp;
  }
  xb += beta(X.rows());
  
  temp = 1/(1 + exp(-xb)) - y;

  //update beta
  for (InIterMat i_(X, index); i_; ++i_){
    beta(i_.index()) = beta(i_.index()) - step * i_.value() * temp;
  }
  beta(X.rows()) =  beta(X.rows()) - step * temp;

  if(iteration == max_iteration - 1){
    for(int i = 0; i< X.rows() ;i++){
      beta(i) = beta(i) * pow(1 - lamda, iteration - lastUpdate(i));//optimize by precompute power of 1 - lamda
    }
  }
}

double predict(const MSpMat& X, const Vd& Y, const Vd& beta){
  double temp = 0.0;
  int correct = 0;
  double truth = 0.0;
  for(int column = 0; column < X.cols(); column++){
    double xb = 0.0;
    for (InIterMat i_(X, column); i_; ++i_){
       xb += i_.value() * beta(i_.index());
    }
    xb += beta(X.rows());
    temp = 1/(1 + exp(-xb));
    truth = Y(column);
    if(temp < 0.5&& Y(column) < 0.5 || temp>0.5 && Y(column) > 0.5){
      correct++;
    }
  }
  return correct * 1.0/X.cols();
}

// [[Rcpp::export]]
void mymain(const MSpMat& X, const Vd& Y, const MSpMat& testX, const Vd& testY){
  Vd beta(X.rows()+1);
  Vi lastUpdate(X.rows()+1);
  
  for(int i=0; i<X.rows()+1;i++){
    beta(i) = 0.0;
    lastUpdate(i) = 0;
  }
  
  for(int iteration = 0; iteration < max_iteration; iteration++){
    int index = iteration%(X.cols());
    SGD(X, Y, beta, index, lastUpdate, iteration);
  }

  double trainAcc = predict(X,Y,beta);
  double testAcc = predict(testX,testY,beta);
  Rcout<< trainAcc << " " << testAcc<< std::endl;
  //return beta;
}


/*** R
library(Matrix)

#sourceCpp("/Users/dgy/Desktop/SDS385/R/Exercise4.2.cpp")
require(stream)
X = readRDS("/Users/dgy/Desktop/SDS385/binaryURL/url_X.rds")
Y = readRDS("~/Desktop/SDS385/binaryURL/url_y.rds")
X = t(X)

order <- sample.int(dim(X)[2])
test <- order[1:2000]
train <- order[2001:length(order)]

beta = mymain(X[,train],Y[train], X[,test],Y[test])


#Xstream <- DSD_Memory(X, loop= TRUE, class = Y)

#implement AdaGrad
*/
