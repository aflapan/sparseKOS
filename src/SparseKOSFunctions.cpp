#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <cmath>
using namespace Rcpp;
using namespace std;
using namespace arma;

// [[Rcpp::export]]
int signCPP(double x){
  if(x>0){
    return(1);
  }
  else if(x<0){
    return(-1);
  }
  else return(0);
}

// [[Rcpp::export]]
double SoftThreshCPP(double x,double lambda){
  double t=0;
  t=signCPP(x)*fmax(abs(x)-lambda,0);
  return(t);
}


// [[Rcpp::export]]
arma::Col<double> CoordDesCPP(arma::Col<double> w0, arma::Mat<double> Q, arma::Col<double> beta,double Lambda,double Epsilon, int Maxniter=1.0e7){
  int p=w0.size();
  int niter=0; // Declare iteration number.
  double error=1; // Declare error term. This will measure convergence of weigh
  arma::Col<double> w=w0; //Declare updated weight and initialize it to old weights.
  double normDiff=0;
  while((error>Epsilon)&&(niter<Maxniter)){
    normDiff=0;
    for(int i=0;i<p;i++){
      if (w0 [i] != 0){
        if(beta[i] == 0){
          w[i]=0;
        }
        else if(Q(i,i)==0){
          w[i]=signCPP(beta[i]);
         // std::cout<<i<<"-th Q Diagonal Term is zero without previous weight being zero." <<"\n";
        }
        else{
          double a=beta(i)-dot(Q.row(i),w)+Q(i,i)*w(i);
          double b=Lambda/2;
          w[i]=SoftThreshCPP(a,b);
          w[i]=w[i]/Q(i,i);
          if(w[i]>1){
            w[i]=1;
          }
          else if(-w[i]>1){
            w[i]=-1;
          }
        }

      }
      normDiff=normDiff+pow(w[i]-w0[i],2);
    }
    w0=w;
    niter++;
    error=sqrt(normDiff);
  }
          return(w);
}

// [[Rcpp::export]]
arma::Col<double> SolveKOSCPP(arma::Mat<double> YTheta,arma::Mat<double> K,double Gamma,double Epsilon=1e-5){
   int n=K.n_rows;
   arma::Mat<double> C=arma::eye(n,n)-(1/n)*arma::ones(n,n); //Create column-centering matrix
   arma::Mat<double> M=(C*K)*C; // This doubly-centered kernel matrix appears often. Compute it once to save time
   arma::Mat<double> Mat=(M*M+ n*Gamma*(M+Epsilon*arma::diagmat(arma::ones(n)))); //Add the epsilon*I term for stability
   arma::Mat<double> A=M*YTheta;
   arma::Col<double> Dvec=arma::solve(Mat,A);
   return(Dvec);
 }

// [[Rcpp::export]]
arma::Mat<double> DerivCPP(arma::Row<double> x, arma::Mat<double> Data, arma::Col<double> w0, double sigmaD){
  arma::Mat<double> Diff=Data.each_row()-x;
  Diff=square(Diff);

  arma::Mat<double> DiffPart = (-2/(sigmaD))*(Diff*diagmat(w0));
  int k=Data.n_cols;
  arma::Col<double> RowSums=Diff*arma::ones(k);
  arma::Col<double> KernPart=exp(-RowSums*(1/(sigmaD)));
  arma::Mat<double> result=arma::diagmat(KernPart)*DiffPart;
  return(result);
}

// [[Rcpp::export]]
arma::Mat<double> TMatCPP(arma::Mat<double> Data, arma::Col<double> A, arma::Col<double> w0, double sigmaTm){
  int p=Data.n_cols;
  int n=A.n_rows;
  arma::Mat<double> C=arma::eye(n,n)-(1/n)*arma::ones(n,n);
  A=C*A;
  arma::Mat<double> T=arma::zeros(n,p);
  for(int j=0; j<n; j++){
    T.row(j)=A.t()*DerivCPP(Data.row(j),Data,w0,sigmaTm);
  }
  T=C*T;
  return(T);
}

// [[Rcpp::export]]
double ObjectiveFuncCPP(arma::Col<double> w, arma::Mat<double> KwOF, arma::Mat<double> Data, arma::Col<double> DVectors, arma::Mat<double> YTheta, double LambdaOF, double GammaOF, double EpsilonOF=1e-5){
  int n = KwOF.n_rows;
  arma::Mat<double> C=arma::eye(n,n)-(1/n)*arma::ones(n,n);
  arma::Mat<double> M=C*(KwOF*C); //Form the doubly-centered kernel matrix.
  double Regression = (1/n)*accu(square(YTheta-(M*DVectors))); //Least Squares component
  double LASSO = LambdaOF*sum(abs(w)); //LASSO Component
  arma::Mat<double> RidgeMat=DVectors.t()*(M*DVectors+EpsilonOF*DVectors);
  double Ridge = GammaOF*accu(arma::diagmat(RidgeMat)); //Ridge component
  return(Regression+LASSO+Ridge); //The sum of all three gives the objective function
}
