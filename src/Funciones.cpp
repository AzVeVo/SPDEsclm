// [[Rcpp::depends(RcppArmadillo, roptim, RcppProgress)]]

#include <RcppArmadillo.h>
#include <roptim.h>
#include <progress.hpp>
#include <progress_bar.hpp>

using namespace Rcpp;
using namespace arma;
using namespace roptim;


// Optimizer for phi and gamma
//-------------------------------------------------------------------------------------------------------------
class optimQ : public Functor {
private:
  const List mats;    // Matrices for Q
  const arma::mat B;  // EYY matrix
  const arma::vec C;  // EY vector
  const arma::vec D;  // Vector of means
  const double sig;   // Estimate for sigma2
public:
  optimQ(const List a, const arma::mat b, const arma::vec c, const arma::vec d, const double sig)
    : mats(a), B(b), C(c), D(d), sig(sig) {}

  double operator()(const arma::vec& x) override {
    double phie = x(0);
    double gammae = x(1);

    arma::sp_mat D1  = as<arma::sp_mat>(mats["D"]);
    arma::sp_mat G1 = as<arma::sp_mat>(mats["G1"]);
    arma::sp_mat G2 = as<arma::sp_mat>(mats["G2"]);
    arma::sp_mat Aa  = as<arma::sp_mat>(mats["A"]);
    arma::sp_mat Q  = (pow(phie,2.0)/(4.0*M_PI))*(D1/pow(phie,4.0) + 2.0*G1/pow(phie,2.0) + G2);
    Q = 0.5*(Q + Q.t());
    arma::uword p = B.n_rows;
    arma::mat In  = eye(p, p);
    arma::mat K   = arma::mat((1.0-gammae)*Q + gammae*Aa.t()*Aa);
    arma::mat PsiInv1 = (1.0/(1.0-gammae))*In - (gammae/(1.0-gammae))*(arma::mat(Aa))*K.i()*(arma::mat(Aa)).t();
    double ldet, sign;
    log_det(ldet, sign, PsiInv1);
    const double f = 0.5*(-ldet + (trace(B*PsiInv1) - as_scalar(D.t()*PsiInv1*(2.0*C - D)))/sig);

    return f;
  }
};

arma::vec optimlL(arma::vec rhoG, List mats, arma::mat yy1, arma::vec y1, arma::vec media1, double sigma2,
                  arma::vec lower2, arma::vec upper2) {
  arma::vec optR(2, arma::fill::zeros);
  optimQ fun(mats, yy1, y1, media1, sigma2);
  Roptim<optimQ> opt("L-BFGS-B");
  opt.set_lower(lower2);
  opt.set_upper(upper2);
  opt.control.trace = 0;
  arma::vec x = rhoG;
  opt.minimize(fun, x);
  optR = opt.par();
  return optR;
}


// [[Rcpp::export]]
arma::mat crossdist(arma::mat m1){
  arma::uword nrow1 = m1.n_rows;
  arma::mat out(nrow1, nrow1, fill::zeros);

  for (uword r1 = 0; r1<nrow1; r1++) {
    for (uword r2 = (r1+1); r2<nrow1; r2++) {
      out(r1,r2) = sqrt(pow(m1(r1,0)-m1(r2,0),2.0) + pow(m1(r1,1)-m1(r2,1),2.0));
      out(r2,r1) = out(r1,r2);
    }
  }
  return out;
}

// Simulate from a Truncated Normal Distribution
//-------------------------------------------------------------------------------------------------------------
arma::mat rtnormal(int n, arma::vec mu, arma::vec s, arma::mat Rinv, arma::vec a, arma::vec b, int burn, int lag){
  int m = lag*n + burn;
  int p = Rinv.n_cols;
  arma::mat X(n, p, fill::zeros);

  Rcpp::NumericVector l1 = Rcpp::wrap((a - mu)/s);
  Rcpp::NumericVector u1 = Rcpp::wrap((b - mu)/s);
  arma::vec pa = Rcpp::pnorm(l1,0,1,1,0);
  arma::vec pb = Rcpp::pnorm(u1,0,1,1,0);
  arma::vec x0 = randu<arma::vec>(p);
  Rcpp::NumericVector x1 = Rcpp::wrap(pa + (pb - pa)%x0);
  arma::colvec x = Rcpp::qnorm(x1,0,1,1,0);
  arma::vec lower = as<arma::vec>(l1);
  arma::vec upper = as<arma::vec>(u1);

  arma::uvec q1 = find_nonfinite(x);
  x.elem(q1) = lower.elem(q1);
  q1 = find_nonfinite(x);
  x.elem(q1) = upper.elem(q1);

  umat minusj(p-1, p, fill::zeros);
  for(int j=0; j<p; j++){
    int k=0;
    for(int l=0; l<p; l++){
      if(l!=j){
        minusj(k,j) = l;
        k++;
      }
    }
  }
  double delta, kap, mj, tj, lv, rv, xij;
  arma::uvec pj; arma::rowvec a1; arma::vec xj;
  int count = 1;
  for(int i=0; i<m; i++){
    delta = as_scalar(x.t()*Rinv*x);
    kap = -2.0*log(arma::randu<double>()) + delta;
    for(int j=0; j<p; j++){
      pj = minusj.col(j);
      xj = x(pj);
      a1 = xj.t()*Rinv.rows(pj);
      mj = -a1(j)/Rinv(j,j);
      tj = sqrt(mj*mj + (kap-as_scalar(a1.cols(pj)*xj))/Rinv(j,j));
      lv = std::max(lower(j),(mj-tj));
      rv = std::min(upper(j),(mj+tj));
      xij = lv + (rv - lv)*arma::randu<double>();
      x(j) = xij;
    }
    if (i==(burn + count*lag - 1)){
      X.row(count-1) = x.t();
      count++;
    }
  }
  X = X.t();
  X = X.each_col()%s;
  X = (X.each_col() + mu).t();
  X.replace(arma::datum::inf,arma::datum::nan);
  X.replace(-arma::datum::inf,arma::datum::nan);
  return X;
}


// Computes the Correlation matrix approximation by SPDE approach
//-------------------------------------------------------------------------------------------------------------
List corrSpatial(double phi, List mats) {
  arma::sp_mat D  = as<arma::sp_mat>(mats["D"]);
  arma::sp_mat G1 = as<arma::sp_mat>(mats["G1"]);
  arma::sp_mat G2 = as<arma::sp_mat>(mats["G2"]);
  arma::sp_mat A  = as<arma::sp_mat>(mats["A"]);

  // Precision Matrix
  arma::sp_mat Q = (pow(phi,2.0)/(4.0*M_PI))*(D/pow(phi,4.0) + 2.0*G1/pow(phi,2.0) + G2);
  Q = 0.5*(Q + Q.t());
  arma::mat tempQ = arma::mat(Q);

  //Correlation Matrix
  arma::mat QInv = tempQ.i();
  arma::mat R    = A*QInv*A.t();
  R = 0.5*(R + R.t());

  List output;
  output["R"] = R;
  output["Q"] = Q;
  output["Qinv"] = QInv;
  return output;
}


// Compute the inverse of Psi
//-------------------------------------------------------------------------------------------------------------
List PsiInverse(double gamma, List corr, arma::sp_mat A, arma::uvec ind0) {
  arma::mat R    = as<arma::mat>(corr["R"]);
  arma::sp_mat Q = as<arma::sp_mat>(corr["Q"]);
  arma::uword p  = R.n_rows;
  arma::mat In   = eye(p, p);

  // Compute Psi
  arma::mat Psi = gamma*R + (1.0 - gamma)*In;
  Psi = 0.5*(Psi + Psi.t());

  // Compute PsiInv
  arma::mat K   = arma::mat((1.0-gamma)*Q + gamma*(A.t()*A));
  arma::mat PsiInv = (1.0/(1.0 - gamma))*In - (gamma/(1.0-gamma))*A*K.i()*A.t();
  PsiInv = 0.5*(PsiInv + PsiInv.t());

  // Compute PsiInv_00
  arma::mat Koo    = (1.0-gamma)*arma::mat(Q) + gamma*(arma::mat(A).rows(ind0)).t()*(arma::mat(A).rows(ind0));
  arma::mat PsiInvoo = (1.0/(1.0 - gamma))*In(ind0,ind0) - (gamma/(1.0-gamma))*(arma::mat(A).rows(ind0))*Koo.i()*(arma::mat(A).rows(ind0)).t();
  PsiInvoo = 0.5*(PsiInvoo + PsiInvoo.t());

  List output;
  output["Psi"] = Psi;
  output["PsiInv"] = PsiInv;
  output["PsiInvoo"] = PsiInvoo;
  return output;
}


// SAEM - Estimate parameters in Gaussian spatial model using SPDE approach
// ------------------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
List saemEsp(arma::vec y, arma::mat X, arma::vec cc, arma::mat coords, double init_phi, double init_gamma,
             arma::vec lower, arma::vec upper,  arma::vec lowerp, arma::vec upperp,
             arma::uword M, arma::uword Maxiter, double pc, double tol, List mats, bool infMat){

  // Auxiliary variables
  Progress time(Maxiter, true);
  arma::uword p = y.size();
  arma::uword q = X.n_cols;
  arma::uvec indexs = arma::regspace<arma::uvec>(0,1,q-1);
  arma::uvec ind0   = find(cc==0);  // Index of uncensored observations
  arma::uvec ind1   = find(cc==1);  // Index of censored observations
  arma::uword p1    = ind1.size();  // Number of censored observations
  arma::vec lower1  = lower(ind1);  // Lower bound of censored observations
  arma::vec upper1  = upper(ind1);  // Upper bound of censored observations
  arma::vec mu21(p1, fill::zeros);  // Conditional mean
  arma::mat Sc(p1, p1, fill::zeros);// Conditional variance matrix
  arma::mat derv(p, p, fill::zeros);
  arma::mat aux1(p, p, fill::zeros);
  arma::mat aux2(p, p, fill::zeros);
  arma::vec aux3(p, fill::zeros);

  // Initial values
  arma::vec beta(q);  arma::vec media(p);  double sigma2 = 0.0;
  if (y.has_nan()) {
    arma::uvec indFin = find_finite(y);
    beta   = (((X.rows(indFin)).t()*X.rows(indFin)).i())*(X.rows(indFin)).t()*y(indFin);
    media  = X*beta;
    sigma2 = as_scalar(sum(pow(y(indFin)-media(indFin),2.0)))/(indFin.size());
  } else {
    beta   = ((X.t()*X).i())*X.t()*y;
    media  = X*beta;
    sigma2 = as_scalar(sum(pow(y-media,2.0)))/p;
  }
  double phi   = init_phi;
  double gamma = init_gamma;

  // Estimates
  arma::vec theta(q+3, fill::zeros);
  arma::vec theta1(q+3, fill::zeros);
  theta(indexs)=beta; theta(q)=sigma2; theta(q+1)=phi; theta(q+2)=gamma;
  arma::mat Theta = theta.t();

  // Variance-covariance matrix
  arma::mat In = eye(p, p);
  arma::sp_mat A  = as<arma::sp_mat>(mats["A"]);
  arma::sp_mat D  = as<arma::sp_mat>(mats["D"]);
  arma::sp_mat G2 = as<arma::sp_mat>(mats["G2"]);
  List corr       = corrSpatial(phi, mats);
  arma::mat R     = as<arma::mat>(corr["R"]);
  arma::mat Qinv  = as<arma::mat>(corr["Qinv"]);
  arma::sp_mat Q  = as<arma::sp_mat>(corr["Q"]);
  List covInver   = PsiInverse(gamma, corr, A, ind0);
  arma::mat Psi   = as<arma::mat>(covInver["Psi"]);
  arma::mat PsiInv   = as<arma::mat>(covInver["PsiInv"]);
  arma::mat PsiInvoo = as<arma::mat>(covInver["PsiInvoo"]);

  // Initial values for the optimization function
  arma::vec optP(2, fill::zeros);
  optP(0) = phi;
  optP(1) = gamma;

  // Stopping criteria
  double criterio   = 10.0;
  arma::uword count = 0;

  // SAEM ALGORITHM ------------------------------------------------------------
  arma::mat gibbs(1, p1, fill::zeros);
  arma::vec SAEMY(p, fill::zeros);      // Estimate of first conditional moment
  arma::mat SAEMYY(p, p, fill::zeros);  // Estimate of second conditional moment
  arma::vec EY(p, fill::zeros);
  arma::mat EYY(p, p, fill::zeros);
  arma::vec auxY = y;

  arma::mat ESS(q+3, q+3, fill::zeros);
  arma::mat SAEM_SS(q+3, q+3, fill::zeros);
  arma::vec score(q+3, fill::zeros);
  arma::vec auxMean(p, fill::zeros);
  arma::mat IF(q+3, q+3, fill::zeros);
  arma::vec s(p1, fill::zeros);
  arma::mat Rc(p1, p1, fill::zeros);
  arma::mat Rinv(p1, p1, fill::zeros);

  // lambda sequence
  arma::vec lambda(Maxiter, fill::ones);
  if (pc < 1) {
    arma::vec nonM  = 1.0/arma::regspace<arma::vec>(1, 1, Maxiter*(1-pc));
    arma::uvec nonI = arma::regspace<arma::uvec>(Maxiter-nonM.n_elem, 1, Maxiter-1);
    lambda(nonI)    = nonM;
  }

  while (criterio>tol) {
    time.increment();
    count += 1;

    // Update auxiliary parameters
    mu21 = media(ind1) + Psi(ind1,ind0)*PsiInvoo*(y(ind0) - media(ind0));
    Sc   = sigma2*(Psi(ind1,ind1) - Psi(ind1,ind0)*PsiInvoo*Psi(ind0,ind1));
    Sc   = 0.50*(Sc + Sc.t());
    if (infMat) {
      derv = gamma*arma::mat(A*Qinv*(D/(2.0*M_PI*pow(phi, 3.0)) - phi/(2.0*M_PI)*G2)*Qinv*A.t());
      aux1 = PsiInv*derv;
      aux2 = PsiInv*(R - In);
    }

    // Simulation E-step
    s  = sqrt(Sc.diag());
    Rc = Sc%(1.0/(s * s.t()));
    Rinv = Rc.i();
    for (uword i=0; i<M; i++){
      //int MO = round(M*0.2);
      gibbs  = rtnormal(1, mu21, s, Rinv, lower1, upper1, 3, 1); //MO, 5);
      auxY(ind1) = (gibbs).t();
      EY  = EY + auxY;
      EYY = EYY + auxY*auxY.t();

      if (infMat) {
        auxMean = auxY - media;
        aux3    = PsiInv*auxMean;
        score(indexs) = (X.t()*aux3)/sigma2;
        score(q)   = (-0.5*p)/sigma2 + (0.5/pow(sigma2, 2.0))*as_scalar(auxMean.t()*aux3);
        score(q+1) = (-0.5)*trace(aux1) + (0.5/sigma2)*as_scalar(auxMean.t()*aux1*aux3);
        score(q+2) = (-0.5)*trace(aux2) + (0.5/sigma2)*as_scalar(auxMean.t()*aux2*aux3);
        ESS = ESS + score*score.t();
      }
    }

    // Stochastic Approximation E-step
    SAEMY  = SAEMY  + lambda(count-1)*(EY/M - SAEMY);
    SAEMYY = SAEMYY + lambda(count-1)*(EYY/M - SAEMYY);
    SAEMYY = 0.50*(SAEMYY + SAEMYY.t());
    if (infMat) {
      SAEM_SS = SAEM_SS + lambda(count-1)*(ESS/M - SAEM_SS);
    }

    // Maximization M-step
    beta   = (X.t()*PsiInv*X).i()*X.t()*PsiInv*SAEMY;
    media  = X*beta;
    sigma2 = (trace(SAEMYY*PsiInv) - as_scalar(media.t()*PsiInv*(2.0*SAEMY - media)))/p;
    optP   = optimlL(optP, mats, SAEMYY, SAEMY, media, sigma2, lowerp, upperp);
    phi    = optP(0);
    gamma  = optP(1);

    // Update variance-covariance matrices
    corr = corrSpatial(phi, mats);
    R    = as<arma::mat>(corr["R"]);
    Qinv = as<arma::mat>(corr["Qinv"]);
    Q    = as<arma::sp_mat>(corr["Q"]);
    covInver = PsiInverse(gamma, corr, A, ind0);
    Psi      = as<arma::mat>(covInver["Psi"]);
    PsiInv   = as<arma::mat>(covInver["PsiInv"]);
    PsiInvoo = as<arma::mat>(covInver["PsiInvoo"]);

    // Stopping criteria
    theta1(indexs) = beta;  theta1(q) = sigma2;  theta1(q + 1) = phi;  theta1(q + 2) = gamma;
    criterio = as_scalar(sqrt(sum((theta1/theta - 1)%(theta1/theta - 1))));
    if (count == Maxiter) criterio = 1e-12;
    theta = theta1;
    Theta = join_vert(Theta, theta.t());

    // Reset auxiliary variables
    EY.zeros(p);
    EYY.zeros(p,p);
    if (infMat) { ESS.zeros(q+3, q+3); }
  }

  // Results
  List output;
  output["Theta"] = Theta;   output["theta"] = theta;   output["beta"]  = beta;
  output["sigma2"]= sigma2;  output["phi"]   = phi;     output["gamma"] = gamma;
  output["EY"]    = SAEMY;   output["EYY"]   = SAEMYY;  output["Psi"] = Psi;

  // Information Matrix
  if (infMat) {
    auxMean = SAEMY - media;
    arma::vec diff = 2.0*SAEMY - media;
    arma::mat d1   = Qinv*(D/(2.0*M_PI*pow(phi, 3.0)) - phi*G2/(2.0*M_PI))*Qinv;
    derv = gamma*A*d1*A.t();
    aux1 = PsiInv*derv;
    aux2 = PsiInv*(R - In);
    arma::mat aux4 = aux1*PsiInv;
    arma::mat aux5 = aux2*PsiInv;
    arma::mat d2   = d1*(D/(2.0*M_PI*pow(phi, 3.0)) - phi*G2/(2.0*M_PI))*Qinv - Qinv*(1.5*D/(M_PI*pow(phi, 4.0)) + 0.5*G2/M_PI)*Qinv + Qinv*(D/(2.0*M_PI*pow(phi, 3.0)) - phi*G2/(2.0*M_PI))*d1;
    arma::mat aux6 = gamma*PsiInv*A*(d2)*A.t();

    // Q'(theta)
    score(indexs) = X.t()*PsiInv*auxMean/sigma2;
    score(q)   = (-0.5*p)/sigma2 + (0.5/pow(sigma2, 2.0))*as_scalar(trace(SAEMYY*PsiInv) - media.t()*PsiInv*diff);
    score(q+1) = (-0.5)*trace(aux1) + (0.5/sigma2)*as_scalar(trace(SAEMYY*aux4) - media.t()*aux4*diff);
    score(q+2) = (-0.5)*trace(aux2) + (0.5/sigma2)*as_scalar(trace(SAEMYY*aux5) - media.t()*aux5*diff);

    // Q''(theta)
    arma::mat Hessian(q+3, q+3, fill::zeros);
    Hessian(indexs, indexs)   = (-X.t()*PsiInv*X)/sigma2;
    Hessian.submat(0,q,q-1,q) = (-X.t()*PsiInv*auxMean)/pow(sigma2, 2.0);
    Hessian.submat(q,0,q,q-1) = (Hessian.submat(0,q,q-1,q)).t();
    Hessian.submat(0,q+1,q-1,q+1) = (-X.t()*aux4*auxMean)/sigma2;
    Hessian.submat(q+1,0,q+1,q-1) = (Hessian.submat(0,q+1,q-1,q+1)).t();
    Hessian.submat(0,q+2,q-1,q+2) = (-X.t()*aux5*auxMean)/sigma2;
    Hessian.submat(q+2,0,q+2,q-1) = (Hessian.submat(0,q+2,q-1,q+2)).t();

    Hessian(q,q)   = (0.5*p/pow(sigma2, 2.0)) - as_scalar(trace(SAEMYY*PsiInv) - media.t()*PsiInv*diff)/pow(sigma2, 3.0);
    Hessian(q,q+1) = Hessian(q+1,q) = (-0.5/pow(sigma2, 2.0))*as_scalar(trace(SAEMYY*aux4) - media.t()*aux4*diff);
    Hessian(q,q+2) = Hessian(q+2,q) = (-0.5/pow(sigma2, 2.0))*as_scalar(trace(SAEMYY*aux5) - media.t()*aux5*diff);

    Hessian(q+1,q+1) = 0.5*trace(aux1*aux1 - aux6) - (0.5/sigma2)*as_scalar(trace(SAEMYY*(2.0*aux1*aux4 - aux6*PsiInv)) - media.t()*(2.0*aux1*aux4 - aux6*PsiInv)*diff);
    Hessian(q+1,q+2) = Hessian(q+2,q+1) = 0.5*trace(aux2*aux1 - aux1/gamma) - (0.5/sigma2)*as_scalar(trace(SAEMYY*(aux2*aux4 - aux4/gamma + aux1*aux5)) - media.t()*(aux2*aux4 - aux4/gamma + aux1*aux5)*diff);
    Hessian(q+2,q+2) = 0.5*trace(aux2*aux2) - as_scalar(trace(SAEMYY*aux2*aux5) - media.t()*(aux2*aux5)*diff)/sigma2;

    IF = score*score.t() - Hessian - SAEM_SS;
    arma::mat IFinv = IF.i();
    output["SE"] = sqrt(IFinv.diag());
    output["InfMat"] = 0.5*(IF + IF.t());
  }
  output["Iterations"] = count;
  return output;
}


List corrSpatialPred(double phi, List mats, arma::sp_mat Apred) {
  arma::sp_mat D  = as<arma::sp_mat>(mats["D"]);
  arma::sp_mat G1 = as<arma::sp_mat>(mats["G1"]);
  arma::sp_mat G2 = as<arma::sp_mat>(mats["G2"]);
  arma::sp_mat A  = as<arma::sp_mat>(mats["A"]);
  A = join_cols(A, Apred);

  // Precision Matrix
  arma::sp_mat Q = (pow(phi,2.0)/(4.0*M_PI))*(D/pow(phi,4.0) + 2.0*G1/pow(phi,2.0) + G2);
  Q = 0.5*(Q + Q.t());
  arma::mat tempQ = arma::mat(Q);

  //Correlation Matrix
  arma::mat QInv = inv_sympd(tempQ);
  arma::mat R    = A*QInv*A.t();
  R = 0.5*(R + R.t());

  List output;
  output["R"] = R;
  output["Q"] = Q;
  output["Qinv"] = QInv;
  output["A"] = A;
  return output;
}

// Prediction
// ------------------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
List prediccion (List model, List mats, arma::sp_mat Apred, arma::mat Xnew, arma::mat coordnew, arma::vec cpred) {

  arma::vec beta    = as<arma::mat>(model["beta"]);
  arma::mat Xo      = as<arma::mat>(model["X"]);
  arma::vec EY      = as<arma::mat>(model["EY"]);
  arma::vec meano   = Xo*beta;
  arma::vec meannew = Xnew*beta;
  arma::mat coord   = as<arma::mat>(model["coord"]);
  arma::uvec ind0   = find(cpred==0);
  arma::uvec ind1   = find(cpred==1);

  double phi    = as<double>(model["phi"]);
  double sigma2 = as<double>(model["sigma2"]);
  double gamma  = as<double>(model["gamma"]);

  List corr      = corrSpatialPred(phi, mats, Apred);
  arma::sp_mat A = as<arma::sp_mat>(corr["A"]);
  List covInver  = PsiInverse(gamma, corr, A, ind0);

  arma::mat Psi      = as<arma::mat>(covInver["Psi"]);
  arma::mat PsiInv   = as<arma::mat>(covInver["PsiInv"]);
  arma::mat PsiInvoo = as<arma::mat>(covInver["PsiInvoo"]);

  arma::vec yPred   = meannew + Psi(ind1,ind0)*PsiInvoo*(EY(ind0) - meano);
  arma::mat varPred = sigma2*(Psi(ind1,ind1) - Psi(ind1,ind0)*PsiInvoo*Psi(ind0,ind1));
  arma::vec sdyPred = sqrt(varPred.diag());

  List output;
  output["coord"]      = coordnew;
  output["predValues"] = yPred;
  output["sdPred"]     = sdyPred;
  return output;
}


