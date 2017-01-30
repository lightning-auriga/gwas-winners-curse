/*
  Copyright 2016 Cameron Palmer

  This file is part of gwas-winners-curse.
  
  gwas-winners-curse is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  gwas-winners-curse is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with gwas-winners-curse.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "gwas_winners_curse/zero_finding.h"


double gwas_winners_curse::ci_log_likelihood(double x, void *params) {
  struct gwas_winners_curse::ci_log_likelihood_params *p = reinterpret_cast<ci_log_likelihood_params *>(params);

  double beta = p->beta_debiased;
  double stderr = p->stderr_debiased;
  double p_cut = p->p_threshold;
  double cut = qnorm(1.0 - p_cut / 2.0);
  double xlog = 0.0;//log(dnorm((x - beta) / stderr) / (stderr * (pnorm(beta / stderr - cut) + pnorm(-beta / stderr - cut))));
  double betalog = 0.0;//log(dnorm(0) / (stderr * (pnorm(beta / stderr - cut) + pnorm(-beta / stderr - cut))));
  if (parameters::get_flag("invalid-approximation")) {
    xlog = log(dnorm((x - beta) / stderr) / (stderr * (pnorm(fabs(beta) / stderr - cut))));
    betalog = log(dnorm(0) / (stderr * (pnorm(fabs(beta) / stderr - cut))));
  } else {
    xlog = log(dnorm((x - beta) / stderr) / (stderr * (pnorm(beta / stderr - cut) + pnorm(-beta / stderr - cut))));
    betalog = log(dnorm(0) / (stderr * (pnorm(beta / stderr - cut) + pnorm(-beta / stderr - cut))));
  }
  return xlog - betalog + qchisq(0.95)/2.0;
}

std::pair<double, double> gwas_winners_curse::zero_finding::calculate_ci(const double &beta_debiased,
									    const double &stderr_debiased,
									    const double &p_threshold) {
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T = 0;
  gsl_root_fsolver *s = 0;
  double factor = 10.0;
  double x_lo = beta_debiased - stderr_debiased * factor, x_hi = beta_debiased;
  double r = 0.0;
  gsl_function F;

  std::pair<double, double> res;
  
  struct gwas_winners_curse::ci_log_likelihood_params params = {beta_debiased, stderr_debiased, p_threshold};
    
  F.function = &gwas_winners_curse::ci_log_likelihood;
  F.params = &params;
  
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);
  try {
    gsl_root_fsolver_set(s, &F, x_lo, x_hi);
    
    do {
      iter++;
      status = gsl_root_fsolver_iterate(s);
      r = gsl_root_fsolver_root(s);
      x_lo = gsl_root_fsolver_x_lower(s);
      x_hi = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(x_lo, x_hi,
				      0, 0.001);
    } while (status == GSL_CONTINUE && iter < max_iter);
    
    gsl_root_fsolver_free(s);
  } catch (...) {
    gsl_root_fsolver_free(s);
    throw;
  }
  res.first = r;

  r = 0.0;
  x_lo = beta_debiased;
  x_hi = beta_debiased + stderr_debiased * factor;
  s = gsl_root_fsolver_alloc(T);
  try {
    gsl_root_fsolver_set(s, &F, x_lo, x_hi);
    
    do {
      iter++;
      status = gsl_root_fsolver_iterate(s);
      r = gsl_root_fsolver_root(s);
      x_lo = gsl_root_fsolver_x_lower(s);
      x_hi = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(x_lo, x_hi,
				      0, 0.001);
    } while (status == GSL_CONTINUE && iter < max_iter);
    
    gsl_root_fsolver_free(s);
  } catch (...) {
    gsl_root_fsolver_free(s);
    throw;
  }
  res.second = r;
  
  return res;
}

double gwas_winners_curse::beta_debiasing_function(double beta_debiased_init, void *params) {
  struct beta_debiasing_params *p = reinterpret_cast<beta_debiasing_params *>(params);
  double beta_debiased = (beta_debiased_init < 1e-16 && beta_debiased_init > -1e-16) ? 1e-16 : beta_debiased_init;
  double f = p->frequency;
  double p_cut = p->p_threshold;
  double beta_biased = p->beta_biased;
  double stderr_biased = p->stderr_biased;
  double beta0 = p->beta0;
  double N = p->N;
  GLM_TYPE type = p->type;
  double precfact = p->precision_factor;
  double vif = p->vif;
  bool invariant_standard_error = p->invariant_standard_error;

  optimized_stderr_calculator stderr_calc;    
  std::pair<double, double> stderr_and_deriv;
  if (invariant_standard_error) {
    stderr_and_deriv.first = stderr_biased;
    stderr_and_deriv.second = 0.0;
  } else {
    stderr_and_deriv = stderr_calc.partial_stderr_expectation(beta_debiased, beta0, f, N, type, precfact, false);
  }
  double stderr = stderr_and_deriv.first * sqrt(vif);
  //std::cout << "beta debiased is " << beta_debiased << std::endl;
  //std::cout << "beta0 is " << beta0 << std::endl;
  //std::cout << "stderr is " << stderr << std::endl;
  //std::cout << "beta biased is " << beta_biased << std::endl;
  //std::cout << "pcut is " << p_cut << std::endl;
  //std::cout << "Q is " << Q << std::endl;
  ////////////////////uncomment me if revert
  //return beta_debiased + stderr * dnorm(Q) / pnorm(Q) - beta_biased;
  if (parameters::get_flag("invalid-approximation")) {
    double J = fabs(beta_debiased) / stderr - qnorm(1.0 - p_cut / 2.0);
    return beta_debiased + stderr * dnorm(J) / pnorm(J) - beta_biased;
  } else {
    double Q = beta_debiased / stderr - qnorm(1.0 - p_cut / 2.0);
    double R = -beta_debiased / stderr - qnorm(1.0 - p_cut / 2.0);
    return beta_debiased + stderr * (dnorm(Q) - dnorm(R)) / (pnorm(Q) + pnorm(R)) - beta_biased;
  }
}

double gwas_winners_curse::beta_debiasing_derivative(double beta_debiased_init, void *params) {
  throw std::domain_error("DO NOT USE THE DEBIASING DERIVATIVE");
  struct beta_debiasing_params *p = reinterpret_cast<beta_debiasing_params *>(params);
  double beta_debiased = (beta_debiased_init < 1e-16 && beta_debiased_init > -1e-16) ? 1e-16 : beta_debiased_init;
  double f = p->frequency;
  double p_cut = p->p_threshold;
  double beta_biased = p->beta_biased;
  double stderr_biased = p->stderr_biased;
  double beta0 = p->beta0;
  double N = p->N;
  GLM_TYPE type = p->type;
  double precfact = p->precision_factor;
  double vif = p->vif;
  bool invariant_standard_error = p->invariant_standard_error;
  if (invariant_standard_error) throw std::domain_error("beta debiasing derivative: invariant standard error not implemented");

  //first the function
  optimized_stderr_calculator stderr_calc;
  std::pair<double, double> stderr_and_deriv = stderr_calc.partial_stderr_expectation(beta_debiased, beta0, f, N, type, precfact, true);
  double stderr = stderr_and_deriv.first * sqrt(vif);
  //std::cout << "beta debiased is " << beta_debiased << std::endl;
  //std::cout << "beta0 is " << beta0 << std::endl;
  //std::cout << "stderr is " << stderr << std::endl;
  //std::cout << "beta biased is " << beta_biased << std::endl;
  //std::cout << "pcut is " << p_cut << std::endl;
  double Q = beta_debiased / stderr - qnorm(1.0 - p_cut / 2.0);
  double R = -beta_debiased / stderr - qnorm(1.0 - p_cut / 2.0);
  double dnormQ = dnorm(Q);
  double dnormR = dnorm(R);
  double pnormQ = pnorm(Q);
  double pnormR = pnorm(R);
  //std::cout << "Q is " << Q << std::endl;
  //*f = beta_debiased + stderr * dnormQ / pnormQ - beta_biased;

  //then its derivative
  //double spartial = beta_debiased / (stderr * (N - 1.0));
  double spartial = stderr_and_deriv.second * sqrt(vif);
  double qpartial = (stderr - beta_debiased * spartial) / (stderr * stderr);
  double rpartial = -qpartial;
  //return 1.0 - stderr * (Q * dnormQ * pnormQ * qpartial - sqrt(2.0) * dnormQ * dnormQ * qpartial) / (pnormQ * pnormQ) + spartial * dnormQ / pnormQ;
  return 1.0 + stderr * ((pnormQ - pnormR) * (Q*dnormQ*qpartial - R*dnormR*rpartial) - (dnormQ - dnormR) * (dnormQ*qpartial - dnormR*rpartial)) / ((pnormQ - pnormR) * (pnormQ - pnormR)) + (dnormQ - dnormR) / (pnormQ - pnormR) * spartial;
}

void gwas_winners_curse::beta_debiasing_bothfunctions(double beta_debiased_init, void *params, double *fxn, double *df) {
  throw std::domain_error("DO NOT USE THE COMBINED DERIVATIVE FUNCTION");
  struct beta_debiasing_params *p = reinterpret_cast<beta_debiasing_params *>(params);
  double beta_debiased = (beta_debiased_init < 1e-16 && beta_debiased_init > -1e-16) ? 1e-16 : beta_debiased_init;
  double f = p->frequency;
  double p_cut = p->p_threshold;
  double beta_biased = p->beta_biased;
  double stderr_biased = p->stderr_biased;
  double beta0 = p->beta0;
  double N = p->N;
  GLM_TYPE type = p->type;
  double precfact = p->precision_factor;
  double vif = p->vif;
  bool invariant_standard_error = p->invariant_standard_error;
  if (invariant_standard_error) throw std::domain_error("beta debiasing bothfunctions: invariant standard error not implemented in derivative");

  //first the function
  optimized_stderr_calculator stderr_calc; 
  std::pair<double, double> stderr_and_deriv = stderr_calc.partial_stderr_expectation(beta_debiased, beta0, f, N, type, precfact, true);
  //std::cout << "beta debiased is " << beta_debiased << std::endl;
  //std::cout << "beta0 is " << beta0 << std::endl;
  //std::cout << "stderr is " << stderr << std::endl;
  //std::cout << "beta biased is " << beta_biased << std::endl;
  //std::cout << "pcut is " << p_cut << std::endl;
  double stderr = stderr_and_deriv.first * sqrt(vif);
  double Q = beta_debiased / stderr - qnorm(1.0 - p_cut / 2.0);
  double R = -beta_debiased / stderr - qnorm(1.0 - p_cut / 2.0);
  double dnormQ = dnorm(Q);
  double dnormR = dnorm(R);
  double pnormQ = pnorm(Q);
  double pnormR = pnorm(R);
  //std::cout << "Q is " << Q << std::endl;
  *fxn = beta_debiased + stderr * (dnormQ - dnormR) / (pnormQ + pnormR) - beta_biased;

  //then its derivative
  //double spartial = beta_debiased / (stderr * (N - 1.0));
  double spartial = stderr_and_deriv.second * sqrt(vif);
  double qpartial = (stderr - beta_debiased * spartial) / (stderr * stderr);
  double rpartial = -qpartial;
  *df = 1.0 + stderr * ((pnormQ - pnormR) * (Q*dnormQ*qpartial - R*dnormR*rpartial) - (dnormQ - dnormR) * (dnormQ*qpartial - dnormR*rpartial)) / ((pnormQ - pnormR) * (pnormQ - pnormR)) + (dnormQ - dnormR) / (pnormQ - pnormR) * spartial;
}


double gwas_winners_curse::zero_finding::debias_beta(gwas_winners_curse::GLM_TYPE type,
							const double &beta_biased,
							const double &stderr_biased,
							const double &frequency,
							const double &N,
							const double &trait_characteristic,
							const double &p_threshold_init,
							const double &precision_factor,
							double &empirical_vif,
							bool use_derivative_mode,
							bool invariant_standard_error) {


  bool verbose = parameters::get_flag("verbose");
  
  double beta_nosign = fabs(beta_biased);


  
  double p_actual = 2.0 * (1.0 - pnorm(beta_nosign / stderr_biased));
  double p_threshold = 0.0;
  if (p_actual > p_threshold_init) {
    p_threshold = p_actual;
    std::ostringstream o;
    o << "WARNING: pair " << beta_biased << "(" << stderr_biased << ") corresponds to p="
      << p_actual << ">" << p_threshold_init << ", adjusting";
    std::cerr << o.str() << std::endl;
    //throw std::domain_error(o.str());
  } else {
    p_threshold = p_threshold_init;
  }


  
  int status;
  int iter = 0, max_iter = 1000;
  const gsl_root_fdfsolver_type *T = 0;
  gsl_root_fdfsolver *s = 0;

  
  const gsl_root_fsolver_type *Talt = 0;
  gsl_root_fsolver *salt = 0;


  double x_lo = -beta_nosign/4.0, x_hi = beta_nosign;
  double x_init = beta_nosign / 2.0;
  double r = x_init, rnext = 0.0;
  gsl_function_fdf F;
  gsl_function Falt;

  double beta0 = calculate_adjusted_trait(type,
					  beta_nosign,
					  frequency,
					  trait_characteristic);

  //do one calibration to estimate effective VIF
  optimized_stderr_calculator stderr_calc;
  double stderr = 0.0;
  if (invariant_standard_error) {
    stderr = stderr_biased;
  } else {
    stderr = stderr_calc.partial_stderr_expectation(beta_biased, beta0, frequency, N, type, precision_factor, false).first;
  }
  double vif = stderr_biased / stderr;
  vif *= vif;
  if (verbose) print_message("effective VIF is " + to_string<double>(vif));

  
  if (vif < 1.0) {
    //std::cout << "WARNING: setting VIF to 1.0" << std::endl;
    //vif = 1.0;
  }

  empirical_vif = vif;
  
  struct gwas_winners_curse::beta_debiasing_params params = {beta_nosign,
								beta0,
								frequency,
								p_threshold,
								vif,
								precision_factor,
								N,
								type,
								stderr_biased,
								invariant_standard_error};


  bool use_mode_fsolver = !use_derivative_mode;

  
  if (verbose) {
    std::ostringstream o;
    o << "function evaluated at boundaries is {" << beta_debiasing_function(x_lo, reinterpret_cast<void *>(&params))
      << "," << beta_debiasing_function(x_hi, reinterpret_cast<void *>(&params)) << "}";
    print_message(o.str());
  }//std::cout << "function evaluated at guess is {" << beta_debiasing_function(x_init, reinterpret_cast<void *>(&params)) << std::endl;
  if (use_mode_fsolver) {
    Falt.function = &gwas_winners_curse::beta_debiasing_function;
    Falt.params = &params;
    Talt = gsl_root_fsolver_brent;
    salt = gsl_root_fsolver_alloc(Talt);
    try {
      gsl_root_fsolver_set(salt, &Falt, x_lo, x_hi);
      do {
	iter++;
	status = gsl_root_fsolver_iterate(salt);
	rnext = gsl_root_fsolver_root(salt);
	x_lo = gsl_root_fsolver_x_lower(salt);
	x_hi = gsl_root_fsolver_x_upper(salt);
	status = gsl_root_test_interval(x_lo, x_hi, 0, 0.001);
	r = rnext;
      } while (status == GSL_CONTINUE && iter < max_iter);
      if (iter == max_iter) {
	throw std::domain_error("fsolver: reached max iteration without convergence");
      }
    } catch (...) {
      gsl_root_fsolver_free(salt);
      throw;
    }
  } else {
    F.f = &gwas_winners_curse::beta_debiasing_function;
    F.params = &params;
    F.df = &gwas_winners_curse::beta_debiasing_derivative;
    F.fdf = &gwas_winners_curse::beta_debiasing_bothfunctions;
    T = gsl_root_fdfsolver_steffenson;
    s = gsl_root_fdfsolver_alloc(T);
    try {
      gsl_root_fdfsolver_set(s, &F, x_init);
      //    std::cout << "starting estimate is " << std::endl;
      do {
	iter++;
	status = gsl_root_fdfsolver_iterate(s);
	rnext = gsl_root_fdfsolver_root(s);
	//x_lo = gsl_root_fsolver_x_lower(s);
	//x_hi = gsl_root_fsolver_x_upper(s);
	status = gsl_root_test_delta(rnext, r,
				     0, 0.001);
	r = rnext;
	//std::cout << "updated estimate is " << r << std::endl;
      } while (status == GSL_CONTINUE && iter < max_iter);
      if (iter == max_iter) {
	//throw std::domain_error("fdfsolver: reached max iteration without convergence");
	std::cerr << "fdfsolver warning: no convergence" << std::endl;
      }
      gsl_root_fdfsolver_free(s);
    } catch (...) {
      gsl_root_fdfsolver_free(s);
      throw;
    }
  }

  double beta_mle = beta_biased >= 0.0 ? r : -r;
  return beta_mle;
}




double gwas_winners_curse::zero_finding::calculate_adjusted_trait(gwas_winners_curse::GLM_TYPE type,
								     const double &beta,
								     const double &freq,
								     const double &trait_characteristic) const {

  if (type == logistic) {
    //predicted logistic additive intercept beta is not closed; choose positive real root
    
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T = 0;
    gsl_root_fsolver *s = 0;
    double x_lo = 0.0, x_hi = 1.0;
    double r = 0.0;
    gsl_function F;
    
    struct gwas_winners_curse::binomial_root_params params = {trait_characteristic, beta, freq};
    
    if (beta < 0.0) {
      x_lo = 0.0;
      x_hi = (1.0 - trait_characteristic) / trait_characteristic;
    } else {
      x_lo = (1.0 - trait_characteristic) / trait_characteristic;
      x_hi = trait_characteristic / 8.0;
      x_hi = (1.0 - x_hi) / x_hi;
    }
    
    F.function = &gwas_winners_curse::binomial_root;
    F.params = &params;
    
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);
    try {
      gsl_root_fsolver_set(s, &F, x_lo, x_hi);
      
      do {
	iter++;
	status = gsl_root_fsolver_iterate(s);
	r = gsl_root_fsolver_root(s);
	x_lo = gsl_root_fsolver_x_lower(s);
	x_hi = gsl_root_fsolver_x_upper(s);
	status = gsl_root_test_interval(x_lo, x_hi,
					0, 0.001);
      } while (status == GSL_CONTINUE && iter < max_iter);
      
      gsl_root_fsolver_free(s);
    } catch (...) {
      gsl_root_fsolver_free(s);
      throw;
    }
    
    //r = 1.0 / (1.0 + r);
    return -log(r);
  } else if (type == linear) {
    //adjust mean trait by per-snp contribution
    return trait_characteristic - 2.0 * freq * beta;
  } else if (type == poisson) {
    //beta0 is log(count) for all other predictors 0; thus can take difference of class logs to calculate
    return log(trait_characteristic) - log(freq*freq*exp(2.0*beta) + 2.0*freq*(1.0-freq)*exp(beta) + (1.0-freq)*(1.0-freq));
  } else {
    throw std::domain_error("gwas_winners_curse::zero_finding::calculate_adjusted_trait: unrecognized GLM type");
  }
}
