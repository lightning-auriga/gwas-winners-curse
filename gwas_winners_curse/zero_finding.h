/*
  Copyright 2016 Lightning Auriga

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

#ifndef __GWAS_WINNERS_CURSE_ZERO_FINDING_H__
#define __GWAS_WINNERS_CURSE_ZERO_FINDING_H__

#include <stdexcept>
#include <cstdio>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "gwas_winners_curse/binomial_zero_finding.h"
#include "gwas_winners_curse/multinomial.h"
#include "gwas_winners_curse/parameters.h"
#include "gwas_winners_curse/utilities.h"

namespace gwas_winners_curse {
  struct beta_debiasing_params {
    double beta_biased, beta0, frequency, p_threshold, vif, precision_factor;
    double N;
    GLM_TYPE type;
    double stderr_biased;
    bool invariant_standard_error;
  };

  struct ci_log_likelihood_params {
    double beta_debiased, stderr_debiased, p_threshold;
  };
  
  double beta_debiasing_function(double beta_debiased, void *params);
  double beta_debiasing_derivative(double beta_debiased, void *params);
  void beta_debiasing_bothfunctions(double beta_debiased, void *params, double *f, double *df);

  double ci_log_likelihood(double x, void *params);
  
  class zero_finding {
  public:
    zero_finding() {}
    ~zero_finding() throw() {}

    double calculate_adjusted_trait(GLM_TYPE type,
				    const double &beta,
				    const double &freq,
				    const double &trait_characteristic) const;
    
    double debias_beta(GLM_TYPE type,
		       const double &beta_biased,
		       const double &stderr_biased,
		       const double &frequency,
		       const double &N,
		       const double &trait_characteristic,
		       const double &p_threshold,
		       const double &precision_factor,
		       double &vif,
		       bool use_derivative_mode,
		       bool invariant_standard_error);

    std::pair<double, double> calculate_ci(const double &beta_debiased,
					   const double &stderr_debiased,
					   const double &p_threshold);
  private:
  };
}

#endif //__GWAS_WINNERS_CURSE_ZERO_FINDING_H__
