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

#ifndef __GWAS_WINNERS_CURSE_MULTINOMIAL_H__
#define __GWAS_WINNERS_CURSE_MULTINOMIAL_H__

#include <vector>
#include <iostream>
#include <algorithm>
#include <utility>
#include <cstdio>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "gwas_winners_curse/utilities.h"
#include "gwas_winners_curse/variance_weights.h"

namespace gwas_winners_curse {
  inline double get_multinomial_probability(const std::vector<unsigned> &counts,
					    const std::vector<double> &probabilities) {
    if (counts.size() != probabilities.size() || counts.empty()) {
      throw std::domain_error("get_multinomial_probability: invalid dimensions: \"" +
			      to_string<unsigned>(counts.size()) + "\" \"" +
			      to_string<unsigned>(probabilities.size()) + "\"");
    }
    return gsl_ran_multinomial_pdf(counts.size(), probabilities.data(), counts.data());
  }

  class optimized_stderr_calculator {
  public:
    optimized_stderr_calculator() {}
    ~optimized_stderr_calculator() throw() {}

    double complete_stderr_expectation(const double &beta,
				       const double &beta0,
				       const double &freq,
				       unsigned N,
				       GLM_TYPE type) const;
    
    std::pair<double, double> partial_stderr_expectation(const double &beta,
							 const double &beta0,
							 const double &freq,
							 const double &N,
							 GLM_TYPE type,
							 const double &precision_adj_factor,
							 bool actually_compute_derivative) const;
  private:
    double stderr_expectation_term(const double &freq,
				   const double &multinom_prob,
				   const double &w0,
				   const double &w1,
				   const double &w2,
				   unsigned i,
				   unsigned j,
				   unsigned k) const;
    double stderr_expectation_derivative_term(const double &freq,
					      const double &multinom_prob,
					      const double &w0,
					      const double &w1,
					      const double &w2,
					      const double &dw0,
					      const double &dw1,
					      const double &dw2,
					      unsigned i,
					      unsigned j,
					      unsigned k) const;
  };
}

#endif //__GWAS_WINNERS_CURSE_MULTINOMIAL_H__
