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

#ifndef __GWAS_WINNERS_CURSE_VARIANCE_WEIGHTS_H__
#define __GWAS_WINNERS_CURSE_VARIANCE_WEIGHTS_H__

#include <cmath>

namespace gwas_winners_curse {
  inline double linear_variance_weight(const double &beta0,
				       const double &beta1,
				       unsigned genotype) {
    return 1.0;
  }

  inline double linear_variance_derivative(const double &beta0,
					   const double &beta1,
					   unsigned genotype) {
    return 0.0;
  }
  inline double logistic_variance_weight(const double &beta0,
					 const double &beta1,
					 unsigned genotype) {
    double expXB = exp(beta0 + beta1 * genotype);
    return expXB/pow(1.0 + expXB, 2);
  }
  inline double logistic_variance_derivative(const double &beta0,
					     const double &beta1,
					     unsigned genotype) {
    double expXB = exp(beta0 + beta1 * genotype);
    return (genotype * expXB * (1.0 - expXB)) / (expXB * expXB * expXB);
  }
  inline double poisson_variance_weight(const double &beta0,
					const double &beta1,
					unsigned genotype) {
    return exp(beta0 + beta1 * genotype);
  }
  inline double poisson_variance_derivative(const double &beta0,
					    const double &beta1,
					    unsigned genotype) {
    return genotype * exp(beta0 + beta1 * genotype);
  }
}

#endif //__GWAS_WINNERS_CURSE_VARIANCE_WEIGHTS_H__
