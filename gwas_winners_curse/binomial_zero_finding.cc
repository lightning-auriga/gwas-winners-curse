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

#include "gwas_winners_curse/binomial_zero_finding.h"

double gwas_winners_curse::binomial_root(double y, void *params) {
  struct binomial_root_params *p = reinterpret_cast<binomial_root_params *>(params);

  double x = p->x;
  double beta = p->beta;
  double f = p->f;

  double expbeta = exp(-beta);
  double exp2beta = exp(-2.0*beta);
  double exp3beta = exp(-3.0*beta);

  double hprob = f*f;
  double hetprob = 2.0*f*(1.0-f);
  double otherhprob = (1.0-f)*(1.0-f);
  
  double constant_coeff = x - 1.0;
  double y_coeff = x*(1.0 + expbeta + exp2beta) - hprob*(1.0 + expbeta) - hetprob*(1.0 + exp2beta) - otherhprob*(expbeta + exp2beta);
  double y2_coeff = x*(expbeta + exp2beta + exp3beta) - hprob*expbeta - hetprob*exp2beta - otherhprob*exp3beta;
  double y3_coeff = x*exp3beta;

  return constant_coeff + y_coeff * y + y2_coeff * y * y + y3_coeff * y * y * y;
}
