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

#ifndef __GWAS_WINNERS_CURSE_BINOMIAL_ZERO_FINDING_H__
#define __GWAS_WINNERS_CURSE_BINOMIAL_ZERO_FINDING_H__

#include <iostream>
#include <cmath>

namespace gwas_winners_curse {
  struct binomial_root_params {
    double x, beta, f;
  };

  double binomial_root(double y, void *params);
}

#endif //__GWAS_WINNERS_CURSE_BINOMIAL_ZERO_FINDING_H__
