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

#ifndef __GWAS_WINNERS_CURSE_UTILITIES_H__
#define __GWAS_WINNERS_CURSE_UTILITIES_H__

#include <sstream>
#include <string>
#include <iostream>
#include <stdexcept>
#include <cmath>

#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"

namespace gwas_winners_curse {

  template <class value_type>
    std::string to_string(const value_type &obj) {
    std::ostringstream o;
    if (!(o << obj))
      throw std::domain_error("to_string: conversion failed");
    return o.str();
  }

  template <class value_type>
    value_type from_string(const std::string &str) {
    std::istringstream strm1(str);
    value_type res;
    if (!(strm1 >> res))
      throw std::domain_error("from_string: conversion failed: \"" + str + "\"");
    return res;
  }

  inline void print_message(const std::string &str) {
    std::clog << str << std::endl;
  }
  
  typedef enum {
    linear, logistic, poisson
  } GLM_TYPE;
  
  inline unsigned max(unsigned u1, unsigned u2) throw() {
    return u1 > u2 ? u1 : u2;
  }

  inline double max(const double &d1, const double &d2) throw() {
    return d1 > d2 ? d1 : d2;
  }

  inline unsigned convert_to_tabs(unsigned filled_chars, unsigned max_tabs, unsigned space_per_tab) {
    unsigned filled_tabs = filled_chars / space_per_tab;
    return max_tabs - filled_tabs;
  }

  inline bool cicompare(const std::string &str1, const std::string &str2) {
    if (str1.size() == str2.size()) {
      for (unsigned i = 0; i < str1.size(); ++i) {
	if (tolower(str1.at(i)) != tolower(str2.at(i))) return false;
      }
      return true;
    }
    return false;
  }

  inline double dnorm(const double &x) {
    return gsl_ran_ugaussian_pdf(x);
  }

  inline double pnorm(const double &x) {
    return gsl_cdf_ugaussian_P(x);
  }

  inline double qnorm(const double &p) {
    return gsl_cdf_ugaussian_Pinv(p);
  }

  inline double dchisq(const double &x) {
    return gsl_ran_chisq_pdf(x, 1.0);
  }

  inline double pchisq(const double &x, bool lower_tail = true) {
    if (lower_tail) 
      return gsl_cdf_chisq_P(x, 1.0);
    else
      return gsl_cdf_chisq_Q(x, 1.0);
  }
  
  inline double qchisq(const double &p, bool lower_tail = true) {
    if (lower_tail)
      return gsl_cdf_chisq_Pinv(p, 1.0);
    else
      return gsl_cdf_chisq_Qinv(p, 1.0);
  }
}

#endif //__GWAS_WINNERS_CURSE_UTILITIES_H__
