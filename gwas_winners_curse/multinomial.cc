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

#include "gwas_winners_curse/multinomial.h"

double gwas_winners_curse::optimized_stderr_calculator::complete_stderr_expectation(const double &beta,
										       const double &beta0,
										       const double &freq,
										       unsigned N,
										       GLM_TYPE type) const {
  double total = 0.0;
  double w0 = 0.0, w1 = 0.0, w2 = 0.0, multinom_prob = 0.0;
  std::vector<unsigned> counts(3, 0);
  std::vector<double> probabilities;
  probabilities.push_back((1.0 - freq) * (1.0 - freq));
  probabilities.push_back(2.0 * freq * (1.0 - freq));
  probabilities.push_back(freq * freq);

  //set variance weights based on current GLM linker
  if (type == linear) {
    w0 = linear_variance_weight(beta0, beta, 0);
    w1 = linear_variance_weight(beta0, beta, 1);
    w2 = linear_variance_weight(beta0, beta, 2);
  } else if (type == logistic) {
    w0 = logistic_variance_weight(beta0, beta, 0);
    w1 = logistic_variance_weight(beta0, beta, 1);
    w2 = logistic_variance_weight(beta0, beta, 2);
  } else if (type == poisson) {
    w0 = poisson_variance_weight(beta0, beta, 0);
    w1 = poisson_variance_weight(beta0, beta, 1);
    w2 = poisson_variance_weight(beta0, beta, 2);
  } else {
    throw std::domain_error("optimized_stderr_calculator::complete_stderr_expectation: unrecognized regression type");
  }


  unsigned n_skipped = 0;
  
  for ( ; counts.at(0) <= (N-1); ++counts.at(0)) {
    for (counts.at(1) = 0; counts.at(1) <= (N-counts.at(0)-1); ++counts.at(1)) {
      counts.at(2) = N - counts.at(0) - counts.at(1) - 1;
      if ((!counts.at(0) && !counts.at(1)) ||
	  (!counts.at(0) && !counts.at(2)) ||
	  (!counts.at(1) && !counts.at(2))) continue;

      multinom_prob = get_multinomial_probability(counts, probabilities);
      if (multinom_prob < 0) {
	++n_skipped;
	continue;
      }
      total += stderr_expectation_term(freq,
				       multinom_prob,
				       w0,
				       w1,
				       w2,
				       counts.at(0),
				       counts.at(1),
				       counts.at(2));
    }
  }
  //std::cout << "calculation complete, skipped " << n_skipped << std::endl;
  return sqrt(N * total);
}

std::pair<double, double> gwas_winners_curse::optimized_stderr_calculator::partial_stderr_expectation(const double &beta,
													 const double &beta0,
													 const double &freq,
													 const double &N_init,
													 GLM_TYPE type,
													 const double &precision_adj_factor,
													 bool actually_compute_derivative) const {
  unsigned N = static_cast<unsigned>(N_init);
  std::vector<double> all_points, all_deriv_points;
  double total = 0.0, deriv_total = 0.0, sorted_total = 0.0, sorted_deriv_total = 0.0;
  double w0 = 0.0, w1 = 0.0, w2 = 0.0, multinom_prob = 0.0, value = 0.0;
  double dw0 = 0.0, dw1 = 0.0, dw2 = 0.0;
  std::vector<unsigned> counts(3, 0);
  std::vector<double> probabilities;
  probabilities.push_back((1.0 - freq) * (1.0 - freq));
  probabilities.push_back(2.0 * freq * (1.0 - freq));
  probabilities.push_back(freq * freq);

  //set variance weights based on current GLM linker
  if (type == linear) {
    w0 = linear_variance_weight(beta0, beta, 0);
    w1 = linear_variance_weight(beta0, beta, 1);
    w2 = linear_variance_weight(beta0, beta, 2);
    if (actually_compute_derivative) {
      dw0 = linear_variance_derivative(beta0, beta, 0);
      dw1 = linear_variance_derivative(beta0, beta, 1);
      dw2 = linear_variance_derivative(beta0, beta, 2);
    }
  } else if (type == logistic) {
    w0 = logistic_variance_weight(beta0, beta, 0);
    w1 = logistic_variance_weight(beta0, beta, 1);
    w2 = logistic_variance_weight(beta0, beta, 2);
    if (actually_compute_derivative) {
      dw0 = logistic_variance_derivative(beta0, beta, 0);
      dw1 = logistic_variance_derivative(beta0, beta, 1);
      dw2 = logistic_variance_derivative(beta0, beta, 2);
    }
  } else if (type == poisson) {
    w0 = poisson_variance_weight(beta0, beta, 0);
    w1 = poisson_variance_weight(beta0, beta, 1);
    w2 = poisson_variance_weight(beta0, beta, 2);
    if (actually_compute_derivative) {
      dw0 = poisson_variance_derivative(beta0, beta, 0);
      dw1 = poisson_variance_derivative(beta0, beta, 1);
      dw2 = poisson_variance_derivative(beta0, beta, 2);
    }
  } else {
    throw std::domain_error("optimized_stderr_calculator::partial_stderr_expectation: unrecognized regression type");
  }


  unsigned n_accepted = 0;
  long long n_total = (N-1)*(N-2)/2;

  double min_accepted_probability = 1.0 / n_total / precision_adj_factor;
  
  for ( ; counts.at(0) <= (N-1); ++counts.at(0)) {

    //start at the most likely j
    unsigned most_likely_j_index = (N - counts.at(0) - 1.0) * probabilities.at(1) / (probabilities.at(1) + probabilities.at(2));

    for (counts.at(1) = most_likely_j_index + 1; counts.at(1) <= (N-counts.at(0)-1); ++counts.at(1)) {
      counts.at(2) = N - counts.at(0) - counts.at(1) - 1;
      if ((!counts.at(0) && !counts.at(1)) ||
	  (!counts.at(0) && !counts.at(2)) ||
	  (!counts.at(1) && !counts.at(2))) continue;

      multinom_prob = get_multinomial_probability(counts, probabilities);
      if (multinom_prob < min_accepted_probability) break;
      ++n_accepted;
      value = stderr_expectation_term(freq,
				      multinom_prob,
				      w0,
				      w1,
				      w2,
				      counts.at(0),
				      counts.at(1),
				      counts.at(2));
      total += value;
      all_points.push_back(value);

      if (actually_compute_derivative) {
	value = stderr_expectation_derivative_term(freq,
						   multinom_prob,
						   w0,
						   w1,
						   w2,
						   dw0,
						   dw1,
						   dw2,
						   counts.at(0),
						   counts.at(1),
						   counts.at(2));
	deriv_total += value;
	all_deriv_points.push_back(value);
      }
    }

    int jcount = most_likely_j_index;
    for ( ; jcount >= 0; --jcount) {
      counts.at(1) = static_cast<unsigned>(jcount);
      counts.at(2) = N - counts.at(0) - counts.at(1) - 1;
      if ((!counts.at(0) && !counts.at(1)) ||
	  (!counts.at(0) && !counts.at(2)) ||
	  (!counts.at(1) && !counts.at(2))) continue;

      multinom_prob = get_multinomial_probability(counts, probabilities);
      if (multinom_prob < min_accepted_probability) break;
      ++n_accepted;
      value = stderr_expectation_term(freq,
				      multinom_prob,
				      w0,
				      w1,
				      w2,
				      counts.at(0),
				      counts.at(1),
				      counts.at(2));
      total += value;
      all_points.push_back(value);

      if (actually_compute_derivative) {
	value = stderr_expectation_derivative_term(freq,
						   multinom_prob,
						   w0,
						   w1,
						   w2,
						   dw0,
						   dw1,
						   dw2,
						   counts.at(0),
						   counts.at(1),
						   counts.at(2));
	deriv_total += value;
	all_deriv_points.push_back(value);
      }
    }
  }
  std::sort(all_points.begin(), all_points.end());
  std::sort(all_deriv_points.begin(), all_deriv_points.end());
  //std::cout << "calculation complete, skipped " << (n_total - n_accepted) << " of " << n_total << "; result was " << sqrt(N * total) << "/";

  for (std::vector<double>::const_iterator iter = all_points.begin(); iter != all_points.end(); ++iter) {
    sorted_total += *iter;
  }
  for (std::vector<double>::const_iterator iter = all_deriv_points.begin(); iter != all_deriv_points.end(); ++iter) {
    sorted_deriv_total += *iter;
  }
  //std::cout << sqrt(N * sorted_total) << std::endl;
  return std::pair<double, double>(sqrt(N * sorted_total), sqrt(N * sorted_deriv_total));
}


double gwas_winners_curse::optimized_stderr_calculator::stderr_expectation_term(const double &freq,
										   const double &multinom_prob,
										   const double &w0,
										   const double &w1,
										   const double &w2,
										   unsigned i,
										   unsigned j,
										   unsigned k) const {
  double first_term = w0 * (1.0 - freq) * (1.0 - freq) * multinom_prob / ((i + 1.0) * j * w0 * w1 + j * k * w1 * w2 + 4.0 * (i + 1.0) * k * w0 * w2);
  double second_term = 2.0 * w1 * freq * (1.0 - freq) * multinom_prob / (i * (j + 1.0) * w0 * w1 + (j + 1.0) * k * w1 * w2 + 4.0 * i * k * w0 * w2);
  double third_term = w2 * freq * freq * multinom_prob / (i * j * w0 * w1 + j * (k + 1.0) * w1 * w2 + 4.0 * i * (k + 1.0) * w0 * w2);
  return first_term + second_term + third_term;
}


double gwas_winners_curse::optimized_stderr_calculator::stderr_expectation_derivative_term(const double &freq,
											      const double &multinom_prob,
											      const double &w0,
											      const double &w1,
											      const double &w2,
											      const double &dw0,
											      const double &dw1,
											      const double &dw2,
											      unsigned i,
											      unsigned j,
											      unsigned k) const {
  double first_term_denom = (i+1.0)*j*w0*w1 + j*k*w1*w2 + 4.0*(i+1.0)*k*w0*w2;
  double first_term = first_term_denom * dw0*(1.0-freq)*(1.0-freq)*multinom_prob -
    w0*(1.0-freq)*(1.0-freq)*multinom_prob * ((i+1.0)*j*(w0*dw1+dw0*w1) + j*k*(dw1*w2+w1*dw2) + 4.0*(i+1.0)*k*(dw0*w2+w0*dw2));


  double second_term_denom = i*(j+1.0)*w0*w1 + (j+1.0)*k*w1*w2 + 4.0*i*k*w0*w2;
  double second_term = second_term_denom * 2.0*dw1*freq*(1.0-freq)*multinom_prob -
    2.0*w1*freq*(1.0-freq)*multinom_prob * (i*(j+1.0)*(w0*dw1+dw0*w1) + (j+1.0)*k*(dw1*w2+w1*dw2) + 4.0*i*k*(dw0*w2+w0*dw2));


  double third_term_denom = i*j*w0*w1 + j*(k+1.0)*w1*w2 + 4.0*i*(k+1.0)*w0*w2;
  double third_term = third_term_denom * dw2*freq*freq*multinom_prob -
    w2*freq*freq*multinom_prob * (i*j*(w0*dw1+dw0*w1) + j*(k+1.0)*(dw1*w2+w1*dw2) + 4.0*i*(k+1.0)*(dw0*w2+w0*dw2));
  
  
  return first_term / (first_term_denom * first_term_denom) +
    second_term / (second_term_denom * second_term_denom) +
    third_term / (third_term_denom * third_term_denom);
}
