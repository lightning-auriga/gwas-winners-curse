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

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <utility>
#include <stdexcept>
#include <cmath>
#include <cfloat>
#include "gwas_winners_curse/binomial_zero_finding.h"
#include "gwas_winners_curse/cargs.h"
#include "gwas_winners_curse/multinomial.h"
#include "gwas_winners_curse/parameters.h"
#include "gwas_winners_curse/utilities.h"
#include "gwas_winners_curse/zero_finding.h"
#include "fileinterface/fileinterface.h"

using namespace gwas_winners_curse;

int main(int argc, char **argv) {
  fileinterface_reader *input = 0;
  fileinterface_writer *output = 0;
  std::string line = "", catcher = "";
  double beta_disc = 0.0, stderr_disc = 0.0, freq_disc = 0.0,
    beta_repl = 0.0, stderr_repl = 0.0, freq_repl = 0.0,
    p_thresh = 0.0;
  double n_disc = 0.0, n_repl = 0.0;

  zero_finding solver;
  double debiased_beta = 0.0;

  bool verbose = false;
  
  try {
    cargs arg_parser(argc, argv);
    arg_parser.parse_args(true);

    if (argc == 1 ||
	parameters::get_flag("help")) {
      parameters::print_help(std::cerr);
      return 0;
    }

    verbose = parameters::get_flag("verbose");
    
    std::string user_type = parameters::get_parameter("regression-type");
    GLM_TYPE type;
    if (cicompare(user_type, "linear") ||
	cicompare(user_type, "normal")) {
      type = linear;
      if (verbose) print_message("detected linear regression");
    } else if (cicompare(user_type, "logistic") ||
	       cicompare(user_type, "binomial") ||
	       cicompare(user_type, "logit")) {
      type = logistic;
      if (verbose) print_message("detected logistic regression");
    } else if (cicompare(user_type, "poisson") ||
	       cicompare(user_type, "log")) {
      type = poisson;
      if (verbose) print_message("detected poisson regression");
    } else {
      throw std::domain_error("cannot parse regression type from user-specified \"" + user_type + "\"");
    }

    std::string input_filename = parameters::get_parameter("input-file");
    std::string output_filename = parameters::get_parameter("output-file");

    

    double stderr_precision_factor = from_string<double>(parameters::get_parameter("stderr-precision-factor"));


    bool per_line_mean = false;
    double trait_mean = 0.0;
    if (parameters::get_parameter("trait-mean").empty()) {
      per_line_mean = true;
    } else {
      trait_mean = from_string<double>(parameters::get_parameter("trait-mean"));
    }
    if (parameters::get_parameter("discovery-threshold").empty()) {
      if (!per_line_mean) {
	throw std::domain_error("with a global mean specified, the global discovery threshold must also be specified");
      }
    } else {
      p_thresh = from_string<double>(parameters::get_parameter("discovery-threshold"));      
      if (p_thresh <= 0.0 || p_thresh >= 1.0)
	throw std::domain_error("invalid discovery p-value threshold \"" + parameters::get_parameter("discovery-threshold") + "\"");
    }


    
    if (parameters::get_flag("scan-beta")) {
      if (per_line_mean)
	throw std::domain_error("in scan beta mode, a single shared trait mean must be specified");
      double freq_init = 0.01, freq_interval = 0.01, freq_max = 0.5;
      int N_init = 1000, N_interval = 1000, N_max = 100000;
      double beta_init = 0.01, beta_interval = 0.01, beta_max = 1.0;
      optimized_stderr_calculator stderr_calc;
      //iterate over frequencies
      for (double freq = freq_init; freq <= freq_max; freq += freq_interval) {
	//iterate over sample sizes
	std::cout << "starting frequency " << freq << std::endl;
	double abort_N = false;
	for (unsigned N = N_max; N >= N_init && !abort_N; N -= N_interval) {
	  //iterate over betas
	  for (double beta = beta_init; beta <= beta_max; beta += beta_interval) {
	    //std::cerr << "computing beta0" << std::endl;
	    double beta0 = solver.calculate_adjusted_trait(type,
							   fabs(beta),
							   freq,
							   trait_mean);
	    std::cerr << "beta0 is " << beta0 << std::endl;
	    //std::cerr << "computing stderr" << std::endl;
	    std::pair<double, double> stderr_not_deriv = stderr_calc.partial_stderr_expectation(beta, beta0,
												freq, N,
												type, stderr_precision_factor, false);
	    //std::cerr << "stderr is " << stderr_not_deriv.first << std::endl;
	    double expected_p = 0.0;
	    if ((expected_p = pchisq(beta * beta / (stderr_not_deriv.first * stderr_not_deriv.first), false)) < (p_thresh * 1.0e-9)) {
	      break;
	    }
	    if (expected_p > p_thresh) {
	      //std::cout << "completed a beta loop (beta too large)" << std::endl;
	      if (beta + beta_interval > beta_max) {
		std::cerr << "warning: last beta in series is too small (p=" << expected_p << "), all smaller N will also be too small so skipping" << std::endl;
		abort_N = true;
	      }

	      continue;
	    }
	    //valid pair
	    std::cout << beta << ' ' << stderr_not_deriv.first << ' ' << N << ' ' << freq << ' ' << trait_mean << ' ' << p_thresh << ' ' << expected_p << std::endl;
	  }
	  std::cerr << "completed a beta loop " << N << ' ' << freq << std::endl;
	}
	std::cerr << "completed an N loop" << std::endl;
      }
      std::cerr << "completed freq loop" << std::endl;
      return 0;
    }
    
    if (verbose) print_message("opening files");
    input = reconcile_reader(input_filename);
    output = reconcile_writer(output_filename);

    //get line count the first pass
    unsigned linecount = 0;
    while (input->getline(line)) {
      ++linecount;
    }
    input->close();
    delete input;
    input = reconcile_reader(input_filename);

    if (parameters::get_flag("input-has-header")) input->getline(line);



    while (input->getline(line)) {


      std::istringstream strm1(line);


      if (parameters::get_flag("estimate-stderr")) {
	if (!(strm1 >> beta_disc >> stderr_disc >> n_disc >> freq_disc))
	  throw std::domain_error("could not parse input line \"" + line + "\"");
	if (per_line_mean) {
	  if (!(strm1 >> trait_mean >> p_thresh))
	    throw std::domain_error("could not parse trait mean/discovery threshold from input line \"" + line + "\"");
	}
	optimized_stderr_calculator stderr_calc;
	if (verbose) print_message("repeating beta0 computation");
	double beta0 = solver.calculate_adjusted_trait(type,
						       fabs(beta_disc),
						       freq_disc,
						       trait_mean);
	if (verbose) print_message("repeating stderr computation");
	std::pair<double, double> stderr_not_deriv = stderr_calc.partial_stderr_expectation(beta_disc, beta0,
											    freq_disc, n_disc,
											    type, stderr_precision_factor, false);
	std::cout << beta_disc << ' ' << stderr_disc << ' ' << n_disc << ' ' << freq_disc << ' ' << stderr_not_deriv.first << ' ' << pow(stderr_disc / stderr_not_deriv.first,2) << std::endl;
	continue;
      }




      

      if (!(strm1 >> beta_disc >> stderr_disc >> n_disc >> freq_disc
	    >> beta_repl >> stderr_repl >> n_repl >> freq_repl)) {
	throw std::domain_error("could not parse input line \"" + line + "\"");
      }
      if (per_line_mean) {
	if (!(strm1 >> trait_mean >> p_thresh))
	  throw std::domain_error("could not parse trait mean/discovery threshold from input line \"" + line + "\"");
      }
      if (stderr_disc <= 0.0 || stderr_repl <= 0.0 ||
	  freq_disc <= 0.0 || freq_disc >= 1.0 ||
	  pchisq(pow(beta_disc / stderr_disc, 2), false) > p_thresh*from_string<double>(parameters::get_parameter("threshold-precision-factor"))) {
	throw std::domain_error("something is horribly wrong with file \"" + input_filename + "\" at this line: \"" + line + "\"");
      }
      if (verbose) {
	std::ostringstream o;
	o << "input data are " << beta_disc << ' ' << stderr_disc << ' ' << n_disc << ' ' << freq_disc;
	print_message(o.str());
      }
      double vif = 0.0;
      try {
	if (verbose) print_message("debiasing");
	//if (beta_disc < beta_repl) {
	//std::cout << "skipping variant with replication stronger than discovery" << std::endl;
	//continue;
	//}
	debiased_beta = solver.debias_beta(type,
					   beta_disc,
					   stderr_disc,
					   freq_disc,
					   n_disc,
					   trait_mean,
					   p_thresh,
					   stderr_precision_factor,
					   vif,
					   false,
					   parameters::get_flag("invariant-standard-error"));
      } catch (const std::domain_error &e) {
	std::cerr << "from debiasing: " << e.what() << std::endl;
	continue;
      } catch (...) {
	throw;
      }
      optimized_stderr_calculator stderr_calc;
      if (verbose) print_message("repeating beta0 computation");
      double beta0 = solver.calculate_adjusted_trait(type,
						     fabs(beta_disc),
						     freq_disc,
						     trait_mean);
      if (verbose) print_message("repeating stderr computation");
      std::pair<double, double> stderr_not_deriv;

      if (parameters::get_flag("invariant-standard-error")) {
	stderr_not_deriv.first = stderr_disc;
	stderr_not_deriv.second = 0.0;
      } else {
	stderr_not_deriv = stderr_calc.partial_stderr_expectation(debiased_beta, beta0,
								  freq_disc, n_disc,
								  type, stderr_precision_factor, false);
	stderr_not_deriv.first *= sqrt(vif);
      }
      std::pair<double, double> ci;
      if (verbose) print_message("computing confidence interval");
      ci = solver.calculate_ci(debiased_beta, stderr_not_deriv.first, p_thresh);
      

      if (vif > 20 && trait_mean < 20) throw std::domain_error("possible serious vif error for \"" + input_filename + "\"");

      //the above debiasing corresponds the beta_mle from original paper
      //compute beta_mse and ci_mse now
      double K = stderr_disc*stderr_disc / (stderr_disc*stderr_disc + pow(beta_disc - debiased_beta,2));
      double K_lower = stderr_disc*stderr_disc / (stderr_disc*stderr_disc + pow(beta_disc-1.96*stderr_disc - ci.first,2));
      double K_upper = stderr_disc*stderr_disc / (stderr_disc*stderr_disc + pow(beta_disc-1.96*stderr_disc - ci.second,2));


      K = K_lower = K_upper = 0.0;
      
      
      double beta_mse = K*beta_disc + (1.0-K)*debiased_beta;
      std::pair<double, double> ci_mse;
      ci_mse.first = K_lower*(beta_disc-1.96*stderr_disc) + (1.0-K_lower)*ci.first;
      ci_mse.second = K_upper*(beta_disc+1.96*stderr_disc) + (1.0-K_upper)*ci.second;
      if (parameters::get_flag("invariant-standard-error")) {
	stderr_not_deriv.first = stderr_disc;
	stderr_not_deriv.second = 0.0;
      } else {
	stderr_not_deriv = stderr_calc.partial_stderr_expectation(beta_mse, beta0,
								  freq_disc, n_disc,
								  type, stderr_precision_factor, false);
	stderr_not_deriv.first *= sqrt(vif);
      }
      std::cout << "\tand resulting debiased beta is " << beta_mse << "[" << ci_mse.first << ',' << ci_mse.second
		<< "] (original was " << beta_disc << ", repl was " << beta_repl << ")\n";
      std::ostringstream o;
      
      std::istringstream strm_repair(line);
      for (unsigned xi = 0; xi < 10; ++xi) {
	strm_repair >> catcher;
	if (xi) o << ' ';
	o << catcher;
      }
      strm_repair >> catcher >> catcher >> catcher >> catcher;
      o << ' ' << catcher;
      o << ' ' << beta_mse << ' ' << stderr_not_deriv.first << ' ' << ci_mse.first << ' ' << ci_mse.second << ' ' << linecount << ' ' << p_thresh;
      output->writeline(o.str());
    }

    input->close();
    delete input;
    input = 0;
    output->close();
    delete output;
    output = 0;
    return 0;
  } catch (const std::domain_error &e) {
    if (input) delete input;
    if (output) delete output;
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  } catch (const std::bad_alloc &e) {
    if (input) delete input;
    if (output) delete output;
    std::cerr << "error: out of memory" << std::endl;
    return 2;
  } catch (...) {
    if (input) delete input;
    if (output) delete output;
    std::cerr << "error: unhandled exception" << std::endl;
    return 3;
  }
}
