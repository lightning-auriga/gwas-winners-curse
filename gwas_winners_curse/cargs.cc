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

#include "gwas_winners_curse/cargs.h"

void gwas_winners_curse::cargs::parse_args(bool warn_unused) throw (std::domain_error) {
  //check for each accepted flag

  //set_flag("license", "l", "Print the license to this program", false);
  set_flag("version", "v", "Print the program version and license", false);
  set_flag("help", "h", "Print this help information", false);
  set_flag("verbose", "vb", "Print additional log information", false);
  set_flag("scan-beta", "sb", "Find valid parameter combinations for simulation", false);
  set_flag("estimate-stderr", "es", "Just compute stderr from first four entries of file", false);
  set_flag("input-has-header", "ihh", "Whether input file has skippable header", false);
  set_flag("invariant-standard-error", "ise", "Whether to fix the standard error at the discovery estimate", false);
  set_flag("invalid-approximation", "ia", "Do not use", false);
  //set_flag("use-derivative", "ud", "Use derivative for beta zero finding", false);
  
  set_parameter("input-file",
		"if",
		"Name of file containing biased variant data",
		"");
  set_parameter("regression-type",
		"rt",
		"linear, logistic, poisson",
		"linear");
  set_parameter("output-file",
		"of",
		"Name of file to which to write results",
		"gwas_winners_curse_results.txt");
  set_parameter("discovery-threshold",
		"dt",
		"p-value cutoff for bringing variants forward to replication",
		"5e-8");
  set_parameter("stderr-precision-factor",
		"spf",
		"factor adjusting effective precision of internal calculations",
		"100.0");
  set_parameter("trait-mean",
		"tm",
		"effective mean of trait used in regression (if null, per-line in input)",
		"");
  set_parameter("threshold-precision-factor",
		"tpf",
		"how many times larger than the discovery threshold a discovery p can be",
		"250");
  
#ifdef HAVE_PTHREAD
  set_parameter("max-threads", 
		"mt", 
		"Number of threads available (when compiled with pthread support)", 
		"2");
#endif //HAVE_PTHREAD
  //open the log file
  //Logfile::openLog(par::_out);
  if (warn_unused)
    print_ignored();
}
void gwas_winners_curse::cargs::set_parameter(const std::string &parameter_long, 
			      const std::string &parameter_short, 
			      const std::string &description, 
			      const std::string &default_value) {
  std::string value = default_value;
  if (parameter_long.empty())
    throw std::domain_error("gwas_winners_curse::cargs::set_parameter: empty parameter name");
  if (find("--" + parameter_long)) {
    value = get_formatted_value<std::string>("--" + parameter_long);
  } else if (!parameter_short.empty() &&
	     find("-" + parameter_short)) {
    value = get_formatted_value<std::string>("-" + parameter_short);
  }
  parameters::set_parameter(parameter_long, parameter_short, description, value);
}
void gwas_winners_curse::cargs::set_flag(const std::string &parameter_long, 
			 const std::string &parameter_short, 
			 const std::string &description, 
			 bool default_value) {
  parameters::set_flag(parameter_long,
		       parameter_short,
		       description,
		       find("--" + parameter_long) || 
		       (!parameter_short.empty() && find("-" + parameter_short)) 
		       ? !default_value 
		       : default_value);
}
void gwas_winners_curse::cargs::set_defaults() throw (std::domain_error) {
  _found.resize(num_args(), false);
}
bool gwas_winners_curse::cargs::find(const std::string &key) throw (std::domain_error) {
  for (int i = 1; i < num_args(); ++i) {
    if (!std::string(arg_array()[i]).compare(key)) {
      if (_found.at(i)) {
	//throw std::domain_error("User: double use of flag \"" + key + "\"");
      }
      return (_found.at(i) = true);
    }
  }
  return false;
}
void gwas_winners_curse::cargs::print_ignored() const throw () {
  for (int i = 1; i < num_args(); ++i) {
    if (!_found.at(i)) {
      std::cout << "unused command line entry: \"" 
		<< std::string(arg_array()[i]) << "\"" << std::endl;
      //Logfile::printLog("Unused command line parameter: \""
      //		+ std::string(arg_array()[i]) + "\"");
    }
  }
}
