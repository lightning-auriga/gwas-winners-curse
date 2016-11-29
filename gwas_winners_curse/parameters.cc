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

#include "gwas_winners_curse/parameters.h"

bool gwas_winners_curse::parameters::has_parameter(const std::string &parameter_name) throw() {
#ifdef HAVE_PTHREAD
  acquire_parameter_lock();
#endif //HAVE_PTHREAD
  bool return_value = _parameters.find(parameter_name) != _parameters.end();
#ifdef HAVE_PTHREAD
  release_parameter_lock();
#endif //HAVE_PTHREAD
  return return_value;
}
std::string gwas_winners_curse::parameters::get_parameter(const std::string &parameter_name) 
  throw(std::domain_error) {
#ifdef HAVE_PTHREAD
  //DO NOT CALL has_parameter WHEN LOCKED!!!
  acquire_parameter_lock();  
#endif //HAVE_PTHREAD
  std::string return_value = "";
  std::map<std::string, parameter_annotation<std::string> >::const_iterator finder;
  if ((finder = _parameters.find(parameter_name)) == _parameters.end())
    throw std::domain_error("get_parameter: no parameter named \"" + parameter_name + "\"");
  return_value = finder->second.get_annotation();
#ifdef HAVE_PTHREAD
  release_parameter_lock();
#endif //HAVE_PTHREAD
  return return_value;
}
void gwas_winners_curse::parameters::set_parameter(const std::string &parameter_long,
			       const std::string &parameter_short,
			       const std::string &parameter_description,
			       const std::string &parameter_value)
  throw(std::domain_error, std::bad_alloc) {
#ifdef HAVE_PTHREAD
  //DO NOT CALL has_parameter WHEN LOCKED!!!
  acquire_parameter_lock();
#endif //HAVE_PTHREAD
  parameter_annotation<std::string> annotation;
  annotation.set_longform(parameter_long);
  annotation.set_shortform(parameter_short);
  annotation.set_description(parameter_description);
  annotation.set_annotation(parameter_value);
  std::map<std::string, parameter_annotation<std::string> >::iterator finder;
  if ((finder = _parameters.find(parameter_long)) == _parameters.end()) {
    finder = _parameters.insert(std::make_pair(parameter_long, annotation)).first;
  } else finder->second = annotation;
#ifdef HAVE_PTHREAD
  release_parameter_lock();
#endif //HAVE_PTHREAD
}
bool gwas_winners_curse::parameters::has_flag(const std::string &flag_name) throw() {
#ifdef HAVE_PTHREAD
  acquire_parameter_lock();
#endif //HAVE_PTHREAD
  bool return_value = _flags.find(flag_name) != _flags.end();
#ifdef HAVE_PTHREAD
  release_parameter_lock();
#endif //HAVE_PTHREAD
  return return_value;
}
bool gwas_winners_curse::parameters::get_flag(const std::string &flag_name) throw(std::domain_error) {
#ifdef HAVE_PTHREAD
  //DO NOT CALL has_flag WHEN LOCKED!!!
  acquire_parameter_lock();
#endif //HAVE_PTHREAD
  bool return_value = false;
  std::map<std::string, parameter_annotation<bool> >::const_iterator finder;
  if ((finder = _flags.find(flag_name)) == _flags.end())
    throw std::domain_error("get_flag: no flag named \"" + flag_name + "\"");
  return_value = finder->second.get_annotation();
#ifdef HAVE_PTHREAD
  release_parameter_lock();
#endif //HAVE_PTHREAD
  return return_value;
}
void gwas_winners_curse::parameters::set_flag(const std::string &flag_long,
			      const std::string &flag_short,
			      const std::string &flag_description,
			      bool               flag_value) 
  throw(std::domain_error, std::bad_alloc) {
#ifdef HAVE_PTHREAD
  //DO NOT CALL has_flag WHEN LOCKED!!!
  acquire_parameter_lock();
#endif //HAVE_PTHREAD
  parameter_annotation<bool> annotation;
  annotation.set_longform(flag_long);
  annotation.set_shortform(flag_short);
  annotation.set_description(flag_description);
  annotation.set_annotation(flag_value);
  std::map<std::string, parameter_annotation<bool> >::iterator finder;
  if ((finder = _flags.find(flag_long)) == _flags.end())
    _flags.insert(std::make_pair(flag_long, annotation));
  else finder->second = annotation;
#ifdef HAVE_PTHREAD
  release_parameter_lock();
#endif //HAVE_PTHREAD
}
void gwas_winners_curse::parameters::print_license(std::ostream &out) {
  out << GWAS_WINNERS_CURSE_PACKAGE_STRING;
  out << "\nCopyright © 2015 Cameron Palmer \n\
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n\
This is free software: you are free to change and redistribute it.\n\
There is NO WARRANTY, to the extent permitted by law." << std::endl;
}
void gwas_winners_curse::parameters::print_help(std::ostream &out) {
  std::string longflag_header = "Full Flag", shortflag_header = "Abbreviation", 
    description_header = "Description", value_header = "Current Value";
  unsigned longflag_max_width = longflag_header.size(), 
    shortflag_max_width = shortflag_header.size(),
    description_max_width = description_header.size(), value_max_width = value_header.size();
  unsigned longflag_tab = 0, shortflag_tab = 0, description_tab = 0;
  std::string no_description_message = "no description available";
  unsigned space_per_tab = 8;
  
  for (std::map<std::string, parameter_annotation<std::string> >::const_iterator 
	 iter = _parameters.begin();
       iter != _parameters.end(); ++iter) {
    longflag_max_width = max(longflag_max_width, iter->second.get_longform().size() + 2);
    shortflag_max_width = max(shortflag_max_width, iter->second.get_shortform().empty() 
			      ? 4 
			      : iter->second.get_shortform().size() + 3);
    description_max_width = max(description_max_width, iter->second.get_description().empty() 
				? no_description_message.size() 
				: iter->second.get_description().size());
    value_max_width = max(value_max_width, iter->second.get_annotation().empty() 
			  ? 4 
			  : iter->second.get_annotation().size());
  }
  for (std::map<std::string, parameter_annotation<bool> >::const_iterator 
	 iter = _flags.begin();
       iter != _flags.end(); ++iter) {
    longflag_max_width = max(longflag_max_width, iter->second.get_longform().size() + 2);
    shortflag_max_width = max(shortflag_max_width, iter->second.get_shortform().empty() 
			      ? 4 
			      : iter->second.get_shortform().size() + 3);
    description_max_width = max(description_max_width, iter->second.get_description().empty() 
				? no_description_message.size() 
				: iter->second.get_description().size());
    value_max_width = max(value_max_width, iter->second.get_annotation() ? 4 : 5);
  }

  ++longflag_max_width;
  longflag_max_width = longflag_max_width / space_per_tab + 
    (longflag_max_width % space_per_tab ? 1 : 0);
  ++shortflag_max_width;
  shortflag_max_width = shortflag_max_width / space_per_tab + 
    (shortflag_max_width % space_per_tab ? 1 : 0);
  ++description_max_width;
  description_max_width = description_max_width / space_per_tab + 
    (description_max_width % space_per_tab ? 1 : 0);
  //++value_max_width;



  out << longflag_header;
  for (unsigned i = 0;
       i < convert_to_tabs(longflag_header.size(), longflag_max_width, space_per_tab); 
       ++i) out << '\t';
  out << shortflag_header;
  for (unsigned i = 0; 
       i < convert_to_tabs(shortflag_header.size(), shortflag_max_width, space_per_tab);
       ++i) out << '\t';
  out << description_header;
  for (unsigned i = 0;
       i < convert_to_tabs(description_header.size(), description_max_width, space_per_tab);
       ++i) out << '\t';
  out << value_header << std::endl;

  for (std::map<std::string, parameter_annotation<std::string> >::const_iterator 
	 iter = _parameters.begin();
       iter != _parameters.end(); ++iter) {
    longflag_tab = iter->second.get_longform().size() + 2;
    longflag_tab = convert_to_tabs(longflag_tab, longflag_max_width, space_per_tab);
    shortflag_tab = iter->second.get_shortform().empty() 
      ? 4 
      : iter->second.get_shortform().size() + 3;
    shortflag_tab = convert_to_tabs(shortflag_tab, shortflag_max_width, space_per_tab);
    description_tab = iter->second.get_description().empty() 
      ? no_description_message.size() 
      : iter->second.get_description().size();
    description_tab = convert_to_tabs(description_tab, description_max_width, space_per_tab);

    out << "--" << iter->second.get_longform();
    for (unsigned i = 0; i < longflag_tab; ++i) out << '\t';
    out << '[' << (iter->second.get_shortform().empty() 
		   ? std::string("NA") 
		   : "-" + iter->second.get_shortform()) << "]";
    for (unsigned i = 0; i < shortflag_tab; ++i) out << '\t';
    out << (iter->second.get_description().empty() 
	    ? no_description_message 
	    : iter->second.get_description());
    for (unsigned i = 0; i < description_tab; ++i) out << '\t';
    out << (iter->second.get_annotation().empty() 
	    ? std::string("null") 
	    : iter->second.get_annotation()) << std::endl;
  }
  for (std::map<std::string, parameter_annotation<bool> >::const_iterator 
	 iter = _flags.begin();
       iter != _flags.end(); ++iter) {
    longflag_tab = iter->second.get_longform().size() + 2;
    longflag_tab = convert_to_tabs(longflag_tab, longflag_max_width, space_per_tab);
    shortflag_tab = iter->second.get_shortform().empty() 
      ? 4 
      : iter->second.get_shortform().size() + 3;
    shortflag_tab = convert_to_tabs(shortflag_tab, shortflag_max_width, space_per_tab);
    description_tab = iter->second.get_description().empty() 
      ? no_description_message.size() 
      : iter->second.get_description().size();
    description_tab = convert_to_tabs(description_tab, description_max_width, space_per_tab);

    out << "--" << iter->second.get_longform();
    for (unsigned i = 0; i < longflag_tab; ++i) out << '\t';
    out << '[' << (iter->second.get_shortform().empty() 
		   ? std::string("NA") 
		   : "-" + iter->second.get_shortform()) << "]";
    for (unsigned i = 0; i < shortflag_tab; ++i) out << '\t';
    out << (iter->second.get_description().empty() 
	    ? no_description_message 
	    : iter->second.get_description());
    for (unsigned i = 0; i < description_tab; ++i) out << '\t';
    out << (iter->second.get_annotation() ? "true" : "false") << std::endl;
  }
}
#ifdef HAVE_PTHREAD
pthread_mutex gwas_winners_curse::parameters::_mutex;
#endif //HAVE_PTHREAD
std::map<std::string, gwas_winners_curse::parameter_annotation<std::string> > gwas_winners_curse::parameters::_parameters;
std::map<std::string, gwas_winners_curse::parameter_annotation<bool> > gwas_winners_curse::parameters::_flags;
