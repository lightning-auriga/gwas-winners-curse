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

#ifndef __GWAS_WINNERS_CURSE_PARAMETERS_H__
#define __GWAS_WINNERS_CURSE_PARAMETERS_H__
#include <iostream>
#include <vector>
#include <map>
#include <utility>
#include <string>
#include <stdexcept>
#include "gwas_winners_curse/config.h"
#include "gwas_winners_curse/utilities.h"
#ifdef HAVE_PTHREAD
#include <pthread.h>
#include "gwas_winners_curse/pthread_mutex_wrapper.h"
#endif //HAVE_PTHREAD

namespace gwas_winners_curse {
  template <class value_type>
    class parameter_annotation {
  public:
    parameter_annotation() throw() :
    _full_flag(""),
      _truncated_flag(""),
      _description("") {}
    parameter_annotation(const parameter_annotation &obj) throw(std::bad_alloc) :
    _full_flag(obj._full_flag),
      _truncated_flag(obj._truncated_flag),
      _description(obj._description),
      _annotation(obj._annotation) {}
    ~parameter_annotation() throw() {}
    
    std::string get_longform() const throw(std::bad_alloc) {return _full_flag;}
    void set_longform(const std::string &str) throw(std::bad_alloc) {_full_flag = str;}
    std::string get_shortform() const throw(std::bad_alloc) {return _truncated_flag;}
    void set_shortform(const std::string &str) throw(std::bad_alloc) {_truncated_flag = str;}
    std::string get_description() const throw(std::bad_alloc) {return _description;}
    void set_description(const std::string &str) throw(std::bad_alloc) {_description = str;}
    value_type get_annotation() const {return _annotation;}
    void set_annotation(const value_type &obj) {_annotation = obj;}
  private:
    std::string _full_flag;
    std::string _truncated_flag;
    std::string _description;
    value_type _annotation;
  };
  
  
  class parameters {
  public:
    parameters() throw() {}
    ~parameters() throw() {}
    
    static bool has_parameter(const std::string &parameter_nam) throw();
    static std::string get_parameter(const std::string &parameter_name) 
      throw(std::domain_error);
    static void set_parameter(const std::string &parameter_long,
			      const std::string &parameter_short,
			      const std::string &parameter_description,
			      const std::string &parameter_value)
      throw(std::domain_error, std::bad_alloc);
    static bool has_flag(const std::string &flag_name) throw();
    static bool get_flag(const std::string &flag_name) throw(std::domain_error);
    static void set_flag(const std::string &flag_long,
			 const std::string &flag_short,
			 const std::string &flag_description,
			 bool               flag_value) 
      throw(std::domain_error, std::bad_alloc);
    static void print_license(std::ostream &out);
    static void print_help(std::ostream &out);
#ifdef HAVE_PTHREAD
    static void acquire_parameter_lock() throw(std::domain_error) {_mutex.lock();}
    static void release_parameter_lock() throw(std::domain_error) {_mutex.unlock();}
#endif //HAVE_PTHREAD
  private:
#ifdef HAVE_PTHREAD
    static pthread_mutex _mutex;
#endif //HAVE_PTHREAD
    //static std::vector<std::string> _parameters;
    static std::map<std::string, parameter_annotation<std::string> > _parameters;
    static std::map<std::string, parameter_annotation<bool> > _flags;
  };
}
#endif //__GWAS_WINNERS_CURSE_PARAMETERS_H__
