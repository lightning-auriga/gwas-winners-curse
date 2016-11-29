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

#ifndef __FILEINTERFACE_HELPER_H__
#define __FILEINTERFACE_HELPER_H__

#include <string>
#include <sstream>
#include <stdexcept>

inline std::string get_newline() {
#ifdef _WIN64
  return "\r\n";
#elif  _WIN32
  return "\r\n";
#else //if __APPLE__  // linux
  //#elif __linux || __unix
  return "\n";
#endif
}

template <class value_type>
std::string fi_to_string(const value_type &obj) {
  std::ostringstream o;
  if (!(o << obj))
    throw std::domain_error("to_string: cannot convert object to string");
  return o.str();
}

#endif //__FILEINTERFACE_HELPER_H__
