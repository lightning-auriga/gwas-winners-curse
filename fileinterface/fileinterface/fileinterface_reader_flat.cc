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

#include "fileinterface/fileinterface_reader_flat.h"

void gwas_winners_curse::fileinterface_reader_flat::open(const char *filename) {
  _input.open(filename);
  if (!_input.is_open())
    throw std::domain_error("gwas_winners_curse::fileinterface_reader_flat::open: cannot open file \""
			    + std::string(filename) + "\"");
}

bool gwas_winners_curse::fileinterface_reader_flat::getline(std::string &line) {
  line = "";
  if (_input.peek() == EOF)
    return false;
  std::getline(_input, line);
  return true;
}

void gwas_winners_curse::fileinterface_reader_flat::read(char *buf, std::streamsize n) {
  _input.read(buf, n);
}
