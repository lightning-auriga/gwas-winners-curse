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

#ifndef __gwas_winners_curse_FILEINTERFACE_FILEINTERFACE_WRITER_BZIP2_H__
#define __gwas_winners_curse_FILEINTERFACE_FILEINTERFACE_WRITER_BZIP2_H__

#include "fileinterface/config.h"
#ifdef FILEINTERFACE_HAVE_LIBBZ2

#include "fileinterface/fileinterface_writer_parent.h"

#include <string>
#include <stdexcept>
#include <cstdlib>
#include <cstdio>
#include <bzlib.h>

namespace gwas_winners_curse {
  class fileinterface_writer_bzip2 : public fileinterface_writer {
  public:
    fileinterface_writer_bzip2()
      : fileinterface_writer(),
      _raw_output(0), _bz_output(0) {}
    ~fileinterface_writer_bzip2() throw() {close();}

    void open(const char *filename);
    void close();
    void clear();
    bool is_open() const;
    void put(char c);
    void writeline(const std::string &);
    void write(char *, std::streamsize);
    bool eof() const {return false;}
    bool good() const {return _good;}
    bool fail() const {return _fail;}
    bool bad() const {return _bad;}
  private:
    FILE *_raw_output;
    BZFILE *_bz_output;
  };
}


#endif //FILEINTERFACE_HAVE_LIBBZ2

#endif //__gwas_winners_curse_FILEINTERFACE_FILEINTERFACE_WRITER_BZIP2_H__
