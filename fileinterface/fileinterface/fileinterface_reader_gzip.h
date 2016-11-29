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

#ifndef __gwas_winners_curse_FILEINTERFACE_FILEINTERFACE_READER_GZIP_H__
#define __gwas_winners_curse_FILEINTERFACE_FILEINTERFACE_READER_GZIP_H__

#include "fileinterface/config.h"
#ifdef FILEINTERFACE_HAVE_LIBZ
#include "fileinterface/fileinterface_reader_parent.h"

#include <string>
#include <stdexcept>
#include <cstdio>
#include <zlib.h>

#include "fileinterface/helper.h"

namespace gwas_winners_curse {
  class fileinterface_reader_gzip : public fileinterface_reader {
  public:
    fileinterface_reader_gzip();
    ~fileinterface_reader_gzip() throw() {close(); delete [] _buf;}
    
    void open(const char *filename);
    void close();
    void clear();
    bool is_open() const;
    char get();
    bool getline(std::string &);
    bool eof() const {return _eof;}
    bool good() const {return _good;}
    bool bad() const {return _bad;}
    void read(char *, std::streamsize);
  private:
    gzFile _gz_input;
    bool _eof;
    char *_buf;
    unsigned _buf_max;
  };
}


#endif //HAVE_LIBZ

#endif //__gwas_winners_curse_FILEINTERFACE_FILEINTERFACE_READER_GZIP_H__
