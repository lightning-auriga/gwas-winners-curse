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

#ifndef __gwas_winners_curse_FILEINTERFACE_PLINKBED_H__
#define __gwas_winners_curse_FILEINTERFACE_PLINKBED_H__

#include <string>
#include <stdexcept>
#include <cstdlib>
#include <iostream>

#include "fileinterface/fileinterface_reader.h"
#include "fileinterface/fileinterface_writer.h"
#include "fileinterface/helper.h"

namespace gwas_winners_curse {
  typedef enum {
    NONE, GZIP, BZIP2
  } COMPRESSION_TYPE;

  inline COMPRESSION_TYPE interpret_compression(const std::string &s) {
    if (!s.compare("gzip"))
      return GZIP;
    if (!s.compare("bzip2"))
      return BZIP2;
    return NONE;
  }

  inline char plink_packed_to_char(unsigned val) {
    if (val == 0) return '2';
    if (val == 1) return '9';
    if (val == 2) return '1';
    return '0';
  }

  inline unsigned char char_to_plink_packed(char c) {
    if (c == '2') return (unsigned char)0;
    if (c == '1') return (unsigned char)2;
    if (c == '0') return (unsigned char)3;
    return (unsigned char)1;
  }

  class plinkbed_reader : public gwas_winners_curse::fileinterface_reader {
  public:
    plinkbed_reader(unsigned, const std::string &cl = "");
    ~plinkbed_reader() throw() {close();}

    void open(const char *filename);
    void close();
    void clear();
    bool is_open() const;
    char get();
    bool getline(std::string &);
    bool eof() const;
    bool good() const;
    bool bad() const;
    void read(char *, std::streamsize);

  private:
    gwas_winners_curse::fileinterface_reader *_input;
    gwas_winners_curse::COMPRESSION_TYPE _compress;
    char *_buf;
    unsigned _line_length;
    unsigned _line_length_packed;
  };

  class plinkbed_writer : public gwas_winners_curse::fileinterface_writer {
  public:
    plinkbed_writer(const std::string &cl = "");
    ~plinkbed_writer() throw() {close();}

    void open(const char *filename);
    void close();
    void clear();
    bool is_open() const;
    void put(char c);
    void writeline(const std::string &);
    void write(char *, std::streamsize);
    bool eof() const;
    bool good() const;
    bool fail() const;
    bool bad() const;

  private:
    gwas_winners_curse::fileinterface_writer *_output;
    gwas_winners_curse::COMPRESSION_TYPE _compress;
    char *_buf;
    unsigned _line_length;
    unsigned _line_length_packed;
  };
}

#endif //__gwas_winners_curse_FILEINTERFACE_PLINKBED_H__
