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

#ifndef __gwas_winners_curse_FILEINTERFACE_FILEINTERFACE_H__
#define __gwas_winners_curse_FILEINTERFACE_FILEINTERFACE_H__

#include <string>
#include <iostream>

#include "fileinterface/fileinterface_reader.h"
#include "fileinterface/fileinterface_writer.h"
#include "fileinterface/plinkbed.h"

namespace gwas_winners_curse {
  inline fileinterface_reader *reconcile_reader(const std::string &filename) {
    gwas_winners_curse::fileinterface_reader *ptr = 0;
    /*if (filename.find(".bed.gz") == filename.size() - 7) {
      ptr = new gwas_winners_curse::plinkbed_reader("gzip");
      ptr->open(filename);
    } else if (filename.find(".bed.bz2") == filename.size() - 8) {
      ptr = new gwas_winners_curse::plinkbed_reader("bzip2");
      ptr->open(filename);
      } else if (filename.find(".bed") == filename.size() - 4) {
      ptr = new gwas_winners_curse::plinkbed_reader("none");
      ptr->open(filename);
      } else */if (filename.find(".gz") == filename.size() - 3) {
#ifdef FILEINTERFACE_HAVE_LIBZ
      ptr = new gwas_winners_curse::fileinterface_reader_gzip;
      ptr->open(filename);
#else
      throw std::domain_error("zlib support not compiled into software, cannot open \"" + filename + "\"");
#endif //FILEINTERFACE_HAVE_LIBZ
    } else if (filename.find(".bz2") == filename.size() - 4) {
#ifdef FILEINTERFACE_HAVE_LIBBZ2
      ptr = new gwas_winners_curse::fileinterface_reader_bzip2;
      ptr->open(filename);
#else
      throw std::domain_error("libbz2 support not compiled into software, cannot open \"" + filename + "\"");
#endif //FILEINTERFACE_HAVE_LIBBZ2
    } else {
      ptr = new gwas_winners_curse::fileinterface_reader_flat;
      ptr->open(filename);
    }
    return ptr;
  }
  
  inline fileinterface_writer *reconcile_writer(const std::string &filename) {
    gwas_winners_curse::fileinterface_writer *ptr = 0;
    if (filename.find(".bed.gz") == filename.size() - 7) {
      ptr = new gwas_winners_curse::plinkbed_writer("gzip");
      ptr->open(filename);
    } else if (filename.find("bed.bz2") == filename.size() - 8) {
      ptr = new gwas_winners_curse::plinkbed_writer("bzip2");
      ptr->open(filename);
    } else if (filename.find(".bed") == filename.size() - 4) {
      ptr = new gwas_winners_curse::plinkbed_writer("none");
      ptr->open(filename);
    } else if (filename.find(".gz") == filename.size() - 3) {
#ifdef FILEINTERFACE_HAVE_LIBZ
      ptr = new gwas_winners_curse::fileinterface_writer_gzip;
      ptr->open(filename);
#else
      throw std::domain_error("zlib support not compiled into software, cannot write \"" + filename + "\"");
#endif //FILEINTERFACE_HAVE_LIBZ
    } else if (filename.find(".bz2") == filename.size() - 4) {
#ifdef FILEINTERFACE_HAVE_LIBBZ2
      ptr = new gwas_winners_curse::fileinterface_writer_bzip2;
      ptr->open(filename);
#else
      throw std::domain_error("libbz2 support not compiled into software, cannot write \"" + filename + "\"");
#endif //FILEINTERFACE_HAVE_LIBBZ2
    } else {
      ptr = new gwas_winners_curse::fileinterface_writer_flat;
      ptr->open(filename);
    }
    return ptr;
  }
}

#endif //__gwas_winners_curse_FILEINTERFACE_FILEINTERFACE_H__
