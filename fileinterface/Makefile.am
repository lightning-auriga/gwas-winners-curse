# Automake file

ACLOCAL_AMFLAGS = -I m4
# Build rules for libraries.
lib_LTLIBRARIES = lib/libfileinterface.la

lib_libfileinterface_la_SOURCES = fileinterface/fileinterface_reader_bzip2.cc \
  fileinterface/fileinterface_reader_flat.cc \
  fileinterface/fileinterface_writer_bzip2.cc \
  fileinterface/fileinterface_writer_flat.cc \
  fileinterface/fileinterface_reader_gzip.cc \
  fileinterface/fileinterface_writer_gzip.cc \
  fileinterface/plinkbed.cc

library_includedir = $(includedir)/libfileinterface/include

library_include_HEADERS = \
  fileinterface/fileinterface_reader.h \
  fileinterface/fileinterface_reader_parent.h \
  fileinterface/fileinterface_reader_bzip2.h \
  fileinterface/fileinterface_reader_flat.h \
  fileinterface/fileinterface_reader_gzip.h \
  fileinterface/fileinterface_writer.h \
  fileinterface/fileinterface_writer_parent.h \
  fileinterface/fileinterface_writer_bzip2.h \
  fileinterface/fileinterface_writer_flat.h \
  fileinterface/fileinterface_writer_gzip.h \
  fileinterface/helper.h \
  fileinterface/plinkbed.h \
  fileinterface/fileinterface.h \
  fileinterface/config.h

lib_libfileinterface_la_libincludedir = $(libdir)/libfileinterface/include/fileinterface

#nodist_lib_libfileinterface_la_libinclude_HEADERS = fileinterface/config.h

lib_libfileinterface_la_LDFLAGS = -version-info 1:0:0

dist_doc_DATA = README
#ACLOCAL_AMFLAGS = -I m4

#lib_libfileinterface_la_libADD = lib/libfileinterface.la
