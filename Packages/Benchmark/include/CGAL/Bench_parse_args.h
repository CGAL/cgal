// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>
#ifndef BENCH_PARSE_ARGS_H
#define BENCH_PARSE_ARGS_H

#include <stdlib.h>
#include <string>
#include <iostream>
#include "getopt.h"
#include "CGAL/Dir_search.h"

CGAL_BEGIN_NAMESPACE

/*!
 */
class Bench_parse_args {
public:
  /*!
   */
  enum IOId  {
    IO_FORMAT = 0,
    IO_F
  };

  /*!
   */
  enum FormatId  {
    FORMAT_RAT = 0,
    FORMAT_INT,
    FORMAT_FLT,
    FORMAT_R,
    FORMAT_I,
    FORMAT_F
  };

  /*!
   */
  enum BenchId {
    BENCH_TYPE_NAME = 0,
    BENCH_TYPE_MASK,
    BENCH_STRATEGY_NAME,
    BENCH_STRATEGY_MASK,
    BENCH_HEADER,
    BENCH_NAME_LENGTH,
    BENCH_TN,
    BENCH_TM,
    BENCH_SN,
    BENCH_SM,
    BENCH_H,
    BENCH_NL
  };

  /*!
   */
  enum TypeId {
    TYPE_INCREMENT = 0,
    TYPE_AGGREGATE,
    TYPE_DISPLAY,
    TYPE_POINT_LOCATION,
    TYPE_SUBCURVES,
    TYPE_POINTS,
    TYPE_I,
    TYPE_A,
    TYPE_D,
    TYPE_L,
    TYPE_C,
    TYPE_P
  };

  /*!
   */
  enum StrategyId {
    STRATEGY_TRAPEZOIDAL = 0,
    STRATEGY_NAIVE,
    STRATEGY_WALK,
    STRATEGY_SIMPLE,
    STRATEGY_TRIANGLE,
    STRATEGY_DUMMY,
    STRATEGY_T,
    STRATEGY_N,
    STRATEGY_W,
    STRATEGY_S,
    STRATEGY_G,
    STRATEGY_D
  };
    
  /*!
   */
  enum BoolId {
    BOOL_TRUE = 0,
    BOOL_FALSE,
    BOOL_T,
    BOOL_F
  };

  enum MaxFilesNumber {
    MAX_FILES_NUMBER = 20
  };

public:    
  Bench_parse_args(int argc, char * argv[]) :
      m_argc(argc), m_argv(argv), 
      m_verbose(false),
      m_postscript(false),
      m_type_mask(0xffffffff), m_strategy_mask(0xffffffff),
      m_input_format(FORMAT_RAT),
      m_print_header(true), m_files_number(0),
      m_name_length(32), m_seconds(1), 
      m_samples(0),
#if (defined _MSC_VER)
      m_iterations(10)
#else
      m_iterations(0)
#endif
  {
    m_prog_name = strrchr(argv[0], '\\');
    m_prog_name = (m_prog_name) ? m_prog_name+1 : argv[0];

    m_dirs.add(".");
    const char * root = getenv("ROOT");
    if (root) m_dirs.add(std::string(root) + "/data/Segments_2");
    if (root) m_dirs.add(std::string(root) + "/data/Conics_2");
    if (root) m_dirs.add(std::string(root) + "/data/Polylines_2");

    for (int j=0; j<MAX_FILES_NUMBER; j++) {
      m_filename[j] = 0;
    }
  }

  /*!
   */
  unsigned int get_type_mask() const { return m_type_mask; }
  unsigned int get_strategy_mask() const { return m_strategy_mask; }
  bool get_verbose() const { return m_verbose; }
  bool get_postscript() const { return m_postscript; }
  FormatId get_input_format() const { return m_input_format; }
  int get_samples() const { return m_samples; }
  int get_iterations() const { return m_iterations; }
  int get_seconds() const { return m_seconds; }
  int get_files_number() const { return m_files_number; }
  bool get_print_header()const { return m_print_header; }
  int get_name_length()const { return m_name_length; }
  const char * get_type_name(TypeId id) const { return s_type_opts[id]; }
  const char * get_strategy_name(StrategyId id) const
  { return s_strategy_opts[id]; }

  const char * get_full_name(int file_index = 0) const 
  { 
    if (file_index >= MAX_FILES_NUMBER) {
      std::cerr << "too large file index" << std::endl;
      return NULL;
    }
    return m_fullname[file_index].c_str(); 
  }
  /*!
   */
  const char * get_file_name(int file_index = 0) const
  {
    if ((file_index >= MAX_FILES_NUMBER) ||
	(!m_filename[file_index])) {
      std::cerr << "Data file missing!" << std::endl;
      return NULL;
    }
    return m_filename[file_index];
  }
  
  /*!
   */
  void printHelp(void)
  {
    printf("Usage: %s [options] data-file\n\
  -b <options>\tset bench options\n\
  \t\ttype_name=<type>\n\
  \t\ttn=<type>\tset bench type to <type> (default all)\n\
  \t\t\t\t<type> is one of:\n\
  \t\t\t\t\ti[ncrement]\t(0x1)\n\
  \t\t\t\t\ta[ggregate]\t(0x2)\n\
  \t\t\t\t\td[isplay]\t(0x4)\n\
  \t\t\t\t\t[point_]l[ocation]\t(0x8)\n\
  \t\ttype_mask=<mask>\n\
  \t\ttm=<mask>\tset bench type mask to <mask>\n\
  \t\tstrategy_name=<strategy>\n\
  \t\tsn=<strategy>\tset bench strategy to <strategy> (default all)\n\
  \t\t\t\t<strategy> is one of:\n\
  \t\t\t\t\tt[rapezoidal]\t(0x1)\n\
  \t\t\t\t\tn[aive]\t\t(0x2)\n\
  \t\t\t\t\tw[alk]\t\t(0x4)\n\
  \t\t\t\t\ts[imple]\t(0x8)\n\
  \t\t\t\t\t[trian]g[le]\t(0x10)\n\
  \t\t\t\t\td[ummy]\t\t(0x20)\n\
  \t\tstrategy_mask=<mask>\n\
  \t\tsm=<mask>\tset bench strategy mask to <mask>\n\
  \t\th[eader]=<bool>\tprint header (default true)\n\
  \t\tname_length=<length>\n\
  \t\tnl=<length>\tset the length of the name field to <length\n\
  -d <dir>\tadd directory <dir> to list of search directories\n\
  -h\t\tprint this help message\n\
  -i <iters>\tset number of iterations to <iters> (default 0)\n\
  -I <options>\tset input options\n\
  \t\tf[ormat]=<format>\tset format to <format> (default rat)\n\
  \t\t\t\t<format is one of {i[nt], f[lt], r[at]}\n\
  -s <samples>\tset number of samples to <samples> (default 10)\n\
  -t <seconds>\tset number of seconds to <seconds> (default 1)\n\
  -v\t\ttoggle verbosity (default false)\n",
           m_prog_name);
  }

  int parse()
  {
    int c;
    while ((c = getopt(m_argc, m_argv, s_option_str)) != EOF) {
      switch (c) {
        case 'b': if (get_bench_param(optarg) < 0) return -1; break;
        case 'd': m_dirs.add(optarg); break;
        case 'h': printHelp(); return 1;
        case 'I': if (get_io_param(optarg)) return -1; break;
        case 'i': m_iterations = atoi(optarg); break;
        case 'p': m_postscript = !m_postscript; break;
        case 's': m_samples = atoi(optarg); break;
        case 't': m_seconds = atoi(optarg); break;
        case 'v': m_verbose = !m_verbose; break;         
	  //case 'f': m_files_number = atoi(optarg); break;
        default:
          std::cerr << m_prog_name << ": invalid option -- "
                    << static_cast<char>(c) << std::endl;
          std::cerr << "Try `" << m_prog_name << " -h' for more information."
                    << std::endl;
          return -1;
      }
    }

    while (optind < m_argc) {
      m_filename[m_files_number] = m_argv[optind];
      if (!m_dirs.find(m_filename[m_files_number], m_fullname[m_files_number])) {
        std::cerr << "Cannot find file " << m_filename << "!" << std::endl;
        return -1;
      }
      m_files_number++;
      optind++;
      if (m_files_number > MAX_FILES_NUMBER) {
	std::cerr << "Too many input files !" << std::endl;
        return -1;
      }
    }
     
    return 0;
  }
    
private:
  /*!
   */
  int get_io_param(char * optarg)
  {
    char * options = optarg;
    char * value = 0;
    if (*options == '\0') return 0;
    while (*options != '\0') {
      switch (getsubopt(&options, s_io_opts, &value)) {
        case IO_FORMAT:
        case IO_F: if (get_format_param(value) < 0) return -1; break;
        default:
          std::cerr << "Unrecognized IO option '" << optarg << "'!"
                    << std::endl;
          std::cerr << "Try `" << m_prog_name << " -h' for more information."
                    << std::endl;
          return -1;
      }
    }
    return 0;
  }

  /*!
   */
  int get_format_param(char * value)
  {
    if (value) {
      char * subvalue;
      switch (getsubopt(&value, s_format_opts, &subvalue)) {
        case FORMAT_INT:
        case FORMAT_I: m_input_format = FORMAT_INT; return 0;
        case FORMAT_RAT:
        case FORMAT_R: m_input_format = FORMAT_RAT; return 0;
        case FORMAT_FLT:
        case FORMAT_F: m_input_format = FORMAT_FLT; return 0;
      }
    }
    std::cerr << "Unrecognized Format option '" << optarg << "'!" << std::endl;
    std::cerr << "Try `" << m_prog_name << " -h' for more information."
              << std::endl;
    return -1;
  }

  /*!
   */
  int get_bool_param(bool & param, char * value)
  {
    if (value) {
      char * subvalue;
      switch (getsubopt(&value, s_bool_opts, &subvalue)) {
        case BOOL_TRUE:
        case BOOL_T: param = true; return 0;
        case BOOL_FALSE:
        case BOOL_F: param = false; return 0;
      }
    }
    std::cerr << "Unrecognized Bool option '" << optarg << "'!" << std::endl;
    std::cerr << "Try `" << m_prog_name << " -h' for more information."
              << std::endl;
    return -1;
  }    
    
  /*!
   */
  int get_bench_param(char * optarg)
  {
    char * options = optarg;
    char * value = 0;
    if (*options == '\0') return 0;
    while (*options != '\0') {
      switch(getsubopt(&options, s_bench_opts, &value)) {
        case BENCH_TYPE_NAME:
        case BENCH_TN: if (getType_name_param(value) < 0) return -1; break;

        case BENCH_STRATEGY_NAME:
        case BENCH_SN:
         if (get_strategy_name_param(value) < 0) return -1; break;
          
        case BENCH_TYPE_MASK:
        case BENCH_TM:
          if (!value) goto err;
          m_type_mask = strtoul(value, 0, 0);
          break;

        case BENCH_STRATEGY_MASK:
        case BENCH_SM:
          if (!value) goto err;
          m_strategy_mask = strtoul(value, 0, 0);
          break;

        case BENCH_HEADER:
        case BENCH_H:
          if (get_bool_param(m_print_header, value) < 0) return -1;
          break;

        case BENCH_NAME_LENGTH:
        case BENCH_NL:
          if (!value) goto err;
          m_name_length = atoi(value);
          break;
        default:
      err:
          std::cerr << "Unrecognized Bench option '" << optarg << "'!"
                    << std::endl;
          std::cerr << "Try `" << m_prog_name << " -h' for more information."
                    << std::endl;
          return -1;
      }
    }
    return 0;
  }

  /*
   */
  int getType_name_param(char * value)
  {
    if (value) {
      char * subvalue;
      switch (getsubopt(&value, s_type_opts, &subvalue)) {
        case TYPE_INCREMENT:
        case TYPE_I:
          m_type_mask = 0x1 << TYPE_INCREMENT; return 0;
        case TYPE_AGGREGATE:
        case TYPE_A:
          m_type_mask = 0x1 << TYPE_AGGREGATE; return 0;
        case TYPE_DISPLAY:
        case TYPE_D:
          m_type_mask = 0x1 << TYPE_DISPLAY; return 0;
        case TYPE_POINT_LOCATION:
        case TYPE_L:
	  m_type_mask = 0x1 << TYPE_POINT_LOCATION; return 0;
        case TYPE_SUBCURVES:
        case TYPE_C:
	  m_type_mask = 0x1 << TYPE_SUBCURVES; return 0;
        case TYPE_POINTS:
        case TYPE_P:
	  m_type_mask = 0x1 << TYPE_POINTS; return 0;
      }
    }
    std::cerr << "Unrecognized Bench name option '" << optarg << "'!"
              << std::endl;
    std::cerr << "Try `" << m_prog_name << " -h' for more information."
              << std::endl;
    return -1;
  }

  /*
   */
  int get_strategy_name_param(char * value)
  {
    if (value) {
      char * subvalue;
      switch (getsubopt(&value, s_strategy_opts, &subvalue)) {
        case STRATEGY_TRAPEZOIDAL:
        case STRATEGY_T:
          m_strategy_mask = 0x1 << STRATEGY_TRAPEZOIDAL; return 0;
        case STRATEGY_NAIVE:
        case STRATEGY_N:
          m_strategy_mask = 0x1 << STRATEGY_NAIVE; return 0;
        case STRATEGY_WALK:
        case STRATEGY_W:
          m_strategy_mask = 0x1 << STRATEGY_WALK; return 0;
        case STRATEGY_SIMPLE:
        case STRATEGY_S:
          m_strategy_mask = 0x1 << STRATEGY_SIMPLE; return 0;
        case STRATEGY_TRIANGLE:
        case STRATEGY_G:
          m_strategy_mask = 0x1 << STRATEGY_TRIANGLE; return 0;
        case STRATEGY_DUMMY:
        case STRATEGY_D:
          m_strategy_mask = 0x1 << STRATEGY_DUMMY; return 0;
      }
    }
    std::cerr << "Unrecognized Bench Name option '" << optarg << "'!"
              << std::endl;
    std::cerr << "Try `" << m_prog_name << " -h' for more information."
              << std::endl;
    return -1;
  }

private:
  static char * s_format_opts[];
  static char * s_io_opts[];
  static char * s_bench_opts[];
  static char * s_type_opts[];
  static char * s_strategy_opts[];
  static char * s_bool_opts[];
    
private:
  static char s_option_str[];
  const char * m_prog_name;

  std::string m_fullname[MAX_FILES_NUMBER];
  Dir_search m_dirs;

  int m_argc;
  char ** m_argv;
  const char * m_filename[MAX_FILES_NUMBER];
  bool m_verbose;
  bool m_postscript;
  unsigned int m_type_mask;
  unsigned int m_strategy_mask;
  FormatId m_input_format;
  bool m_print_header;
  int m_files_number;
  int m_name_length;
  int m_seconds;
  int m_samples;
  int m_iterations;
};

CGAL_END_NAMESPACE

#endif
