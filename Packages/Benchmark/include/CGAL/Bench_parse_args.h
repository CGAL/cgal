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
    TYPE_I,
    TYPE_A,
    TYPE_D
  };

  /*!
   */
  enum StrategyId {
    STRATEGY_TRAPEZOIDAL = 0,
    STRATEGY_NAIVE,
    STRATEGY_WALK,
    STRATEGY_T,
    STRATEGY_N,
    STRATEGY_W
  };
    
  /*!
   */
  enum BoolId {
    BOOL_TRUE = 0,
    BOOL_FALSE,
    BOOL_T,
    BOOL_F
  };

public:    
  Bench_parse_args(int argc, char * argv[]) :
      m_argc(argc), m_argv(argv), m_filename(0),
      m_verbose(false), m_typeMask(0xffffffff), m_strategyMask(0xffffffff),
      m_inputFormat(FORMAT_RAT),
      m_printHeader(true), m_nameLength(32), m_seconds(10), m_samples(0),
#if (defined _MSC_VER)
      m_iterations(10)
#else
      m_iterations(0)
#endif
  {
    m_progName = strrchr(argv[0], '\\');
    m_progName = (m_progName) ? m_progName+1 : argv[0];

    m_dirs.add(".");
    const char * root = getenv("ROOT");
    if (root) m_dirs.add(std::string(root) + "/data/Segments_2");
  }

  /*!
   */
  unsigned int getTypeMask() const { return m_typeMask; }
  unsigned int getStrategyMask() const { return m_strategyMask; }
  bool getVerbose() const { return m_verbose; }
  FormatId getInputFormat() const { return m_inputFormat; }
  int getSamples() const { return m_samples; }
  int getIterations() const { return m_iterations; }
  int getSeconds() const { return m_seconds; }
  bool getPrintHeader()const { return m_printHeader; }
  int getNameLength()const { return m_nameLength; }
  const char * getFilename() const { return m_filename; }
  const char * getFullname() const { return m_fullname.c_str(); }
  const char * getTypeName(TypeId id) const { return s_typeOpts[id]; }
  const char * getStrategyName(StrategyId id) const { return s_strategyOpts[id]; }

  /*!
   */
  void printHelp(void)
  {
    printf("Usage: %s [options]\n\
  -b <options>\tset bench options\n\
  \t\ttype_name=<type>\n\
  \t\ttn=<type>\tset bench type to <type> (default all)\n\
  \t\t\t\t<type> is one of:\n\
  \t\t\t\t\ti[ncrement]\t(0x1)\n\
  \t\t\t\t\ta[ggregate]\t(0x2)\n\
  \t\t\t\t\td[isplay]\t(0x4)\n\
  \t\ttype_mask=<mask>\n\
  \t\ttm=<mask>\tset bench type mask to <mask>\n\
  \t\tstrategy_name=<strategy>\n\
  \t\tsn=<strategy>\tset bench strategy to <strategy> (default all)\n\
  \t\t\t\t<strategy> is one of:\n\
  \t\t\t\t\tt[rapezoidal]\t(0x1)\n\
  \t\t\t\t\tn[aive]\t\t(0x2)\n\
  \t\t\t\t\tw[alk]\t\t(0x4)\n\
  \t\tstrategy_mask=<mask>\n\
  \t\tsm=<mask>\tset bench strategy mask to <mask>\n\
  \t\th[eader]=<bool>\tprint header (default true)\n\
  \t\tname_length=<length>\n\
  \t\tnl=<length>\tset the length of the name field to <length\n\
  -h\t\tprint this help message\n\
  -i <iters>\tset number of iterations to <iters> (default 0)\n\
  -I <options>\tset input options\n\
  \t\tf[ormat]=<format>\tset format to <format> (default rat)\n\
  \t\t\t\t<format is one of {i[nt], f[lt], r[at]}\n\
  -s <samples>\tset number of samples to <samples> (default 10)\n\
  -t <seconds>\tset number of seconds to <seconds> (default 1)\n\
  -v\t\ttoggle verbosity (default false)\n",
           m_progName);
  }

  int parse()
  {
    int c;
    while ((c = getopt(m_argc, m_argv, s_optionStr)) != EOF) {
      switch (c) {
        case 'b': if (getBenchParam(optarg) < 0) return -1; break;
        case 'd': m_dirs.add(optarg); break;
        case 'h': printHelp(); return 1;
        case 'I': if (getIOParm(optarg)) return -1; break;
        case 'i': m_iterations = atoi(optarg); break;
        case 's': m_samples = atoi(optarg); break;
        case 't': m_seconds = atoi(optarg); break;
        case 'v': m_verbose = !m_verbose; break;
        default:
          std::cerr << m_progName << ": invalid option -- "
                    << static_cast<char>(c) << std::endl;
          std::cerr << "Try `" << m_progName << " -h' for more information."
                    << std::endl;
          return -1;
      }
    }
  
    if (optind >= m_argc) {
      std::cerr << "Data file missing!" << std::endl;
      return -1;
    }

    m_filename = m_argv[optind];
    if (!m_dirs.find(m_filename, m_fullname)) {
      std::cerr << "Cannot find file " << m_filename << "!" << std::endl;
      return -1;
    }
    return 0;
  }
    
private:
  /*!
   */
  int getIOParm(char * optarg)
  {
    char * options = optarg;
    char * value = 0;
    if (*options == '\0') return 0;
    while (*options != '\0') {
      switch (getsubopt(&options, s_IOOpts, &value)) {
        case IO_FORMAT:
        case IO_F: if (getFormatParam(value) < 0) return -1; break;
        default:
          std::cerr << "Unrecognized IO option '" << optarg << "'!"
                    << std::endl;
          std::cerr << "Try `" << m_progName << " -h' for more information."
                    << std::endl;
          return -1;
      }
    }
    return 0;
  }

  /*!
   */
  int getFormatParam(char * value)
  {
    if (value) {
      char * subvalue;
      switch (getsubopt(&value, s_formatOpts, &subvalue)) {
        case FORMAT_INT:
        case FORMAT_I: m_inputFormat = FORMAT_INT; return 0;
        case FORMAT_RAT:
        case FORMAT_R: m_inputFormat = FORMAT_RAT; return 0;
        case FORMAT_FLT:
        case FORMAT_F: m_inputFormat = FORMAT_FLT; return 0;
      }
    }
    std::cerr << "Unrecognized Format option '" << optarg << "'!" << std::endl;
    std::cerr << "Try `" << m_progName << " -h' for more information."
              << std::endl;
    return -1;
  }

  /*!
   */
  int getBoolParam(bool & param, char * value)
  {
    if (value) {
      char * subvalue;
      switch (getsubopt(&value, s_boolOpts, &subvalue)) {
        case BOOL_TRUE:
        case BOOL_T: param = true; return 0;
        case BOOL_FALSE:
        case BOOL_F: param = false; return 0;
      }
    }
    std::cerr << "Unrecognized Bool option '" << optarg << "'!" << std::endl;
    std::cerr << "Try `" << m_progName << " -h' for more information."
              << std::endl;
    return -1;
  }    
    
  /*!
   */
  int getBenchParam(char * optarg)
  {
    char * options = optarg;
    char * value = 0;
    if (*options == '\0') return 0;
    while (*options != '\0') {
      switch(getsubopt(&options, s_benchOpts, &value)) {
        case BENCH_TYPE_NAME:
        case BENCH_TN: if (getTypeNameParam(value) < 0) return -1; break;

        case BENCH_STRATEGY_NAME:
        case BENCH_SN: if (getStrategyNameParam(value) < 0) return -1; break;
          
        case BENCH_TYPE_MASK:
        case BENCH_TM:
          if (!value) goto err;
          m_typeMask = strtoul(value, 0, 0);
          break;

        case BENCH_STRATEGY_MASK:
        case BENCH_SM:
          if (!value) goto err;
          m_strategyMask = strtoul(value, 0, 0);
          break;

        case BENCH_HEADER:
        case BENCH_H:
          if (getBoolParam(m_printHeader, value) < 0) return -1;
          break;

        case BENCH_NAME_LENGTH:
        case BENCH_NL:
          if (!value) goto err;
          m_nameLength = atoi(value);
          break;
        default:
      err:
          std::cerr << "Unrecognized Bench option '" << optarg << "'!"
                    << std::endl;
          std::cerr << "Try `" << m_progName << " -h' for more information."
                    << std::endl;
          return -1;
      }
    }
    return 0;
  }

  /*
   */
  int getTypeNameParam(char * value)
  {
    if (value) {
      char * subvalue;
      switch (getsubopt(&value, s_typeOpts, &subvalue)) {
        case TYPE_INCREMENT:
        case TYPE_I:
          m_typeMask = 0x1 << TYPE_INCREMENT; return 0;
        case TYPE_AGGREGATE:
        case TYPE_A:
          m_typeMask = 0x1 << TYPE_AGGREGATE; return 0;
        case TYPE_DISPLAY:
        case TYPE_D:
          m_typeMask = 0x1 << TYPE_DISPLAY; return 0;
      }
    }
    std::cerr << "Unrecognized Bench Name option '" << optarg << "'!"
              << std::endl;
    std::cerr << "Try `" << m_progName << " -h' for more information."
              << std::endl;
    return -1;
  }

  /*
   */
  int getStrategyNameParam(char * value)
  {
    if (value) {
      char * subvalue;
      switch (getsubopt(&value, s_strategyOpts, &subvalue)) {
        case STRATEGY_TRAPEZOIDAL:
        case STRATEGY_T:
          m_strategyMask = 0x1 << STRATEGY_TRAPEZOIDAL; return 0;
        case STRATEGY_NAIVE:
        case STRATEGY_N:
          m_strategyMask = 0x1 << STRATEGY_NAIVE; return 0;
        case STRATEGY_WALK:
        case STRATEGY_W:
          m_strategyMask = 0x1 << STRATEGY_WALK; return 0;
      }
    }
    std::cerr << "Unrecognized Bench Name option '" << optarg << "'!"
              << std::endl;
    std::cerr << "Try `" << m_progName << " -h' for more information."
              << std::endl;
    return -1;
  }

private:
  static char * s_formatOpts[];
  static char * s_IOOpts[];
  static char * s_benchOpts[];
  static char * s_typeOpts[];
  static char * s_strategyOpts[];
  static char * s_boolOpts[];
    
private:
  static char s_optionStr[];
  const char * m_progName;

  std::string m_fullname;
  Dir_search m_dirs;

  int m_argc;
  char ** m_argv;
  const char * m_filename;
  bool m_verbose;
  unsigned int m_typeMask;
  unsigned int m_strategyMask;
  FormatId m_inputFormat;
  bool m_printHeader;
  int m_nameLength;
  int m_seconds;
  int m_samples;
  int m_iterations;
};

CGAL_END_NAMESPACE

#endif
