#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>

#include "Option_parser.hpp"

char * Option_parser::s_type_opts[] = {
  "increment", "aggregate", "pointLocation", "display", "i", "a", "l", "d"
};

char * Option_parser::s_strategy_opts[] = {
  "RIC", "naive", "walk", "simple", "triangle", "landmarks", 
  "R",   "n",     "w",    "s",      "t",        "l"
};

template <class MyId>
void Option_parser::my_validate(boost::any & v,
                                const std::vector<std::string> & values)
{
  typedef std::vector<MyId>     Vector_id;
  MyId tag;
  for (unsigned int i = 0; i < get_number_opts(tag); ++i) {
    if (compare_opt(i, values[0].c_str(), tag)) {
      if (v.empty()) {
        Vector_id vec;
        vec.push_back(MyId(i));
        v = boost::any(vec);
      } else {
        Vector_id vec = boost::any_cast<Vector_id>(v);
        vec.push_back(MyId(i));
        v = boost::any(vec);
      }
      return;
    }
  }
  throw po::validation_error("invalid value");
}

/* Overload the 'validate' function for the user-defined class */
void validate(boost::any & v, const std::vector<std::string> & values,
              Option_parser::Vector_type_id * target_type, int)
{
  Option_parser::my_validate<Option_parser::Type_id>(v, values);
}

/* Overload the 'validate' function for the user-defined class */
void validate(boost::any & v, const std::vector<std::string> & values,
              Option_parser::Vector_strategy_id * target_type, int)
{
  Option_parser::my_validate<Option_parser::Strategy_id>(v, values);
}

/*! Constructor */
Option_parser::Option_parser() :
  m_generic_opts("Generic options"),
  m_config_opts("Configuration options"),
  m_hidden_opts("Hidden options"),
  m_verbose_level(0),
  m_win_width(DEF_WIN_WIDTH),
  m_win_height(DEF_WIN_HEIGHT),
  m_type_mask(0xf),
  m_strategy_mask(0x3f),
  m_postscript(false),
  m_number_files(0)
{
  m_dirs.add(".");

  const char * root = getenv("ROOT");
  if (root) {
    m_dirs.add(std::string(root) + "/data/Segments_2");
    m_dirs.add(std::string(root) + "/data/Conics_2");
    m_dirs.add(std::string(root) + "/data/Polylines_2");
  }

  // Generic options:
  m_generic_opts.add_options()
    ("author,a", "print author name(s)")
    ("help,h", "print help message")
    ("license,l", "print licence information")
    ("version,v", "print version string")
    ;
  
  typedef std::vector<std::string> vs;
  
  // Options allowed on the command line, config file, or env. variables
  m_config_opts.add_options()
    ("input-path,P", po::value<vs>()->composing(), "input path")
    ("verbose,V", po::value<unsigned int>(&m_verbose_level)->default_value(0),
     "verbose level")
    ("type", po::value<std::vector<Type_id> >()->composing(),
     "Type\n"
     "\t\t\t\t\ti[ncrement]\t(0x1}\n"
     "\t\t\t\t\ta[ggregate]\t(0x2)\n"
     "\t\t\t\t\t[point_]l[ocation]\t(0x4)\n"
     "\t\t\t\t\td[isplay]\t(0x8)\n"
     )
    ("type-mask,T", po::value<unsigned int>(&m_type_mask)->default_value(0x3f),
     "type mask")
    ("strategy", po::value<std::vector<Strategy_id> >()->composing(),
     "Strategy\n"
     "\t\t\t\t\tR[IC]\t(0x1}\n"
     "\t\t\t\t\tn[aive]\t(0x2)\n"
     "\t\t\t\t\tw[alk]\t(0x4)\n"
     "\t\t\t\t\ts[imple]\t(0x8)\n"
     "\t\t\t\t\tt[riangle]\t(0x10)\n"
     "\t\t\t\t\tl[enmarks]\t(0x20)\n"
     )
    ("strategy-mask,S",
     po::value<unsigned int>(&m_strategy_mask)->default_value(0x3f),
     "strategy mask")
    ("width,W",
     po::value<unsigned int>(&m_win_width)->default_value(DEF_WIN_WIDTH ),
     "window width")
    ("height,H",
     po::value<unsigned int>(&m_win_height)->default_value(DEF_WIN_HEIGHT),
     "window height")
    ;
  
  // Options hidden to the user. Allowed only on the command line:
  m_hidden_opts.add_options()
    ("input-file", po::value<vs>()->composing(), "input file")
    ;

  m_visible_opts.add(m_generic_opts).add(m_bench_opts).add(m_config_opts);
  m_cmd_line_opts.add(m_generic_opts).add(m_bench_opts).add(m_config_opts).add(m_hidden_opts);
  m_config_file_opts.add(m_bench_opts).add(m_config_opts);
  m_environment_opts.add(m_bench_opts).add(m_config_opts);

  m_positional_opts.add("input-file", -1);  
}

/*! Parse the options */
void Option_parser::operator()(int argc, char * argv[])
{
  po::store(po::command_line_parser(argc, argv).
            options(m_cmd_line_opts).positional(m_positional_opts).run(),
            m_variable_map);
  
  std::ifstream ifs(".bench.cfg");
  po::store(parse_config_file(ifs, m_config_file_opts), m_variable_map);
  po::notify(m_variable_map);

  if (m_variable_map.count("help")) {
    std::cout << m_visible_opts << std::endl;
    throw Generic_option_exception(HELP);
    return;
  }

  if (m_variable_map.count("version")) {
    std::cout << "CGAL_VERSION" << std::endl;
    throw Generic_option_exception(VERSION);
    return;
  }

  if (m_variable_map.count("author")) {
    std::cout << "author Efi Fogel" << std::endl;
    throw Generic_option_exception(AUTHOR);
    return;
  }

  if (m_variable_map.count("license")) {
    std::cout << "license ???" << std::endl;
    throw Generic_option_exception(LICENSE);
    return;
  }

  if (m_variable_map.count("type")) {
    m_type_mask = 0;
    Vector_type_id types = m_variable_map["type"].as<Vector_type_id>();
    for (Vector_type_id_iter it = types.begin(); it != types.end(); ++it) {
      unsigned int size = sizeof(Option_parser::s_type_opts) / 8;
      m_type_mask |= 0x1 << ((*it).m_id % size);
    }
  }
  
  if (m_variable_map.count("strategy")) {
    Vector_strategy_id strategys =
      m_variable_map["strategy"].as<Vector_strategy_id>();
    for (Vector_strategy_id_iter it = strategys.begin();
         it != strategys.end(); ++it)
    {
      unsigned int size = sizeof(Option_parser::s_strategy_opts) / 8;
      m_strategy_mask |= 0x1 << ((*it).m_id % size);
    }
  }
  
  typedef std::vector<std::string> vs;
  if (m_variable_map.count("input-path")) {
    vs dirs = m_variable_map["input-path"].as<vs>();
    for (vs::iterator it = dirs.begin(); it != dirs.end(); ++it) {
      m_dirs.add((*it).c_str());
    }
  }
  
  if (!m_variable_map.count("input-file")) {
    std::string str("input file missing!");
    throw Input_file_missing_error(str);
    return;
  }

  vs files = m_variable_map["input-file"].as<vs>();
  m_full_names.resize(files.size());
  for (vs::iterator it = files.begin(); it != files.end(); ++it) {
    const char * file_name = (*it).c_str();
    if (!m_dirs.find(file_name, m_full_names[m_number_files++])) {
      std::cerr << "cannot find file " << file_name << "!" << std::endl;
      throw Error_exception(FILE_NOT_FOUND);
      return;
    }
  }
}

/*! Obtain the base file-name */
const std::string & Option_parser::get_file_name(unsigned int i) const
{
  typedef std::vector<std::string> vs;
  return m_variable_map["input-file"].as<vs>()[i];
}

/*! Obtain the full file-name */
const std::string & Option_parser::get_full_name(unsigned int i) const
{ return m_full_names[i]; }

/*! Obtain number of type options */
unsigned int Option_parser::get_number_opts(Type_id &)
{ return sizeof(s_type_opts) / sizeof(char *); }

/*! Obtain number of strategy options */
unsigned int Option_parser::get_number_opts(Strategy_id &)
{ return sizeof(s_strategy_opts) / sizeof(char *); }
