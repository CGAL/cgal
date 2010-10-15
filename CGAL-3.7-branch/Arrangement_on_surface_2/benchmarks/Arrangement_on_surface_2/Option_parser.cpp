#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/exception.hpp>

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
  m_dirs.push_front(".");

  const char * root_c = getenv("ROOT");
  if (root_c) {
    fi::path root(root_c, fi::native);
    fi::path dir;
    dir = root / "data/curves/segments";
    m_dirs.push_back(dir);
    dir = root / "data/curves/conics";
    m_dirs.push_back(dir);
    dir = root / "data/curves/polylines";
    m_dirs.push_back(dir);
  }

  // Generic options:
  m_generic_opts.add_options()
    ("author,a", "print author name(s)")
    ("help,h", "print help message")
    ("license,l", "print licence information")
    ("version,v", "print version string")
    ;
    
  // Options allowed on the command line, config file, or env. variables
  m_config_opts.add_options()
    ("input-path,P", po::value<Input_path>()->composing(), "input path")
    ("verbose,V", po::value<unsigned int>(&m_verbose_level)->default_value(0),
     "verbose level")
    ("type", po::value<std::vector<Type_id> >()->composing(),
     "Type\n"
     "  i[ncrement]        (0x1}\n"
     "  a[ggregate]        (0x2)\n"
     "  [point_]l[ocation] (0x4)\n"
     "  d[isplay]          (0x8)\n"
     )
    ("type-mask,T", po::value<unsigned int>(&m_type_mask)->default_value(0x3f),
     "type mask")
    ("strategy", po::value<std::vector<Strategy_id> >()->composing(),
     "Strategy\n"
     "  R[IC]      (0x1}\n"
     "  n[aive]    (0x2)\n"
     "  w[alk]     (0x4)\n"
     "  s[imple]   (0x8)\n"
     "  t[riangle] (0x10)\n"
     "  l[enmarks] (0x20)\n"
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
    ("input-file", po::value<Input_path>()->composing(), "input file")
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

  // Add directories specified on the command line via the "input-path" opt.:
  Add_dir tmp = for_each_dir(Add_dir(m_dirs));
  
  if (!m_variable_map.count("input-file")) {
    std::string str("input file missing!");
    throw Input_file_missing_error(str);
    return;
  }

  Input_path files = m_variable_map["input-file"].as<Input_path>();
  m_full_names.resize(files.size());
  Input_path_const_iterator it;
  for (it = files.begin(); it != files.end(); ++it) {
    fi::path file_path(*it, fi::native);
    if (file_path.is_complete()) {
      if (fi::exists(file_path))
        m_full_names[m_number_files] = file_path.native_file_string();
    } else {
      for (Path_iter pi = m_dirs.begin(); pi != m_dirs.end(); ++pi) {
        fi::path full_file_path = *pi / file_path;
        if (!fi::exists(full_file_path)) continue;
        m_full_names[m_number_files] = full_file_path.native_file_string();
        break;
      }
    }
    if (m_full_names[m_number_files].empty()) {      
      std::cerr << "cannot find file " << (*it).c_str() << "!" << std::endl;
      throw Error_exception(FILE_NOT_FOUND);
      return;
    }
    ++m_number_files;
  }
}

/*! Obtain the base file-name */
const std::string & Option_parser::get_file_name(unsigned int i) const
{
  return m_variable_map["input-file"].as<Input_path>()[i];
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
