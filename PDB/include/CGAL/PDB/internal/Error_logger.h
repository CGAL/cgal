#ifndef CGAL_DSR_PDB_ERROR_LOGGER_H
#define CGAL_DSR_PDB_ERROR_LOGGER_H
#include <CGAL/PDB/basic.h>
#include <iostream>
#include <map>
#include <string>
#include <CGAL/Tools/Log.h>

namespace CGAL { namespace PDB { namespace internal {
  class Error_logger {
  public:
    Error_logger():enabled_(true){}
    ~Error_logger(){dump();}
    
    bool is_output() const {return enabled_;}
    void set_is_output(bool tf) { enabled_=tf;}

    void set_context(const char *c) {
      context_=c;
    }

    void new_warning(const char *c);

    void new_fatal_error(const char *c);

    void new_internal_error(const char *c);

    void dump(); 
    
    std::map<std::string, int> warnings_;
    bool enabled_;
    std::string context_;
  };

  extern Error_logger error_logger;

}}}
#endif
