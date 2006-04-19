#include <dsrpdb_internal/Error_logger.h>
#include <cassert>
#include <cstdlib>

namespace dsrpdb_internal {
  Error_logger error_logger;

  void Error_logger::new_warning(const char *str) {
    if (enabled_) {
      if (warnings_.find(str) == warnings_.end()) {
	for (std::map<std::string, int>::const_iterator it =warnings_.begin(); 
	     it != warnings_.end(); ++it){
	  assert(it->first != str);
	}
	warnings_[str]=0;
      }
      ++warnings_[str];
    }
  }

  void Error_logger::new_fatal_error(const char *err) {
    std::cerr << "DSRPDB fatal error: " << err << std::endl;
    assert(0);
    exit(EXIT_FAILURE);
  }

  void Error_logger::new_internal_error(const char* err) {
    std::cerr << "DSRPDB internal error: " << err << std::endl;
    std::cerr << "Please report this to the author (and provide a PDB)." << std::endl;
    assert(0);
    exit(EXIT_FAILURE);
  }

  void Error_logger::dump() {
    if (enabled_) {
      if (!context_.empty()) {
	std::cerr << "In PDB file " << context_ << ":\n";
      }
      for (std::map<std::string, int>::const_iterator it =warnings_.begin(); 
	   it != warnings_.end(); ++it){
	if (it->second==1) {
	  std::cerr << "DSRPDB Warning: " << it->first << std::endl;
	} else {
	  std::cerr << "DSRPDB " << it->second << " occurences of Warning: " << it->first << std::endl;
	}
      }
    }
    warnings_.clear();
  }
}
