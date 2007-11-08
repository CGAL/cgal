#include <CGAL/PDB/internal/Error_logger.h>
#include <cstdlib>
#include <cassert>
CGAL_PDB_BEGIN_INTERNAL_NAMESPACE 

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
    CGAL_ERROR("DSRPDB fatal error: " << err);
    CGAL_error();
    //exit(EXIT_FAILURE);
  }

  void Error_logger::new_internal_error(const char* err) {
    CGAL_ERROR( "DSRPDB internal error: " << err);
    CGAL_ERROR( "Please report this to the author (and provide a PDB).");
    CGAL_error();
    //exit(EXIT_FAILURE);
  }

  void Error_logger::dump() {
    if (enabled_) {
      if (!context_.empty()) {
	CGAL_LOG(Log::SOME, "In PDB file " << context_ << ":\n");
      }
      for (std::map<std::string, int>::const_iterator it =warnings_.begin(); 
	   it != warnings_.end(); ++it){
	if (it->second==1) {
	  CGAL_LOG(Log::SOME, "DSRPDB Warning: " << it->first << std::endl);
	} else {
	  CGAL_LOG(Log::SOME, "DSRPDB " << it->second 
		   << " occurences of Warning: " << it->first << std::endl);
	}
      }
    }
    warnings_.clear();
  }
CGAL_PDB_END_INTERNAL_NAMESPACE
