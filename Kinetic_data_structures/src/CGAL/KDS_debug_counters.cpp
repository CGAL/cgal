#include <CGAL/Kinetic/internal/debug_counters.h>

#include <iostream>

CGAL_KINETIC_BEGIN_NAMESPACE
namespace internal {
  unsigned int function_degeneracies__=0;
  unsigned int io_errors__=0;
  unsigned int audit_failures__=0;

  void write_debug_counters(std::ostream &out) {
    out << "Degeneracies " << function_degeneracies__ << std::endl;
    if (io_errors__ != 0) out << "I/O errors " << io_errors__ << std::endl;
    if (audit_failures__ != 0) out << "Audit failures " << audit_failures__ << std::endl;
  }
}
CGAL_KINETIC_END_NAMESPACE
