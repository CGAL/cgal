#include <CGAL/Kinetic/basic.h>

CGAL_KINETIC_BEGIN_NAMESPACE
namespace internal {
  extern unsigned int zero_certificates__;
  extern unsigned int function_degeneracies__;
  extern unsigned int io_errors__;
  extern unsigned int audit_failures__;

  void write_debug_counters(std::ostream &out);
}
CGAL_KINETIC_END_NAMESPACE
