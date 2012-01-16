// CORE LIBRARY FILE
#ifndef _CORE_GMP_H_
#define _CORE_GMP_H_

#include <CGAL/CORE/Impl.h>
#include <gmp.h>

namespace CORE { 

CGAL_CORE_EXPORT std::ostream& io_write (std::ostream &, mpz_srcptr);
CGAL_CORE_EXPORT std::ostream& io_write (std::ostream &, mpq_srcptr);
CGAL_CORE_EXPORT std::istream& io_read (std::istream &, mpz_ptr);
CGAL_CORE_EXPORT std::istream& io_read (std::istream &, mpq_ptr);
//std::ostream& operator<< (std::ostream &, mpz_srcptr);
//std::ostream& operator<< (std::ostream &, mpq_srcptr);
//std::istream& operator>> (std::istream &, mpz_ptr);
//std::istream& operator>> (std::istream &, mpq_ptr);

} //namespace CORE
#endif // _CORE_GMP_H_
