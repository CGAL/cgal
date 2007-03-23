// CORE LIBRARY FILE
#ifndef _CORE_GMP_H_
#define _CORE_GMP_H_

#include <CGAL/CORE/Impl.h>
#include <gmp.h>

CORE_BEGIN_NAMESPACE

std::ostream& io_write (std::ostream &, mpz_srcptr);
std::ostream& io_write (std::ostream &, mpq_srcptr);
std::istream& io_read (std::istream &, mpz_ptr);
std::istream& io_read (std::istream &, mpq_ptr);
//std::ostream& operator<< (std::ostream &, mpz_srcptr);
//std::ostream& operator<< (std::ostream &, mpq_srcptr);
//std::istream& operator>> (std::istream &, mpz_ptr);
//std::istream& operator>> (std::istream &, mpq_ptr);

CORE_END_NAMESPACE
#endif // _CORE_GMP_H_
