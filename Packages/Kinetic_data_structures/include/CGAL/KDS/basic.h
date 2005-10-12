#ifndef CGAL_KDS_BASIC_H
#define CGAL_KDS_BASIC_H

/*!
  \file CGAL/KDS/basic.h The header which defines the standard things for the KDS namespace and module.

  \todo cgal documentation

  \namespace CGAL The CGAL project namespace

  \namespace CGAL::KDS The namespace for classes which involve the Kinetic Data Structures framework.

*/


#include <CGAL/basic.h>

#include <CGAL/KDS/internal/config.h>

#ifdef CGAL_CHECK_EXPENSIVE
#ifndef CGAL_KDS_CHECK_EXPENSIVE
#define CGAL_KDS_CHECK_EXPENSIVE
#endif
#endif

#ifdef CGAL_CHECK_EXACTNESS
#ifndef CGAL_KDS_CHECK_EXACTNESS
#define CGAL_KDS_CHECK_EXACTNESS
#endif
#endif

#define CGAL_KDS_BEGIN_NAMESPACE CGAL_BEGIN_NAMESPACE \
namespace KDS {

#define CGAL_KDS_END_NAMESPACE } \
CGAL_END_NAMESPACE 

#define CGAL_KDS_BEGIN_INTERNAL_NAMESPACE CGAL_KDS_BEGIN_NAMESPACE \
namespace internal {

#define CGAL_KDS_END_INTERNAL_NAMESPACE } \
CGAL_KDS_END_NAMESPACE 



CGAL_KDS_BEGIN_NAMESPACE
//! The types of logs available.
typedef enum {LOG_NONE=0, LOG_SOME=2, LOG_LOTS=3} Log_level;

CGAL_KDS_END_NAMESPACE


#include <CGAL/KDS/internal/Log.h>

#define CGAL_KDS_LOG(level, expr) if (CGAL::KDS::internal::Logs::get().is_output(level))\
 { CGAL::KDS::internal::Logs::get().stream(level) << expr;};
#define CGAL_KDS_LOG_WRITE(level, expr) if (CGAL::KDS::internal::Logs::get().is_output(level))\
{std::ostream &LOG_STREAM= CGAL::KDS::internal::Logs::get().stream(level); expr;}
#define CGAL_KDS_ERROR(expr) std::cerr << expr;
#define CGAL_KDS_SET_LOG_LEVEL(level) CGAL::KDS::internal::Logs::get().set_level(level);



#endif
