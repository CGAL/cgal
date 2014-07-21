#ifndef GLSPLAT_CONFIG_H
#define GLSPLAT_CONFIG_H

#include <CGAL/export/helpers.h>

#ifdef gl_splat_EXPORTS
  #define GLSPLAT_EXPORT CGAL_DLL_EXPORT
#else
  #define GLSPLAT_EXPORT CGAL_DLL_IMPORT
#endif

#endif // GLSPLAT_CONFIG_H
