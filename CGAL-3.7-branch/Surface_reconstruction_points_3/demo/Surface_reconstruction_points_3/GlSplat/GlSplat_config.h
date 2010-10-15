#ifndef GLSPLAT_CONFIG_H
#define GLSPLAT_CONFIG_H

#ifdef WIN32
  #ifdef gl_splat_EXPORTS
    #define GLSPLAT_EXPORT __declspec(dllexport)
  #else
    #define GLSPLAT_EXPORT __declspec(dllimport)
  #endif
#else // if not Windows
    #define GLSPLAT_EXPORT
#endif

#endif // GLSPLAT_CONFIG_H
