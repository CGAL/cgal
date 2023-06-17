
#ifndef COMMON_DEFS_H
#define COMMON_DEFS_H

#include <vector>

#include <qopenglfunctions_3_3_core.h>
//#include <qopenglfunctions_4_5_core.h>
#include <qvector3d.h>


using OpenGLFunctionsBase = QOpenGLFunctions_3_3_Core;


struct Line_strip_approx
{
  std::vector<QVector3D>  points;
  std::vector<GLuint>     offsets;
};


#endif
