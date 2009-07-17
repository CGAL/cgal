// This file is part of GlSplat, a simple splatting C++ library
//
// Copyright (C) 2008-2009 Gael Guennebaud <g.gael@free.fr>
//
// GlSplat is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// GlSplat is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with GlSplat. If not, see <http://www.gnu.org/licenses/>.

#ifndef _GLSPLAT_SPLATRENDERER_H_
#define _GLSPLAT_SPLATRENDERER_H_

#include "GlSplat_config.h"
#include "Shader.h"
#include <QObject>
#include <QAction>
#include <QString>
#include <map>

class QGLFramebufferObject;
class QGLWidget;

namespace GlSplat {

/** \class SplatRenderer
  * \brief Helper class to render a set of points using a splatting alogirthm
  *
  * This class aims to render a set of oriented point with radius (called splats) using an
  * OpenGL based splating algorithm. This class is only responsible for the managing of the
  * OpenGL stats and shaders related to the splatting. The drawing of the geometry, i.e.,
  * sending the point data to the GPU has to be done by the user.
  *
  * Here is an example:
  * \code
  *  SplatRenderer renderer;
  *  renderer.init();
  *
  *  renderer.beginVisibilityPass();
  *  drawpoints();
  *  renderer.beginAttributePass();
  *  drawpoints();
  *  renderer.finalize();
  * \endcode
  *
  * Have a look at the demo to see a complete example based on QGLviewer.
  */
class GLSPLAT_EXPORT SplatRenderer
{
  bool mIsSupported;
  enum {
    DEFERRED_SHADING_BIT	= 0x000001,
    DEPTH_CORRECTION_BIT	= 0x000002,
    OUTPUT_DEPTH_BIT			= 0x000004,
    BACKFACE_SHADING_BIT	= 0x000008,
    FLOAT_BUFFER_BIT			= 0x000010
  };
  int mFlags;
  int mCachedFlags;
  int mRenderBufferMask;
  int mSupportedMask;

  int mCurrentPass;
  int mBindedPass;
  GLuint mDummyTexId; // on ATI graphics card we need to bind a texture to get point sprite working !
  bool mWorkaroundATI;
  bool mBuggedAtiBlending;
  bool mIsInitialized;
  GLuint mNormalTextureID;
  GLuint mDepthTextureID;
  Shader mShaders[3];
  QString mShaderSrcs[6];
  QGLFramebufferObject* mRenderBuffer;
  float mCachedMV[16];    // modelview matrix
  float mCachedProj[16];  // projection matrix
  GLint mCachedVP[4];     // viewport

  struct UniformParameters
  {
    float radiusScale;
    float preComputeRadius;
    float depthOffset;
    float oneOverEwaRadius;
    float halfVp[2];
    float rayCastParameter1[3];
    float rayCastParameter2[3];
    float depthParameterCast[2];

    void loadTo(Shader& prg);
    void update(float* mv, float* proj, GLint* vp);
  };

  UniformParameters mParams;

  QString loadSource(const QString& func,const QString& file);
  void configureShaders();
  void updateRenderBuffer();
  void enablePass(int n);

public:

  SplatRenderer();

  /** Must be called once an OpenGL context has been activated.
    * The main OpenGL context must be enabled, or, if you are using a QGLwiget,
    * you can pass it to this function without caring about the OpenGL context.
    */
  void init(QGLWidget *qglw = 0);

  /** \returns true is the hardware is supported
    * Must be called after init.
    */
  bool isSupported() { return mIsSupported; }

  /** Starts the first rendering pass
    * \returns false if an error occured
    */
  bool beginVisibilityPass();
  /** Starts the (optional) second rendering pass
    * \returns false if an error occured
    */
  bool beginAttributePass();
  /** Draw the rendered splats inside the main render target.
    * \returns false if an error occured
    */
  bool finalize();

  /** Sets a global scale factor for the splat radii
    * Default value is 1
    */
  void setRadiusScale(float v);

};

} // namepsace GlSplat

#endif // _GLSPLAT_SPLATRENDERER_H_

