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

#include <QtGui>

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <GL/glew.h>

#include "GlSplat.h"
#include "Shader.h"

#include <QGLWidget>
#include <QTextStream>
#include <QGLFramebufferObject>


namespace GlSplat {

SplatRenderer::SplatRenderer()
{
  mNormalTextureID = 0;
  mDepthTextureID = 0;
  mIsSupported = false;
  mRenderBuffer = 0;
  mWorkaroundATI = false;
  mBuggedAtiBlending = false;
  mDummyTexId = 0;
  mIsInitialized = false;

  mFlags = DEFERRED_SHADING_BIT | DEPTH_CORRECTION_BIT | FLOAT_BUFFER_BIT | OUTPUT_DEPTH_BIT;
  mCachedFlags = ~mFlags;
  // union of bits which controls the render buffer
  mRenderBufferMask = DEFERRED_SHADING_BIT | FLOAT_BUFFER_BIT;
}

QString SplatRenderer::loadSource(const QString& func,const QString& filename)
{
  QString res;
  QFile f(":/SplatRenderer/shaders/" + filename);
  if (!f.open(QFile::ReadOnly))
  {
    std::cerr << "failed to load shader file " << filename.toAscii().data() << "\n";
    return res;
  }
  else qDebug("Succesfully loaded shader func '%s' in file '%s'",qPrintable(func),qPrintable(filename));
  QTextStream stream(&f);
  res = stream.readAll();
  f.close();
  res = QString("#define __%1__ 1\n").arg(func)
      + QString("#define %1 main\n").arg(func)
      + res;
  return res;
}

void SplatRenderer::configureShaders()
{
  //  const char* passNames[3] = {"Visibility","Attribute","Finalization"};
  QString defines = "";
  if (mFlags & DEFERRED_SHADING_BIT)
    defines += "#define ES_DEFERRED_SHADING\n";
  if (mFlags & DEPTH_CORRECTION_BIT)
    defines += "#define ES_DEPTH_CORRECTION\n";
  if (mFlags & OUTPUT_DEPTH_BIT)
    defines += "#define ES_OUTPUT_DEPTH 1\n";
  if (mFlags & BACKFACE_SHADING_BIT)
    defines += "#define ES_BACKFACE_SHADING\n";
  if (mWorkaroundATI)
    defines += "#define ES_ATI_WORKAROUND\n";

  QString shading =
"vec4 meshlabLighting(vec4 color, vec3 eyePos, vec3 normal)"
"{"
"	normal = normalize(normal);"
"	vec3 lightVec = normalize(gl_LightSource[0].position.xyz);"
"	vec3 halfVec = normalize( lightVec - normalize(eyePos) );"
"	float aux_dot = dot(normal,lightVec);"
"	float diffuseCoeff = clamp(aux_dot, 0.0, 1.0);"
" float specularCoeff = aux_dot>0.0 ? clamp(pow(clamp(dot(halfVec, normal),0.0,1.0),gl_FrontMaterial.shininess), 0.0, 1.0) : 0.0;"
"	return vec4(color.rgb * ( gl_FrontLightProduct[0].ambient.rgb + diffuseCoeff * gl_FrontLightProduct[0].diffuse.rgb) + specularCoeff * gl_FrontLightProduct[0].specular.rgb, 1.0);"
"}\n";

  for (int k=0;k<3;++k)
  {
    QString vsrc = shading + defines + mShaderSrcs[k*2+0];
    QString fsrc = shading + defines + mShaderSrcs[k*2+1];
    if(!mShaders[k].loadSources(mShaderSrcs[k*2+0]!="" ? vsrc.toAscii().data() : 0,
                                mShaderSrcs[k*2+1]!="" ? fsrc.toAscii().data() : 0/*,
                                Shader::Warnings*/))
      mIsSupported = false;
  }
}

void SplatRenderer::init(QGLWidget *qglw)
{
  mIsSupported = true;
  if(qglw)
    qglw->makeCurrent();
  glewInit();

  const char* rs = (const char*)glGetString(GL_RENDERER);
  QString rendererString("");
  if(rs)
    rendererString = QString(rs);
  mWorkaroundATI = rendererString.startsWith("ATI") || rendererString.startsWith("AMD");
  // FIXME: maybe some recent HW correctly supports floating point blending...
  mBuggedAtiBlending = rendererString.startsWith("ATI") || rendererString.startsWith("AMD");

  if (mWorkaroundATI && mDummyTexId==0)
  {
    glActiveTexture(GL_TEXTURE0);
    glGenTextures(1,&mDummyTexId);
    glBindTexture(GL_TEXTURE_2D, mDummyTexId);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE, 4, 4, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, 0);
  }

  // let's check the GPU capabilities
  mSupportedMask = DEPTH_CORRECTION_BIT | BACKFACE_SHADING_BIT;
  if (!QGLFramebufferObject::hasOpenGLFramebufferObjects ())
  {
    std::cout << "SplatRenderer: error OpenGL frame buffer objects are not supported. (please, try to update your drivers)\n";
    mIsSupported = false;
    return;
  }
  if (GLEW_ARB_texture_float)
    mSupportedMask |= FLOAT_BUFFER_BIT;
  else
    std::cout << "SplatRenderer: warning floating point textures are not supported.\n";

  if (GLEW_ARB_draw_buffers && (!mBuggedAtiBlending))
    mSupportedMask |= DEFERRED_SHADING_BIT;
  else
    std::cout << "SplatRenderer: warning deferred shading is not supported.\n";

  if (GLEW_ARB_shadow)
    mSupportedMask |= OUTPUT_DEPTH_BIT;
  else
    std::cerr << "SplatRenderer: warning copy of the depth buffer is not supported.\n";

  mFlags = mFlags & mSupportedMask;

  // load shader source
  mShaderSrcs[0] = loadSource("VisibilityVP","Raycasting.glsl");
  mShaderSrcs[1] = loadSource("VisibilityFP","Raycasting.glsl");
  mShaderSrcs[2] = loadSource("AttributeVP","Raycasting.glsl");
  mShaderSrcs[3] = loadSource("AttributeFP","Raycasting.glsl");
  mShaderSrcs[4] = "";
  mShaderSrcs[5] = loadSource("Finalization","Finalization.glsl");

  mCurrentPass = 2;
  mBindedPass = -1;
  mIsInitialized = true;
  GL_TEST_ERR
}

void SplatRenderer::updateRenderBuffer()
{
  if ( (!mRenderBuffer)
    || (mRenderBuffer->width()!=mCachedVP[2])
    || (mRenderBuffer->height()!=mCachedVP[3])
    || ( (mCachedFlags & mRenderBufferMask) != (mFlags & mRenderBufferMask) ))
  {
    delete mRenderBuffer;
    GLenum fmt = (mFlags&FLOAT_BUFFER_BIT) ? GL_RGBA16F_ARB : GL_RGBA;
    mRenderBuffer = new QGLFramebufferObject(mCachedVP[2], mCachedVP[3],
        (mFlags&OUTPUT_DEPTH_BIT) ? QGLFramebufferObject::NoAttachment : QGLFramebufferObject::Depth,
        GL_TEXTURE_RECTANGLE_ARB, fmt);

    if (!mRenderBuffer->isValid())
    {
      std::cout << "SplatRenderer: invalid FBO\n";
    }

    GL_TEST_ERR
    if (mFlags&DEFERRED_SHADING_BIT)
    {
      // in deferred shading mode we need an additional buffer to accumulate the normals
      if (mNormalTextureID==0)
        glGenTextures(1,&mNormalTextureID);
      glBindTexture(GL_TEXTURE_RECTANGLE_ARB, mNormalTextureID);
      glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, fmt, mCachedVP[2], mCachedVP[3], 0, GL_RGBA, GL_FLOAT, 0);
      glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      mRenderBuffer->bind();
      glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT1_EXT, GL_TEXTURE_RECTANGLE_ARB, mNormalTextureID, 0);
      mRenderBuffer->release();
      GL_TEST_ERR
    }

    if (mFlags&OUTPUT_DEPTH_BIT)
    {
      // to output the depth values to the final depth buffer we need to
      // attach a depth buffer as a texture
      if (mDepthTextureID==0)
        glGenTextures(1,&mDepthTextureID);
      glBindTexture(GL_TEXTURE_RECTANGLE_ARB, mDepthTextureID);
      glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, GL_DEPTH_COMPONENT24_ARB, mCachedVP[2], mCachedVP[3], 0, GL_DEPTH_COMPONENT, GL_FLOAT, 0);
      glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      mRenderBuffer->bind();
      glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_TEXTURE_RECTANGLE_ARB, mDepthTextureID, 0);
      mRenderBuffer->release();
      GL_TEST_ERR
    }
  }
}

bool SplatRenderer::beginVisibilityPass()
{
  if (!mIsInitialized)
  {
    init();
  }
  if (!isSupported())
  {
    std::cerr << "SplatRenderer error: not supported hardware\n";
    return false;
  }
  if (mCurrentPass!=2)
  {
    std::cerr << "SplatRenderer error: programming error when calling beginVisibilityPass\n";
    return false;
  }

  glPushAttrib(GL_ALL_ATTRIB_BITS);

  mCurrentPass = 0;

  // grab projection info
  glGetIntegerv(GL_VIEWPORT, mCachedVP);
  glGetFloatv(GL_MODELVIEW_MATRIX, mCachedMV);
  glGetFloatv(GL_PROJECTION_MATRIX, mCachedProj);

  updateRenderBuffer();
  if (mCachedFlags != mFlags)
    configureShaders();

  // configureShaders may detect that shaders are actually not supported.
  if (!isSupported())
  {
    std::cerr << "SplatRenderer error: not supported hardware\n";
    return false;
  }

  mCachedFlags = mFlags;

  mParams.update(mCachedMV, mCachedProj, mCachedVP);
  mParams.loadTo(mShaders[mCurrentPass]);

  mRenderBuffer->bind();
  if (mFlags & DEFERRED_SHADING_BIT)
  {
    GLenum buf[2] = {GL_COLOR_ATTACHMENT0_EXT,GL_COLOR_ATTACHMENT1_EXT};
    glDrawBuffersARB(2, buf);
  }
  glViewport(mCachedVP[0],mCachedVP[1],mCachedVP[2],mCachedVP[3]);
  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  enablePass(mCurrentPass);
  GL_TEST_ERR;
  return true;
}
bool SplatRenderer::beginAttributePass()
{
  if (!isSupported())
  {
    std::cerr << "SplatRenderer error: not supported hardware\n";
    return false;
  }
  if (mCurrentPass!=0)
  {
    std::cerr << "SplatRenderer error: programming error when calling beginAttributePass (must be called after the visiblity pass)\n";
    return false;
  }

  mCurrentPass = 1;
  mParams.loadTo(mShaders[mCurrentPass]);
  enablePass(mCurrentPass);
  GL_TEST_ERR;
  return true;
}
bool SplatRenderer::finalize()
{
  if (!isSupported())
  {
    std::cerr << "SplatRenderer error: not supported hardware\n";
    return false;
  }

  // this is the last pass: normalization by the sum of weights + deferred shading
  mShaders[mCurrentPass].release();
  mRenderBuffer->release();

  if ( (mCurrentPass!=0) && (mCurrentPass!=1))
  {
    std::cerr << "SplatRenderer error: programming error when calling finalize (must be called after the visiblity or attribute pass)\n";
    return false;
  }

  mCurrentPass = 2;

  if (mFlags&DEFERRED_SHADING_BIT)
    glDrawBuffer(GL_BACK);

  enablePass(mCurrentPass);GL_TEST_ERR

  // switch to normalized 2D rendering mode
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();


  GL_TEST_ERR
  mShaders[2].setUniform("viewport",float(mCachedVP[0]),float(mCachedVP[1]),float(mCachedVP[2]),float(mCachedVP[3]));GL_TEST_ERR
  mShaders[2].setUniform("ColorWeight",0); GL_TEST_ERR // this is a texture unit
  glActiveTexture(GL_TEXTURE0);GL_TEST_ERR
  glBindTexture(GL_TEXTURE_RECTANGLE_ARB,mRenderBuffer->texture());GL_TEST_ERR

  if (mFlags&DEFERRED_SHADING_BIT)
  {
    mShaders[2].setUniform("unproj", mCachedProj[10], mCachedProj[14]);GL_TEST_ERR
    mShaders[2].setUniform("NormalWeight",1);GL_TEST_ERR // this is a texture unit
    glActiveTexture(GL_TEXTURE1);GL_TEST_ERR
    glBindTexture(GL_TEXTURE_RECTANGLE_ARB,mNormalTextureID);GL_TEST_ERR
    GL_TEST_ERR
  }

  if (mFlags&OUTPUT_DEPTH_BIT)
  {
    mShaders[2].setUniform("Depth",2);GL_TEST_ERR // this is a texture unit
    glActiveTexture(GL_TEXTURE2);GL_TEST_ERR
    glBindTexture(GL_TEXTURE_RECTANGLE_ARB,mDepthTextureID);GL_TEST_ERR
    GL_TEST_ERR
  }
  else
  {
    glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);
  }

  // draw a quad covering the whole screen
  float viewVec[] = {1.f/mCachedProj[0], 1.f/mCachedProj[5], -1};

  glBegin(GL_QUADS);
    glColor3f(1, 0, 0);
    glTexCoord3f(viewVec[0],viewVec[1],viewVec[2]);
    glMultiTexCoord2f(GL_TEXTURE1,1.,1.);
    glVertex3f(1,1,0);

    glColor3f(1, 1, 0);
    glTexCoord3f(-viewVec[0],viewVec[1],viewVec[2]);
    glMultiTexCoord2f(GL_TEXTURE1,0.,1.);
    glVertex3f(-1,1,0);

    glColor3f(0, 1, 1);
    glTexCoord3f(-viewVec[0],-viewVec[1],viewVec[2]);
    glMultiTexCoord2f(GL_TEXTURE1,0.,0.);
    glVertex3f(-1,-1,0);

    glColor3f(1, 0, 1);
    glTexCoord3f(viewVec[0],-viewVec[1],viewVec[2]);
    glMultiTexCoord2f(GL_TEXTURE1,1.,0.);
    glVertex3f(1,-1,0);
  glEnd();
  if (!(mFlags&OUTPUT_DEPTH_BIT))
  {
      glEnable(GL_DEPTH_TEST);
      glDepthMask(GL_TRUE);
  }

  mShaders[mCurrentPass].release();

  // restore matrices
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();

  glPopAttrib();
  return true;
}

void SplatRenderer::enablePass(int n)
{
  if (!isSupported())
  {
    return;
  }
  if (mBindedPass!=n)
  {
    if (mBindedPass>=0)
      mShaders[mBindedPass].release();
    mShaders[n].activate();
    mBindedPass = n;

    // set GL states
    if (n==0)
    {
      glDisable(GL_LIGHTING);
// 			glDisable(GL_POINT_SMOOTH);
      glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

      glAlphaFunc(GL_LESS,1);
      glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_FALSE);
      glDepthMask(GL_TRUE);
      glDisable(GL_BLEND);
      glEnable(GL_ALPHA_TEST);
      glEnable(GL_DEPTH_TEST);

// 			glActiveTexture(GL_TEXTURE0);
// 			glTexEnvf(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
// 			glEnable(GL_POINT_SPRITE_ARB);
    }
    if (n==1)
    {
      glDisable(GL_LIGHTING);
      glEnable(GL_POINT_SMOOTH);
      glActiveTexture(GL_TEXTURE0);
      glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

      glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
      glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE, GL_ONE,GL_ONE);
// 			//glBlendFuncSeparate(GL_ONE, GL_ZERO, GL_ONE,GL_ZERO);
// 			glBlendFunc(GL_ONE,GL_ZERO);
      glDepthMask(GL_FALSE);
      glEnable(GL_BLEND);
      glEnable(GL_DEPTH_TEST);
      glDisable(GL_ALPHA_TEST);

// 			glActiveTexture(GL_TEXTURE0);

    }
    if ( (n==0) || (n==1) )
    {
      // enable point sprite rendering mode
      glActiveTexture(GL_TEXTURE0);
      if (mWorkaroundATI)
      {
        glBindTexture(GL_TEXTURE_2D, mDummyTexId);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE, 2, 2, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, 0);
        glPointParameterf(GL_POINT_SPRITE_COORD_ORIGIN, GL_LOWER_LEFT);
        // hm... ^^^^
      }
      glTexEnvf(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
      glEnable(GL_POINT_SPRITE_ARB);
    }
    if (n==2)
    {
      glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
      glDepthMask(GL_TRUE);
      glDisable(GL_LIGHTING);
      glDisable(GL_BLEND);
    }
  }
}

void SplatRenderer::UniformParameters::update(float* mv, float* proj, GLint* vp)
{
  // extract the uniform scale
  float scale = sqrtf(mv[0]*mv[0]+mv[1]*mv[1]+mv[2]*mv[2]);

  radiusScale = scale;
  preComputeRadius = - (std::max)(proj[0]*vp[2], proj[5]*vp[3]);
  depthOffset = 2.0;
  oneOverEwaRadius = 0.70710678118654;
  halfVp[0] = 0.5*vp[2];
  halfVp[1] = 0.5*vp[3];
  rayCastParameter1[0] = 2./(proj[0]*vp[2]);
  rayCastParameter1[1] = 2./(proj[5]*vp[3]);
  rayCastParameter1[2] = 0.0;
  rayCastParameter2[0] = -1./proj[0];
  rayCastParameter2[1] = -1./proj[5];
  rayCastParameter2[2] = -1.0;
  depthParameterCast[0] = 0.5*proj[14];
  depthParameterCast[1] = 0.5-0.5*proj[10];
}

void SplatRenderer::UniformParameters::loadTo(Shader& prg)
{
  prg.activate();
  prg.setUniform("expeRadiusScale",       radiusScale);
  prg.setUniform("expePreComputeRadius",  preComputeRadius);
  prg.setUniform("expeDepthOffset",       depthOffset);
  prg.setUniform("oneOverEwaRadius",      oneOverEwaRadius);
  prg.setUniform2("halfVp",               halfVp);
  prg.setUniform3("rayCastParameter1",    rayCastParameter1);
  prg.setUniform3("rayCastParameter2",    rayCastParameter2);
  prg.setUniform2("depthParameterCast",   depthParameterCast);
}

void SplatRenderer::setRadiusScale(float v)
{
  mParams.radiusScale = v;
}

} // namepsace GlSplat

