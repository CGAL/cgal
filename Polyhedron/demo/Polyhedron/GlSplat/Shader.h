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

#ifndef _GLSPLAT_Shader_h_
#define _GLSPLAT_Shader_h_
#include <string>
#include <map>
#include <QOpenGLFunctions_3_3_Core>
#include "Viewer_interface.h"

namespace GlSplat {

/** Permet de manipuler des shaders en GLSL (OpenGL2.0)
    Exemple d'utilisation:
    \code
    // shader creation:
    Shader* myShader = new Shader();
    // loading from files (compilation + linking):
    myShader->loadFromFiles("myShaderFile.vtx", "myShaderFile.frg");

    // ...

    // at rending time:
    myShader->enable();
    // draw objects
    myShader->disable();
    \endcode
*/

class Shader
{
public:
    Shader(void)
      : mIsValid(false)
    { }
    ~Shader()
    {
    }
    void setViewer(Viewer_interface *);
    /** add a \#define
    */
    void define(const char* name, const char* value);

    /** Compiles and links the shader from 2 source files
        \param fileV vertex shader ("" if no vertex shader)
        \param fileF fragment shader ("" if no fragment shader)
        \return true if no error occurs
    */
//     bool loadFromFiles(const std::string& fileV, const std::string& fileF);

    bool loadSources(const char* vsrc, const char* fsrc);

    /** Enable the shader
    */
    void activate();

    /** Releases the shader
    */
    void release() ;

    /** \return the index of the uniform variable \a name
    */
    int getUniformLocation(const char* name);

    /** Forces a sampler to a given unit
        Example:
        \code
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE2D, myTextureID);
    myShader->setSamplerUnit("mySampler", 2);
        \endcode
    */
    void setSamplerUnit(const char* samplerName, int textureUnit);

    /** \returns the index of the generic attribute \a name
        Tp be used with glVertexAttrib*(...) ou glVertexAttribPointer(...)
        Example:
        \code
    int tangentAttribID = myShader->getAttribLocation("tangent");
    Vector3 tangent(...);
    glVertexAttrib3fv(tangentAttribID, tangent);
    // ou
    Vector3* tangents = new Vector3[...];
    glVertexAttribPointer(tangentAttribID, 3, GL_FLOAT, GL_FALSE, 0, tangents);
    glEnableVertexAttribArray(tangentAttribID);
        \endcode
    */
    int getAttribLocation(const char* name);

    inline void setUniform(const char* name, float a)
    { viewer->glUniform1f(viewer->glGetUniformLocation(mProgramID, name), a); }

    inline void setUniform(const char* name, int a)
    { viewer->glUniform1i(viewer->glGetUniformLocation(mProgramID, name), a); }

    inline void setUniform2(const char* name, float* a)
    { viewer->glUniform2fv(viewer->glGetUniformLocation(mProgramID, name), 1, a); }

    inline void setUniform3(const char* name, float* a)
    { viewer->glUniform3fv(viewer->glGetUniformLocation(mProgramID, name), 1, a); }

    inline void setUniform4(const char* name, float* a)
    { viewer->glUniform4fv(viewer->glGetUniformLocation(mProgramID, name), 1, a); }

    inline void setUniform(const char* name, float a, float b)
    { viewer->glUniform2f(viewer->glGetUniformLocation(mProgramID, name), a, b); }

    inline void setUniform(const char* name, float a, float b, float c)
    { viewer->glUniform3f(viewer->glGetUniformLocation(mProgramID, name), a, b, c); }

    inline void setUniform(const char* name, float a, float b, float c, float d)
    { viewer->glUniform4f(viewer->glGetUniformLocation(mProgramID, name), a, b, c, d); }

     Viewer_interface *viewer;
protected:

    bool mIsValid;
    typedef std::map<std::string,std::string> DefineMap;
    DefineMap mDefines;
    void printInfoLog(GLuint objectID);
    GLuint mProgramID;
};

} // namepsace GlSplat

#endif
