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

#include <GL/glew.h>
#include <CGAL/glu.h>

#ifndef NDEBUG
    #define GL_TEST_ERR\
        {\
            GLenum eCode;\
            if((eCode=glGetError())!=GL_NO_ERROR)\
                std::cerr << "OpenGL error : " <<  gluErrorString(eCode) << " in " <<  __FILE__ << " : " << __LINE__ << std::endl;\
        }
#else
    #define GL_TEST_ERR
#endif

#include <string>
#include <map>

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
    void activate() const;

    /** Releases the shader
    */
    void release() const;

    /** \return the index of the uniform variable \a name
    */
    int getUniformLocation(const char* name) const;

    /** Forces a sampler to a given unit
        Example:
        \code
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE2D, myTextureID);
    myShader->setSamplerUnit("mySampler", 2);
        \endcode
    */
    void setSamplerUnit(const char* samplerName, int textureUnit) const;

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
    int getAttribLocation(const char* name) const;

    inline void setUniform(const char* name, float a) const
    { glUniform1f(glGetUniformLocation(mProgramID, name), a); }

    inline void setUniform(const char* name, int a) const
    { glUniform1i(glGetUniformLocation(mProgramID, name), a); }

    inline void setUniform2(const char* name, float* a) const
    { glUniform2fv(glGetUniformLocation(mProgramID, name), 1, a); }

    inline void setUniform3(const char* name, float* a) const
    { glUniform3fv(glGetUniformLocation(mProgramID, name), 1, a); }

    inline void setUniform4(const char* name, float* a) const
    { glUniform4fv(glGetUniformLocation(mProgramID, name), 1, a); }

    inline void setUniform(const char* name, float a, float b) const
    { glUniform2f(glGetUniformLocation(mProgramID, name), a, b); }

    inline void setUniform(const char* name, float a, float b, float c) const
    { glUniform3f(glGetUniformLocation(mProgramID, name), a, b, c); }

    inline void setUniform(const char* name, float a, float b, float c, float d) const
    { glUniform4f(glGetUniformLocation(mProgramID, name), a, b, c, d); }

protected:

    bool mIsValid;
    typedef std::map<std::string,std::string> DefineMap;
    DefineMap mDefines;
    static void printInfoLog(GLuint objectID);
    GLuint mProgramID;
};

} // namepsace GlSplat

#endif
