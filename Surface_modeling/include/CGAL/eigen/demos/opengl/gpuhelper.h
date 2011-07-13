// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008 Gael Guennebaud <gael.guennebaud@inria.fr>
//
// Eigen is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// Alternatively, you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of
// the License, or (at your option) any later version.
//
// Eigen is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License and a copy of the GNU General Public License along with
// Eigen. If not, see <http://www.gnu.org/licenses/>.

#ifndef EIGEN_GPUHELPER_H
#define EIGEN_GPUHELPER_H

#include <Eigen/Geometry>
#include <GL/gl.h>
#include <vector>

using namespace Eigen;

typedef Vector4f Color;

class GpuHelper
{
  public:

    GpuHelper();

    ~GpuHelper();

    enum ProjectionMode2D { PM_Normalized = 1, PM_Viewport = 2 };
    void pushProjectionMode2D(ProjectionMode2D pm);
    void popProjectionMode2D();

    /** Multiply the OpenGL matrix \a matrixTarget by the matrix \a mat.
        Essentially, this helper function automatically calls glMatrixMode(matrixTarget) if required
        and does a proper call to the right glMultMatrix*() function according to the scalar type
        and storage order.
        \warning glMatrixMode() must never be called directly. If your're unsure, use forceMatrixMode().
        \sa Matrix, loadMatrix(), forceMatrixMode()
    */
    template<typename Scalar, int _Flags>
    void multMatrix(const Matrix<Scalar,4,4, _Flags, 4,4>& mat, GLenum matrixTarget);

    /** Load the matrix \a mat to the OpenGL matrix \a matrixTarget.
        Essentially, this helper function automatically calls glMatrixMode(matrixTarget) if required
        and does a proper call to the right glLoadMatrix*() or glLoadIdentity() function according to the scalar type
        and storage order.
        \warning glMatrixMode() must never be called directly. If your're unsure, use forceMatrixMode().
        \sa Matrix, multMatrix(), forceMatrixMode()
    */
    template<typename Scalar, int _Flags>
    void loadMatrix(const Eigen::Matrix<Scalar,4,4, _Flags, 4,4>& mat, GLenum matrixTarget);

    template<typename Scalar, typename Derived>
    void loadMatrix(
        const Eigen::CwiseNullaryOp<Eigen::internal::scalar_identity_op<Scalar>,Derived>&,
        GLenum matrixTarget);

    /** Make the matrix \a matrixTarget the current OpenGL matrix target.
        Call this function before loadMatrix() or multMatrix() if you cannot guarantee that glMatrixMode()
        has never been called after the last loadMatrix() or multMatrix() calls.
        \todo provides a debug mode checking the sanity of the cached matrix mode.
    */
    inline void forceMatrixTarget(GLenum matrixTarget) {glMatrixMode(mCurrentMatrixTarget=matrixTarget);}

    inline void setMatrixTarget(GLenum matrixTarget);

    /** Push the OpenGL matrix \a matrixTarget and load \a mat.
    */
    template<typename Scalar, int _Flags>
    inline void pushMatrix(const Matrix<Scalar,4,4, _Flags, 4,4>& mat, GLenum matrixTarget);

    template<typename Scalar, typename Derived>
    void pushMatrix(
        const Eigen::CwiseNullaryOp<Eigen::internal::scalar_identity_op<Scalar>,Derived>&,
        GLenum matrixTarget);

    /** Push and clone the OpenGL matrix \a matrixTarget
    */
    inline void pushMatrix(GLenum matrixTarget);

    /** Pop the OpenGL matrix \a matrixTarget
    */
    inline void popMatrix(GLenum matrixTarget);

    void drawVector(const Vector3f& position, const Vector3f& vec, const Color& color, float aspect = 50.);
    void drawVectorBox(const Vector3f& position, const Vector3f& vec, const Color& color, float aspect = 50.);
    void drawUnitCube(void);
    void drawUnitSphere(int level=0);

    /// draw the \a nofElement first elements
    inline void draw(GLenum mode, uint nofElement);

    /// draw a range of elements
    inline void draw(GLenum mode, uint start, uint end);

    /// draw an indexed subset
    inline void draw(GLenum mode, const std::vector<uint>* pIndexes);

protected:

    void update(void);

    GLuint mColorBufferId;
    int mVpWidth, mVpHeight;
    GLenum mCurrentMatrixTarget;
    bool mInitialized;
};

/** Singleton shortcut
*/
extern GpuHelper gpu;


/** \internal
*/
template<bool RowMajor, int _Flags> struct GlMatrixHelper;

template<int _Flags> struct GlMatrixHelper<false,_Flags>
{
    static void loadMatrix(const Matrix<float, 4,4, _Flags, 4,4>&  mat) { glLoadMatrixf(mat.data()); }
    static void loadMatrix(const Matrix<double,4,4, _Flags, 4,4>& mat) { glLoadMatrixd(mat.data()); }
    static void multMatrix(const Matrix<float, 4,4, _Flags, 4,4>&  mat) { glMultMatrixf(mat.data()); }
    static void multMatrix(const Matrix<double,4,4, _Flags, 4,4>& mat) { glMultMatrixd(mat.data()); }
};

template<int _Flags> struct GlMatrixHelper<true,_Flags>
{
    static void loadMatrix(const Matrix<float, 4,4, _Flags, 4,4>&  mat) { glLoadMatrixf(mat.transpose().eval().data()); }
    static void loadMatrix(const Matrix<double,4,4, _Flags, 4,4>& mat) { glLoadMatrixd(mat.transpose().eval().data()); }
    static void multMatrix(const Matrix<float, 4,4, _Flags, 4,4>&  mat) { glMultMatrixf(mat.transpose().eval().data()); }
    static void multMatrix(const Matrix<double,4,4, _Flags, 4,4>& mat) { glMultMatrixd(mat.transpose().eval().data()); }
};

inline void GpuHelper::setMatrixTarget(GLenum matrixTarget)
{
    if (matrixTarget != mCurrentMatrixTarget)
        glMatrixMode(mCurrentMatrixTarget=matrixTarget);
}

template<typename Scalar, int _Flags>
void GpuHelper::multMatrix(const Matrix<Scalar,4,4, _Flags, 4,4>& mat, GLenum matrixTarget)
{
    setMatrixTarget(matrixTarget);
    GlMatrixHelper<_Flags&Eigen::RowMajorBit, _Flags>::multMatrix(mat);
}

template<typename Scalar, typename Derived>
void GpuHelper::loadMatrix(
    const Eigen::CwiseNullaryOp<Eigen::internal::scalar_identity_op<Scalar>,Derived>&,
    GLenum matrixTarget)
{
    setMatrixTarget(matrixTarget);
    glLoadIdentity();
}

template<typename Scalar, int _Flags>
void GpuHelper::loadMatrix(const Eigen::Matrix<Scalar,4,4, _Flags, 4,4>& mat, GLenum matrixTarget)
{
    setMatrixTarget(matrixTarget);
    GlMatrixHelper<(_Flags&Eigen::RowMajorBit)!=0, _Flags>::loadMatrix(mat);
}

inline void GpuHelper::pushMatrix(GLenum matrixTarget)
{
    setMatrixTarget(matrixTarget);
    glPushMatrix();
}

template<typename Scalar, int _Flags>
inline void GpuHelper::pushMatrix(const Matrix<Scalar,4,4, _Flags, 4,4>& mat, GLenum matrixTarget)
{
    pushMatrix(matrixTarget);
    GlMatrixHelper<_Flags&Eigen::RowMajorBit,_Flags>::loadMatrix(mat);
}

template<typename Scalar, typename Derived>
void GpuHelper::pushMatrix(
    const Eigen::CwiseNullaryOp<Eigen::internal::scalar_identity_op<Scalar>,Derived>&,
    GLenum matrixTarget)
{
    pushMatrix(matrixTarget);
    glLoadIdentity();
}

inline void GpuHelper::popMatrix(GLenum matrixTarget)
{
    setMatrixTarget(matrixTarget);
    glPopMatrix();
}

inline void GpuHelper::draw(GLenum mode, uint nofElement)
{
    glDrawArrays(mode, 0, nofElement);
}


inline void GpuHelper::draw(GLenum mode, const std::vector<uint>* pIndexes)
{
    glDrawElements(mode, pIndexes->size(), GL_UNSIGNED_INT, &(pIndexes->front()));
}

inline void GpuHelper::draw(GLenum mode, uint start, uint end)
{
    glDrawArrays(mode, start, end-start);
}

#endif // EIGEN_GPUHELPER_H
