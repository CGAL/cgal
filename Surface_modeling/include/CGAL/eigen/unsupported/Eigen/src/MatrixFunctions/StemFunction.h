// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2010 Jitse Niesen <jitse@maths.leeds.ac.uk>
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

#ifndef EIGEN_STEM_FUNCTION
#define EIGEN_STEM_FUNCTION

/** \ingroup MatrixFunctions_Module 
  * \brief Stem functions corresponding to standard mathematical functions.
  */
template <typename Scalar>
class StdStemFunctions
{
  public:

    /** \brief The exponential function (and its derivatives). */
    static Scalar exp(Scalar x, int)
    {
      return std::exp(x);
    }

    /** \brief Cosine (and its derivatives). */
    static Scalar cos(Scalar x, int n)
    {
      Scalar res;
      switch (n % 4) {
      case 0: 
	res = std::cos(x);
	break;
      case 1:
	res = -std::sin(x);
	break;
      case 2:
	res = -std::cos(x);
	break;
      case 3:
	res = std::sin(x);
	break;
      }
      return res;
    }

    /** \brief Sine (and its derivatives). */
    static Scalar sin(Scalar x, int n)
    {
      Scalar res;
      switch (n % 4) {
      case 0:
	res = std::sin(x);
	break;
      case 1:
	res = std::cos(x);
	break;
      case 2:
	res = -std::sin(x);
	break;
      case 3:
	res = -std::cos(x);
	break;
      }
      return res;
    }

    /** \brief Hyperbolic cosine (and its derivatives). */
    static Scalar cosh(Scalar x, int n)
    {
      Scalar res;
      switch (n % 2) {
      case 0:
	res = std::cosh(x);
	break;
      case 1:
	res = std::sinh(x);
	break;
      }
      return res;
    }
	
    /** \brief Hyperbolic sine (and its derivatives). */
    static Scalar sinh(Scalar x, int n)
    {
      Scalar res;
      switch (n % 2) {
      case 0:
	res = std::sinh(x);
	break;
      case 1:
	res = std::cosh(x);
	break;
      }
      return res;
    }

}; // end of class StdStemFunctions

#endif // EIGEN_STEM_FUNCTION
