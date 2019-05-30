// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>

#include "GraphicsViewCurveInput.h"
#include "AlgebraicCurveExpressionParser.h"

#include <QGraphicsView>

namespace CGAL {
namespace Qt {



template < typename Coefficient_ >
void
GraphicsViewCurveInput<CGAL::Arr_algebraic_segment_traits_2<
                               Coefficient_> > ::
addNewAlgebraicCurve(const std::string& poly_expr_)
{
  this->poly_expr = poly_expr_;

  AlgebraicCurveExpressionParser parser(this->poly_expr);
  std::vector<struct term> terms;

  try {

    if (!parser.extract_poly_terms(terms))
    {
      throw std::invalid_argument("Invalid Expression");
    }

  } catch (std::invalid_argument) {

    QMessageBox msgBox;
    msgBox.setWindowTitle("Wrong Expression");
    msgBox.setIcon(QMessageBox::Critical);
    msgBox.setText(QString::fromStdString(poly_expr_ + " is invalid"));
    msgBox.setStandardButtons(QMessageBox::Ok);

    msgBox.exec();
    return;
  }

  //To create a curve
  Traits::Construct_curve_2 construct_curve
      = traits.construct_curve_2_object();

  Polynomial_2 polynomial;
  Polynomial_2 x = CGAL::shift(Polynomial_2(1),1,0);
  Polynomial_2 y = CGAL::shift(Polynomial_2(1),1,1);

  //extracting coefficients and power
  for (int i=0; i<terms.size(); i++)
  {
    polynomial += terms[i].coefficient 
                *CGAL::ipower(x,terms[i].x_exponent)
                *CGAL::ipower(y,terms[i].y_exponent);
  }

  //adding curve to the arrangement
  Curve_2 cv = construct_curve(polynomial);
  Q_EMIT generate( CGAL::make_object( cv ) );
}



} // namespace Qt
} // namespace CGAL
