// Copyright (c) 2016  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_THREE_EXCEPTIONS_H
#define CGAL_THREE_EXCEPTIONS_H

#include <CGAL/license/Three.h>


#include <exception>
#include <QString>
#include <QApplication>
#include <QScriptable>
#include <QScriptContext>
#include <QScriptEngine>
#include <boost/optional.hpp>
#include <QStringList>

namespace CGAL{
namespace Three{

class Script_exception : public std::runtime_error {
  QStringList bt;
public:
  Script_exception(QString what_arg,
                   QStringList backtrace)
    : std::runtime_error(what_arg.toStdString())
    , bt(backtrace)
  {}

  QStringList backtrace() const {
    return bt;
  }
};

template <typename T>
struct Optional_or_bool {
  typedef boost::optional<T> type;

  template <typename Callable>
  static type invoke(Callable f) { return type(f()); }
};

template <>
struct Optional_or_bool<void> {
  typedef bool type;

  template <typename Callable>
  static type invoke(Callable f) { f(); return true; }
};

enum Context { CURRENT_CONTEXT, PARENT_CONTEXT };

/// This function template wraps a `Callable` that can be called without
/// any argument (such as a lambda expression without arguments), and in
/// case the function call is in a Qt Script context, wraps the call in a
/// try/catch that converts the possible exception throw into a Qt Script
/// exception. That allows a Qt Script to catch the exception and deal
/// with it.
template <typename Callable>
typename Optional_or_bool<typename std::result_of<Callable()>::type>::type
wrap_a_call_to_cpp(Callable f,
                   QScriptable* qs = 0,
                   const char* file = 0,
                   int line = -1,
                   Context c = CURRENT_CONTEXT) {
  typedef typename std::result_of<Callable()>::type Callable_RT;
  typedef Optional_or_bool<Callable_RT> O_r_b;
  typedef typename O_r_b::type Return_type;

  const bool no_try_catch = qApp->property("no-try-catch").toBool();
  if(no_try_catch || qs == 0 || !qs->context()) return O_r_b::invoke(f);
  else
    try {
      return O_r_b::invoke(f);
    }
    catch(const std::exception& e) {
      const Script_exception* se = dynamic_cast<const Script_exception*>(&e);
      QScriptContext* context = qs->context();
      QStringList qt_bt = context->backtrace();
      if(se) qt_bt = se->backtrace();
      std::cerr << "Backtrace:\n";
      Q_FOREACH(QString s, qt_bt)
      {
        std::cerr << "  " << qPrintable(s) << std::endl;
      }
      context = context->parentContext();
      if(c == PARENT_CONTEXT) {
        std::cerr << "> parent";
        context = context->parentContext();
      } else {
        std::cerr << "> current";
      }
      std::cerr << " context: "
                << qPrintable(context->toString()) << std::endl;
      QString error;
      if(se) {
        error = se->what();
      } else {
        error = QObject::tr("Error");
        QString context;
        if(file != 0) context += QObject::tr(" at file %1").arg(file);
        if(line != -1) context += QString(":%1").arg(line);
        if(!context.isNull()) {
          error += context;
          qt_bt.push_front(QObject::tr("<cpp>") + context);
        }
        error += QString(": %1").arg(e.what());
      }
      QScriptValue v = context->throwError(error);
      v.setProperty("backtrace",
                    qScriptValueFromSequence(context->engine(), qt_bt));
      std::cerr << "result after throwError: "
                << qPrintable(v.toString()) << std::endl;
      return Return_type();
    }
}

} // end namespace Three
} // end namespace CGAL

#endif // CGAL_THREE_EXCEPTIONS_H
