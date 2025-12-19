// Copyright (c) 2016  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_THREE_EXCEPTIONS_H
#define CGAL_THREE_EXCEPTIONS_H

#include <CGAL/license/Three.h>


#include <exception>
#include <QString>
#include <QApplication>
#include <optional>
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
  typedef std::optional<T> type;

  template <typename Callable>
  static type invoke(Callable f) { return type(f()); }
};

template <>
struct Optional_or_bool<void> {
  typedef bool type;

  template <typename Callable>
  static type invoke(Callable f) { f(); return true; }
};

/// This function template wraps a `Callable` that can be called without
/// any argument (such as a lambda expression without arguments), and in
/// case the function call is in a Qt Script context, wraps the call in a
/// try/catch that converts the possible exception throw into a Qt Script
/// exception. That allows a Qt Script to catch the exception and deal
/// with it.
template <typename Callable>
typename Optional_or_bool<typename cpp11::result_of<Callable()>::type>::type
wrap_a_call_to_cpp(Callable f,
                   QObject* object = 0,
                   const char* file = 0,
                   int line = -1) {
  typedef decltype(f()) Callable_RT;
  typedef Optional_or_bool<Callable_RT> O_r_b;
  QJSEngine* script_engine = qjsEngine(object);
  const bool no_try_catch = qApp->property("no-try-catch").toBool();
  if(no_try_catch || script_engine == 0) return O_r_b::invoke(f);
  else
    try {
      return O_r_b::invoke(f);
    }
    catch(const std::exception& e) {
      const Script_exception* se = dynamic_cast<const Script_exception*>(&e);
      QStringList qt_bt;
      if(se) {
        qt_bt = se->backtrace();
        if(qt_bt.size() != 0) std::cerr << "Backtrace:\n";
      }
      QJSValue js_bt = script_engine->newArray(qt_bt.size());
      if(qt_bt.size() != 0) {
        quint32 i = 0;
        for(auto s: qt_bt)
        {
          std::cerr << "  " << qPrintable(s) << std::endl;
          js_bt.setProperty(i++, s);
        }
      }
      QString error;
      if(se) {
        error = se->what();
      } else {
        error = QObject::tr("Error");
        QString context;
        if(file != 0) context += QObject::tr(" at file %1").arg(file);
        if(line != -1) context += QString(":%1").arg(line);
        error += QString(": %1").arg(e.what());
      }
      QJSValue error_value = script_engine->newErrorObject(QJSValue::GenericError, error);
      error_value.setProperty("backtrace", js_bt);
      script_engine->throwError(error_value);
      return {};
    }
}

} // end namespace Three
} // end namespace CGAL

#endif // CGAL_THREE_EXCEPTIONS_H
