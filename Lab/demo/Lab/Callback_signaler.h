#ifndef CALLBACK_SIGNALER_H
#define CALLBACK_SIGNALER_H

#include <QObject>

#include <QtCore/qglobal.h>
#ifdef scene_callback_signaler_EXPORTS
#  define SCENE_CALLBACK_SIGNALER_EXPORT Q_DECL_EXPORT
#else
#  define SCENE_CALLBACK_SIGNALER_EXPORT Q_DECL_IMPORT
#endif

class SCENE_CALLBACK_SIGNALER_EXPORT Callback_signaler
  : public QObject
{
  Q_OBJECT

public:

  bool is_canceled;

  Callback_signaler();

  void emit_signal (int value) const;

public Q_SLOTS:

  void cancel()
  {
    is_canceled = true;
  }

Q_SIGNALS:

  void progressChanged(int) const;
};


#endif // CALLBACK_SIGNALER_H
