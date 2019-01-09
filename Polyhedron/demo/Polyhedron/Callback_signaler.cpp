#include "Callback_signaler.h"

Callback_signaler::Callback_signaler()
  : is_canceled(false)
{
}

void Callback_signaler::emit_signal (int value) const
{
  Q_EMIT progressChanged(value);
}
