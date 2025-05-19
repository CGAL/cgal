#include "MainWindow.h"

MainWindow::~MainWindow() {
  process->close();
  delete(process);
  delete(s);
  delete(ui);
}
