
#include"Qt_widget_delete_list_polygon.h"

void Qt_widget_movepolygon_helper::stateChanged(int i)
{
  if(i==2)
    activate();
  else
    if(i == 0)
      deactivate();
}

#include"Qt_widget_delete_list_polygon.moc"
