#include <CGAL/KDS/IO/internal/KDS_pixmaps.h>

namespace CGAL {
  namespace KDS {
    namespace internal {
      namespace pixmaps {
#include "KDS_play.xpm"
#include "KDS_pause.xpm"
#include "KDS_stop.xpm"
#include "KDS_play_to.xpm"
#include "KDS_play_through.xpm"
#include "KDS_reverse.xpm"
#include "KDS_faster.xpm"
#include "KDS_slower.xpm"
      }
      char * const * const play_xpm= pixmaps::play_xpm;
      char * const * const faster_xpm= pixmaps::faster_xpm;
      char * const * const play_through_xpm= pixmaps::play_through_xpm;
      char * const * const slower_xpm= pixmaps::slower_xpm;
      char * const * const pause_xpm= pixmaps::pause_xpm;
      char * const * const play_to_xpm= pixmaps::play_to_xpm;
      char * const * const reverse_xpm= pixmaps::reverse_xpm;
      char * const * const stop_xpm= pixmaps::stop_xpm;

    }
  }
}
