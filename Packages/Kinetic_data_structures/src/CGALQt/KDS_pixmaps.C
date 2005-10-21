#include <CGAL/KDS/IO/internal/KDS_pixmaps.h>

namespace CGAL {
  namespace KDS {
    namespace internal {
      namespace pixmaps {
#include "play.xpm"
#include "pause.xpm"
#include "stop.xpm"
#include "play_to.xpm"
#include "play_through.xpm"
#include "reverse.xpm"
#include "faster.xpm"
#include "slower.xpm"
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
