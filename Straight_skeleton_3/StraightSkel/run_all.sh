#!/bin/sh

cd build

POLYHEDRONS="
$(ls ../res/polyhedrons/events/*.obj)
$(ls ../res/polyhedrons/*.obj)
$(ls ../res/polyhedrons/wolfgang/*.obj)"

POLYHEDRONS="
../res/polyhedrons/cubet.obj
../res/polyhedrons/held_convex.obj
../res/polyhedrons/culver_iron_maiden.obj
../res/polyhedrons/saddle_2.obj
../res/polyhedrons/seastar.obj
../res/polyhedrons/wedge_tabletop.obj
../res/polyhedrons/iron_maiden.obj
../res/polyhedrons/Schoenhardt.obj
../res/polyhedrons/verworrtakelt_fixed.obj
../res/polyhedrons/verworrtakelt.obj
../res/polyhedrons/wolfgang/Armadillo_000096.obj
../res/polyhedrons/wolfgang/Armadillo_000194.obj
../res/polyhedrons/wolfgang/aoi_asteroid.obj
../res/polyhedrons/wolfgang/bunny_small.obj
../res/polyhedrons/wolfgang/chess_bauer.obj
../res/polyhedrons/wolfgang/chinese_lion_174.obj
../res/polyhedrons/wolfgang/convex_piece_2.obj
../res/polyhedrons/wolfgang/hand_small.obj
../res/polyhedrons/wolfgang/small_sphere_shake.obj
../res/polyhedrons/wolfgang/venus_115.obj
../res/polyhedrons/wolfgang/venus_267.obj"

CONFIGS="
StraightSkel.ini
StraightSkel_min_vol.ini"

for CONFIG in $CONFIGS; do
  if [ ! -r "$CONFIG" ]; then
    echo "Error: Config $CONFIG not found."
    exit 1
  fi
  for POLYHEDRON in $POLYHEDRONS; do
    ./StraightSkel 3d load $POLYHEDRON --no-window --save --config $CONFIG
    if [ $? -ne 0 ]; then
      echo "$POLYHEDRON"
    fi
  done
done
