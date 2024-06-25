#pragma once

/*************OPENGL WINDOW PARAMS*************/

#ifndef CGAL_WINDOW_WIDTH_INIT
#define CGAL_WINDOW_WIDTH_INIT 500
#endif

#ifndef CGAL_WINDOW_HEIGHT_INIT
#define CGAL_WINDOW_HEIGHT_INIT 450
#endif

#ifndef CGAL_WINDOW_SAMPLES
#define CGAL_WINDOW_SAMPLES 1
#endif

/*********************************************/
#ifndef CGAL_CLIPPING_PLANE_RENDERING_TRANSPARENCY
#define CGAL_CLIPPING_PLANE_RENDERING_TRANSPARENCY 0.5f
#endif

#ifndef CGAL_SIZE_POINTS
#define CGAL_SIZE_POINTS 7.0f
#endif

#ifndef CGAL_SIZE_EDGES
#define CGAL_SIZE_EDGES 3.1f
#endif

#ifndef CGAL_SIZE_RAYS
#define CGAL_SIZE_RAYS 3.1f
#endif

#ifndef CGAL_SIZE_LINES
#define CGAL_SIZE_LINES 3.1f
#endif

#ifndef CGAL_FACES_MONO_COLOR
#define CGAL_FACES_MONO_COLOR \
    {                    \
        60, 60, 200      \
    }
#endif

#ifndef CGAL_VERTICES_MONO_COLOR
#define CGAL_VERTICES_MONO_COLOR \
    {                       \
        200, 60, 60         \
    }
#endif

#ifndef CGAL_EDGES_MONO_COLOR
#define CGAL_EDGES_MONO_COLOR \
    {                    \
        0, 0, 0          \
    }
#endif

#ifndef CGAL_RAYS_MONO_COLOR
#define CGAL_RAYS_MONO_COLOR \
    {                   \
        0, 0, 0         \
    }
#endif

#ifndef CGAL_LINES_MONO_COLOR
#define CGAL_LINES_MONO_COLOR \
    {                    \
        0, 0, 0          \
    }
#endif

#ifndef CGAL_LIGHT_POSITION
#define CGAL_LIGHT_POSITION         \
    {                          \
        0.0f, 0.0f, 0.0f, 0.0f \
    }
#endif

#ifndef CGAL_AMBIENT_COLOR
#define CGAL_AMBIENT_COLOR         \
    {                         \
        0.6f, 0.5f, 0.5f, 1.f \
    }
#endif

#ifndef CGAL_DIFFUSE_COLOR
#define CGAL_DIFFUSE_COLOR          \
    {                          \
        0.9f, 0.9f, 0.9f, 1.0f \
    }
#endif

#ifndef CGAL_SPECULAR_COLOR
#define CGAL_SPECULAR_COLOR         \
    {                          \
        0.0f, 0.0f, 0.0f, 1.0f \
    }
#endif

#ifndef CGAL_SHININESS
#define CGAL_SHININESS 0.5f
#endif

#ifndef CGAL_CAM_MOVE_SPEED
#define CGAL_CAM_MOVE_SPEED 5.f
#endif

#ifndef CGAL_CAM_ROT_SPEED
#define CGAL_CAM_ROT_SPEED 5.f
#endif

#ifndef CGAL_CLIPPING_PLANE_MOVE_SPEED
#define CGAL_CLIPPING_PLANE_MOVE_SPEED 0.04f
#endif

#ifndef CGAL_CLIPPING_PLANE_ROT_SPEED
#define CGAL_CLIPPING_PLANE_ROT_SPEED 4.f
#endif

#ifndef CGAL_SCENE_ROT_SPEED
#define CGAL_SCENE_ROT_SPEED 0.5f
#endif