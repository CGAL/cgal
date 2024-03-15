
#version 330

uniform vec4 u_plane;

in vec3 v_color;
out vec4 out_color;

void main() {out_color = vec4(v_color, 1); }
