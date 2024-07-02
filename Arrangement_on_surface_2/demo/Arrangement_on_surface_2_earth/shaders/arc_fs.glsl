
#version 330

uniform vec4 u_color;
uniform vec4 u_plane;

in vec3 v_pos;
out vec4 out_color;

void main() {
  if (dot(u_plane, vec4(v_pos, 1)) < 0) discard;

  out_color = u_color;
}
