#version 150
in vec4 color;
uniform vec3 dirView;
uniform vec3 plane_normal;
uniform vec3 plane_pos;
uniform bool is_selected;
out vec4 out_color;

void main(void) {
vec4 t_color = color;
vec3 dir = vec3(plane_pos.x - dirView.x, plane_pos.y - dirView.y, plane_pos.z - dirView.z);
if(dot(dir, plane_normal)>0.0)
  t_color = vec4(1.0-color.r, 1.0-color.g, 1.0-color.b, color.a);
if(is_selected)
  out_color = vec4(t_color.r+70.0/255.0, t_color.g+70.0/255.0, t_color.b+70.0/255.0, 1.0);
else
  out_color = t_color;
}
