
varying highp vec4 color;
uniform highp vec3 dirView;
uniform highp vec3 plane_normal;
uniform highp vec3 plane_pos;
uniform bool is_selected;

void main(void) {
highp vec4 t_color = color;
highp vec3 dir = vec3(plane_pos.x - dirView.x, plane_pos.y - dirView.y, plane_pos.z - dirView.z);
if(dot(dir, plane_normal)>0.0)
  t_color = vec4(1.0-color.r, 1.0-color.g, 1.0-color.b, color.a);
if(is_selected)
  gl_FragColor = vec4(t_color.r+70.0/255.0, t_color.g+70.0/255.0, t_color.b+70.0/255.0, 1.0);
else
  gl_FragColor = t_color;
}
