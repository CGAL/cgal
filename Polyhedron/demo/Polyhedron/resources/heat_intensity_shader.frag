#version 150
in vec4 color;
in vec4 fP;
in vec3 fN;
in float dist[6];
in float out_distance;
uniform highp vec4 light_pos;
uniform highp vec4 light_diff;
uniform highp vec4 light_spec;
uniform highp vec4 light_amb;
uniform highp float spec_power;
uniform int is_two_side;
uniform bool is_selected;
uniform bool is_clipbox_on;
uniform highp float near;
uniform highp float far;
uniform highp float width;
uniform highp float height;
uniform bool comparing;
uniform bool writing;
uniform sampler2D sampler;
uniform highp float alpha;
out vec4 out_color;

// Define this and resources/heat_intensity_shader.f will draw black instead of white lines.
// #define HEAT_METHOD_BLACK_LINES


float depth(float z)
{
  return (2 * near) / (far + near - z * (far - near));
}

void main(void) {

  if(is_clipbox_on)
    if(dist[0]>0.0 ||
      dist[1]>0.0 ||
      dist[2]>0.0 ||
      dist[3]>0.0 ||
      dist[4]>0.0 ||
      dist[5]>0.0)
        discard;

  float d = depth(gl_FragCoord.z);
  float test = texture(sampler, vec2(gl_FragCoord.x/width, gl_FragCoord.y/height)).r;
  if(comparing && d <= test)
    discard;
  vec4 temp_color;
  if(writing)
    temp_color = vec4(d,d,d,1.0);
  else
  {
    vec4 my_color = vec4(color.xyz, 1.0);
    vec3 L = light_pos.xyz - fP.xyz;
    vec3 V = -fP.xyz;
    vec3 N;
    if(fN == vec3(0.0,0.0,0.0))
      N = vec3(0.0,0.0,0.0);
    else
      N = normalize(fN);
    L = normalize(L);
    V = normalize(V);
    vec3 R = reflect(-L, N);
    vec4 diffuse;
    if(is_two_side == 1)
      diffuse = abs(dot(N,L)) * light_diff * color;
    else
      diffuse = max(dot(N,L), 0.0) * light_diff * my_color;
    vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec;
    vec4 ret_color = vec4((my_color*light_amb).xyz + diffuse.xyz + specular.xyz,1);
    if(is_selected)
      temp_color = vec4(ret_color.r+70.0/255.0, ret_color.g+70.0/255.0, ret_color.b+70.0/255.0, alpha);
    else
      temp_color = vec4(ret_color.xyz, alpha);
      
    vec3 c = temp_color.xyz;
    float h = out_distance;
    h = h * 20.;
    h = h - floor(h);
    h = (1./(1.+exp(-100.*(h-.55)))) + (1./(1.+exp(-100.*(-h+.45))));
    h = 1.-h;
#ifdef HEAT_METHOD_BLACK_LINES
      c = h*vec3(0.,0.,0.) + (1.-h)*c;
#else
      c = h*vec3(1.,1.,1.) + (1.-h)*c;
#endif
    out_color.rgb = c;
    out_color.a = alpha;
  }
}
