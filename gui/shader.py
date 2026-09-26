from ursina import Shader

twist_shader = Shader(
    language=Shader.GLSL,
    vertex='''
#version 150
uniform mat4 p3d_ModelViewProjectionMatrix;
in vec4 p3d_Vertex;
in vec4 p3d_Color;
out vec4 vertex_color;

uniform float twist_amount;  // radians, total twist over the whole height
uniform float min_height;
uniform float span;

void main() {
    vec4 vert = p3d_Vertex;
    float t = clamp((vert.y - min_height) / span, 0.0, 1.0);

    float twist_angle = t * twist_amount;
    float tx = vert.x * cos(twist_angle) - vert.z * sin(twist_angle);
    float tz = vert.x * sin(twist_angle) + vert.z * cos(twist_angle);
    vert.x = tx;
    vert.z = tz;

    vertex_color = p3d_Color;
    gl_Position = p3d_ModelViewProjectionMatrix * vert;
}
''',
    fragment='''
#version 150
in vec4 vertex_color;
out vec4 fragColor;
void main() { fragColor = vertex_color; }
''',
)