from ursina import *
from ursina_config import *
from shader import twist_shader
from itertools import zip_longest
from time import perf_counter

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.ProtoDNASequencer import *

app = Ursina(title = TITLE, size=SIZE)

apply_window_config()
apply_camera_config()

txt = Text(text='Rendering Gene.....', orgin=(0, 0), scale=2)

zoom_in_btn = Button(
    text='+', 
    parent=camera.ui, 
    scale=0.1,
    position=(-0.52, -0.4),
    model='quad', 
    radius=0.5, 
    origin=(0, 0), 
    text_origin=(0,0), 
    text_size=1, 
    text_color=color.white, 
    color=color.gray, 
    collider='box', 
    highlight_scale=1, 
    pressed_scale=1, 
    disabled=False)

zoom_out_btn = Button(
    text='-', 
    parent=camera.ui, 
    scale=0.1,
    position=(-0.3, -0.4),
    model='quad', 
    radius=0.1, 
    origin=(1, 0), 
    text_origin=(0,0), 
    text_size=1, 
    text_color=color.white, 
    color=color.gray, 
    collider='box', 
    highlight_scale=1, 
    pressed_scale=1, 
    disabled=False)

color_palette = {
    "A": color.hex('#7ffbaa'),
    "T": color.hex('#d87ef7'),
    "C": color.hex('#f67da2'),
    "G": color.hex('#7f8afa'),
    "N": color.gray,
    " ": color.red,
    "backbone": color.hex('#8790aa')
    }

target_camera_fov = camera.fov

def zoom_in():
    global target_camera_fov
    target_camera_fov = max(5, camera.fov - 5)

def zoom_out():
    global target_camera_fov
    target_camera_fov = min(20, camera.fov + 5)

zoom_in_btn.on_click = zoom_in
zoom_out_btn.on_click = zoom_out

def add_cube(vertices, triangles, colors, center, scale, cube_color): # claude code
    cx, cy, cz = center
    sx, sy, sz = scale
    hx, hy, hz = sx / 2, sy / 2, sz / 2
    base_index = len(vertices)

    corners = [
        (-hx, -hy, -hz), (hx, -hy, -hz), (hx, hy, -hz), (-hx, hy, -hz),
        (-hx, -hy, hz), (hx, -hy, hz), (hx, hy, hz), (-hx, hy, hz),
    ]
    for ox, oy, oz in corners:
        vertices.append((cx + ox, cy + oy, cz + oz))
        colors.append(cube_color)

    faces = [
        (0,1,2),(0,2,3), (4,6,5),(4,7,6),
        (0,4,5),(0,5,1), (3,2,6),(3,6,7),
        (0,3,7),(0,7,4), (1,5,6),(1,6,2),
    ]
    for a, b, c in faces:
        triangles.extend([base_index + a, base_index + b, base_index + c])

def create_dna_helix(sequence, complement, position):
    start_time = perf_counter()

    length = len(sequence)
    spacing = 0.5
    scale_y = spacing * length

    vertices, triangles, colors = [], [], []

    for i in range(length + 1):
        y = -scale_y / 2 + i * spacing
        add_cube(vertices, triangles, colors, (-spacing, y, 0), (0.15, spacing, 0.15), color_palette.get("backbone"))
        add_cube(vertices, triangles, colors, (spacing, y, 0), (0.15, spacing, 0.15), color_palette.get("backbone"))

    for i , (s, c) in enumerate(zip_longest(sequence, complement, fillvalue="")):
        offset_y = (i + 0.5) / length
        y = scale_y / 2 - offset_y * scale_y
        add_cube(vertices, triangles, colors, (-spacing / 2, y, 0), (spacing, 0.1, 0.1), color_palette.get(s.upper(), color.red))
        add_cube(vertices, triangles, colors, (spacing / 2, y, 0), (spacing, 0.1, 0.1), color_palette.get(c.upper(), color.red))

    dna = Entity(model=Mesh(vertices=vertices, triangles=triangles, colors=colors, mode='triangle'), position=position)

    dna.shader = twist_shader
    dna.set_shader_input('min_height', -scale_y / 2)
    dna.set_shader_input('span', scale_y)
    dna.set_shader_input('twist_amount', 1.0)
    dna.set_shader_input('scroll_offset', scroll_offset)
    dna.max_twist = math.radians(36) * length

    elapsed_time = perf_counter() - start_time
    print(f"Time to render DNA helix: {elapsed_time:.6f} seconds")

    return dna

helix = None
offset = 0
window_size = 50
p = None
helix_length = 0

sequence = ''
complement = ''

SPACING = 0.5
SCROLL_STEP = 5
scroll_shift_amount = SCROLL_STEP * SPACING
scroll_offset = 0.0
scroll_decay_speed = 12.0

def refersh_visible_helix(direction=0):
    global helix
    if helix is not None:
        destroy(helix)

    visible_sequence = sequence[offset:offset + window_size]
    visible_complement = complement[offset:offset + window_size]
    helix = create_dna_helix(visible_sequence, visible_complement, (0, 0, 0))
    helix.set_shader_input('scroll_offset', scroll_offset)

def load_data():
    global p, helix_length, sequence, complement

    p = ProtoDNASequencer("data/GCF_000001405.40_GRCh38.p14_genomic.fna")
    p.load_sequence()

    start_time = perf_counter()
    sequence = p.getSeq()
    elapsed_time = perf_counter() - start_time
    print(f"getSeq took {elapsed_time:.6f} seconds")

    start_time = perf_counter()
    complement = p.getComplement()
    elapsed_time = perf_counter() - start_time
    print(f"getComplement took {elapsed_time:.6f} seconds")

    helix_length = len(sequence)
    destroy(txt)
    refersh_visible_helix()

invoke(load_data, delay=1)

twisting = True
twist_amount = 1 # a number from 0 - 1 of a helix's twist that is applied

def update():
    global twist_amount, twisting, scroll_offset

    if mouse.left:
        camera.y -= mouse.velocity[1] * camera.fov

    camera.y = clamp(camera.y, -5, 5)

    if twisting:
        twist_amount = min(twist_amount + 2 * time.dt, 1)
    else:
        twist_amount = max(twist_amount - 2 * time.dt, 0)

    if camera.fov != target_camera_fov:
        if camera.fov < target_camera_fov:
            camera.fov = camera.fov + 20 * time.dt
        else:
            camera.fov = camera.fov - 20 * time.dt

        if abs(camera.fov - target_camera_fov) < 1:
            camera.fov = target_camera_fov
    if helix is not None:
        scroll_offset += (0 - scroll_offset) * min(scroll_decay_speed * time.dt, 1)
        helix.set_shader_input('scroll_offset', scroll_offset)
        helix.set_shader_input('twist_amount', twist_amount * helix.max_twist)

    if camera.fov <= 13:
        twisting = False
    else:
        twisting = True

def input(key):
    global offset, scroll_offset

    if key == 'scroll up' and offset > 0:
        step = min(SCROLL_STEP, offset)
        offset -= step
        scroll_offset -= step * SPACING
        refersh_visible_helix()
    if key == 'scroll down':
        step = min(SCROLL_STEP, helix_length - window_size - offset)
        offset += step
        scroll_offset += step * SPACING
        refersh_visible_helix()

app.run()