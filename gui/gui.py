from ursina import *
from ursina_config import *
from shader import twist_shader
from itertools import zip_longest
from time import perf_counter
import threading

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from src.ProtoDNASequencer import *

app = Ursina(title = TITLE, size=SIZE)

apply_window_config()
apply_camera_config()

color_palette = {
    "A": color.hex('#7ffbaa'),
    "T": color.hex('#d87ef7'),
    "C": color.hex('#f67da2'),
    "G": color.hex('#7f8afa'),
    "N": color.gray,
    " ": color.red,
    "backbone": color.hex('#8790aa')
    }

images = [
    ('../images/c1.png', lambda: print('clicked 1')),
    ('../images/c2.png', lambda: print('clicked 2')),
    ('../images/c3.png', lambda: print('clicked 3')),
    ('../images/c4.png', lambda: print('clicked 4')),
    ('../images/c5.png', lambda: print('clicked 5')),
    ('../images/c6.png', lambda: print('clicked 6')),
    ('../images/c7.png', lambda: print('clicked 7')),
    ('../images/c8.png', lambda: print('clicked 8')),
    ('../images/c9.png', lambda: print('clicked 9')),
    ('../images/c10.png', lambda: print('clicked 10')),
    ('../images/c11.png', lambda: print('clicked 11')),
    ('../images/c12.png', lambda: print('clicked 12')),
    ('../images/c13.png', lambda: print('clicked 13')),
    ('../images/c14.png', lambda: print('clicked 14')),
    ('../images/c15.png', lambda: print('clicked 15')),
    ('../images/c16.png', lambda: print('clicked 16')),
    ('../images/c17.png', lambda: print('clicked 17')),
    ('../images/c18.png', lambda: print('clicked 18')),
    ('../images/c19.png', lambda: print('clicked 19')),
    ('../images/c20.png', lambda: print('clicked 20')),
    ('../images/c21.png', lambda: print('clicked 21')),
    ('../images/c22.png', lambda: print('clicked 22')),
    ('../images/c23.png', lambda: print('clicked 23')),
    ('../images/c24.png', lambda: print('clicked 24')),
]

p = None
chromosomes = []
data_ready = False
pending_selection = None

chromosome_menu = Entity(parent=camera.ui)

menu_cols = 6
menu_spacing_x = 0.09
menu_spacing_y = 0.12
menu_rows = math.ceil(len(images) / menu_cols)
grid_width = (menu_cols - 1) * menu_spacing_x
grid_height = (menu_rows- 1) * menu_spacing_y

menu_title = Text(
    parent=chromosome_menu,
    text='Select a Chromosome',
    origin=(0, 0),
    y=grid_height / 2 + 0.15,
    scale=1.5,
)

def select_chromosome(index):
    global pending_selection
    chromosome_menu.enabled = False
    show_side_panels()
    back_btn.enabled = True
    txt.enabled = True
    txt.text = 'Loading sequence data...' if not data_ready else 'Rendering Gene.....'

    if data_ready:
        invoke(load_data, index, delay=0.1)
    else:
        pending_selection = index  # picked up in update() once parsing finishes

for i, (image, chrom_num) in enumerate(images):
    row = i // menu_cols
    col = i % menu_cols

    x = -grid_width / 2 + col * menu_spacing_x
    y = grid_height / 2 - row * menu_spacing_y

    button = Button(
        parent=chromosome_menu,
        texture=image,
        scale=0.07,
        position=(x, y),
        color=color.white,
        highlight_color=color.white.tint(-.1),
    )
    button.on_click = Func(select_chromosome, i)

left_panel = Entity(
    parent=camera.ui,
    origin=(0, 0),
    x=window.left.x - 5
)

left_background = Entity(
    parent=left_panel,
    model='quad',
    color=color.hex("#4F525C"),
    scale=(0.35, 1),
    origin=(-0.5, 0)
)

right_panel = Entity(
    parent=camera.ui,
    origin=(0, 0),
    x=window.right.x + 5
)

right_background = Entity(
    parent=right_panel,
    model='quad',
    color=color.hex("#4F525C"),
    scale=(0.35, 1),
    origin=(0.5, 0)
)

panel_width = right_background.scale_x

info_title = Text(
    parent=right_panel,
    text='fun facts or som',
    origin=(0, 0),
    x=-panel_width / 2,
    y=0.4,
    scale=1.2
)

info_text = Text(
    parent=right_panel,
    text='Select a chromosome to\nvisualize the DNA!',
    origin=(0, 0),
    x=-panel_width / 2,
    y=0.2,
    scale=0.8
)

txt = Text(text='Rendering Gene.....', orgin=(0, 0), scale=2)
target_camera_fov = camera.fov

def zoom_in():
    global target_camera_fov
    target_camera_fov = max(5, camera.fov - 5)

def zoom_out():
    global target_camera_fov
    target_camera_fov = min(20, camera.fov + 5)

zoom_in_btn = Button(
    text='+', 
    parent=left_panel, 
    scale=0.1,
    position=(.1, -.4),
    model='quad', 
    radius=0.5, 
    text_size=1, 
    text_color=color.white, 
    color=color.gray, 
    highlight_scale=1, 
    pressed_scale=1, 
    disabled=False)

zoom_out_btn = Button(
    text='-', 
    parent=left_panel, 
    scale=0.1,
    position=(.2, -.4),
    model='quad', 
    radius=0.1, 
    text_size=1, 
    text_color=color.white, 
    color=color.gray, 
    collider='box', 
    highlight_scale=1, 
    pressed_scale=1, 
    disabled=False)

zoom_in_btn.on_click = zoom_in
zoom_out_btn.on_click = zoom_out

back_btn = Button(
    text='< Back', parent=camera.ui, scale=(0.15, 0.06),
    position=(window.top_left.x + 0.1, window.top_left.y - 0.05),
    color=color.gray, text_color=color.white, enabled=False,
)

def go_back_to_menu():
    global helix, offset, scroll_offset, show_panels
    if helix is not None:
        destroy(helix)
        helix = None
    hide_side_panels()
    back_btn.enabled = False
    chromosome_menu.enabled = True
    offset = 0
    scroll_offset = 0.0

back_btn.on_click = go_back_to_menu

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
scroll_offset = 0.0
scroll_decay_speed = 12.0

show_panels = False

def refersh_visible_helix(direction=0):
    global helix
    if helix is not None:
        destroy(helix)

    visible_sequence = sequence[offset:offset + window_size]
    visible_complement = complement[offset:offset + window_size]
    helix = create_dna_helix(visible_sequence, visible_complement, (0, 0, 0))
    helix.set_shader_input('scroll_offset', scroll_offset)

def background_load():
    global p, chromosomes, data_ready
    p = ProtoDNASequencer("data/GCF_000001405.40_GRCh38.p14_genomic.fna")
    # p = ProtoDNASequencer("data/gene.fna")
    p.load_sequence()
    chromosomes = p.getChromosomeIds()
    data_ready = True

threading.Thread(target=background_load, daemon=True).start()

def load_data(index):
    global p, helix_length, sequence, complement

    # p = ProtoDNASequencer("data/GCF_000001405.40_GRCh38.p14_genomic.fna")
    # p.load_sequence()

    # chromosomes = p.getChromosomeIds()
    choice = chromosomes[index]

    p.sequenceChromosome(choice)

    sequence = p.getSeq(choice)
    complement = p.getComplement(choice)
    helix_length = len(sequence)

    txt.enabled = False
    refersh_visible_helix()

twisting = True
twist_amount = 1 # a number from 0 - 1 of a helix's twist that is applied

PANEL_SLIDE_DURATION = 0.3

def show_side_panels():
    left_panel.animate_x(window.left.x, duration=PANEL_SLIDE_DURATION, curve=curve.out_quad)
    right_panel.animate_x(window.right.x, duration=PANEL_SLIDE_DURATION, curve=curve.out_quad)

def hide_side_panels():
    left_panel.animate_x(window.left.x - 5, duration=PANEL_SLIDE_DURATION, curve=curve.in_quad)
    right_panel.animate_x(window.right.x + 5, duration=PANEL_SLIDE_DURATION, curve=curve.in_quad)

def update():
    global twist_amount, twisting, scroll_offset, pending_selection, show_panels

    

    if pending_selection is not None and data_ready:
        load_data(pending_selection)
        pending_selection = None

    if helix is None:
        return
    
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

    scroll_offset += (0 - scroll_offset) * min(scroll_decay_speed * time.dt, 1)
    helix.set_shader_input('scroll_offset', scroll_offset)
    helix.set_shader_input('twist_amount', twist_amount * helix.max_twist)

    twisting = camera.fov > 13

def input(key):
    global offset, scroll_offset
    if helix is None:
        return

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