import sys
from pathlib import Path
from itertools import zip_longest
from time import perf_counter
import threading

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
sys.path.insert(0, str(PROJECT_ROOT))

from ursina import *
import parser
from ursina_config import *
from shader import twist_shader
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
    "backbone": color.hex('#8790aa'),
    "bg": color.hex('#456a7d'),
    "fg": color.hex("#588aa3"),
    "highlight": color.hex("#67a2bf"),
    "highlight2": color.hex("#72b2d2"),
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

class AppState:
    def __init__(self):
        self.p = None
        self.chromosomes = []
        self.data_ready = False
        self.pending_selection = None

        self.helix = None
        self.offset = 0
        self.window_size = 50
        self.helix_length = 0
        self.sequence = ''
        self.complement = ''

        self.scroll_offset = 0.0
        self.twisting = True
        self.twist_amount = 1
        self.target_camera_fov = camera.fov

        self.loading_dot_count = 0
        self.loading_dot_timer = 0

        self.last_search_text = ""

state = AppState()

def make_accent_panel(parent, scale, position=(0, 0), panel_color=None, radius=0.06, z=.02):
    return Entity(
        parent=parent,
        model=Quad(radius=radius, segments=8, scale=scale),
        color=panel_color or color_palette.get('bg'),
        position=position,
        z=z,
    )

def background_load(state):
    print("backgrond_load: starting")
    state.p = ProtoDNASequencer("data/GCF_000001405.40_GRCh38.p14_genomic.fna")
    # txt.text = "Indexing DNA sequence...\n(this could take a minute)"
    state.p.load_sequence()
    print("backgrond_load: getting chromosome IDs")
    
    txt.text = ""
    state.chromosomes = state.p.getChromosomeIds()
    state.data_ready = True
    print("backgrond_load: DONE!")

threading.Thread(target=background_load, args=(state,), daemon=True).start()

def refersh_visible_helix(state):
    if state.helix is not None:
        destroy(state.helix)

    visible_sequence = state.sequence[state.offset:state.offset + state.window_size]
    visible_complement = state.complement[state.offset:state.offset + state.window_size]
    state.helix = create_dna_helix(visible_sequence, visible_complement, (0, 0, 0))
    state.helix.set_shader_input('scroll_offset', state.scroll_offset)

    current_pair = state.offset + state.window_size // 2
    info_text.text = (
        f"You are at base pair {current_pair:,}\n"
        f"out of {state.helix_length:,} pairs\n"
        f"This pair is {state.sequence[current_pair].upper()} and {state.complement[current_pair].upper()}"
    )

def load_data(state, index):
    choice = state.chromosomes[index]
    state.p.sequenceChromosome(choice)
    state.sequence = state.p.getSeq(choice)
    state.complement = state.p.getComplement(choice)
    state.helix_length = len(state.sequence)
    state.offset = 20000

    txt.enabled = False
    refersh_visible_helix(state)

def show_side_panels(): 
    left_panel.animate_x(window.left.x, duration=PANEL_SLIDE_DURATION, curve=curve.out_quad)
    right_panel.animate_x(window.right.x, duration=PANEL_SLIDE_DURATION, curve=curve.out_quad)

def hide_side_panels():
    left_panel.animate_x(window.left.x - 5, duration=PANEL_SLIDE_DURATION, curve=curve.in_quad)
    right_panel.animate_x(window.right.x + 5, duration=PANEL_SLIDE_DURATION, curve=curve.in_quad)

def go_back_to_menu():
    if state.helix is not None:
        destroy(state.helix)
        state.helix = None
        hide_side_panels()
        chromosome_menu.enabled = True



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
    dna.max_twist = math.radians(36) * length

    elapsed_time = perf_counter() - start_time
    print(f"Time to render DNA helix: {elapsed_time:.6f} seconds")

    return dna

def zoom_in():
    state.target_camera_fov = max(5, camera.fov - 5)

def zoom_out():
    state.target_camera_fov = min(20, camera.fov + 5)

def handle_back_btn(state, index): # does this work?
    if not state.data_ready:
        return
    select_chromosome(state, index)

def select_chromosome(state, index):
    chromosome_menu.enabled = False
    show_side_panels()
    back_btn.enabled = True
    txt.enabled = True
    txt.text = 'Indexing Sequence Data...' if not state.data_ready else 'Rendering Gene.....'

    if state.data_ready:
        invoke(load_data, state, index, delay=0.5)
    else:
        state.pending_selection = index  # picked up in update() once parsing finishes

def show_side_panels():
    left_panel.animate_x(window.left.x, duration=PANEL_SLIDE_DURATION, curve=curve.out_quad)
    right_panel.animate_x(window.right.x, duration=PANEL_SLIDE_DURATION, curve=curve.out_quad)

def hide_side_panels():
    left_panel.animate_x(window.left.x - 5, duration=PANEL_SLIDE_DURATION, curve=curve.in_quad)
    right_panel.animate_x(window.right.x + 5, duration=PANEL_SLIDE_DURATION, curve=curve.in_quad)

def jump_to_line(state):
    line = int(search_field.text)
    state.offset = clamp(line - state.window_size // 2, 0, max(0, state.helix_length - state.window_size))
    refersh_visible_helix(state)

def jump_to_gene(state, gene=None):
    if gene == None:
        print("Searching for gene:", search_field.text)
        result = find_gene(search_field.text)
    else:
        print("Searching for gene:", gene)
        result = gene

    if result is not None:
        chromosome = result["chromosome"]
        if state.p.sequenceChromosome(chromosome) is None:
            print("Could not load chromosome:", chromosome)
            return

        start = result["start"] - 1
        end = result["end"]
        state.sequence = state.p.getSeq(chromosome)[start:end]
        state.complement = state.p.getComplement(chromosome)[start:end]
        state.helix_length = len(state.sequence)
        state.offset = 0
        state.scroll_offset = 0.0
        refersh_visible_helix(state)
    else:
        print("Gene not found")
SUGGESTED_GENES = ['BRCA1', 'TP53', 'EGFR', 'MYC', 'PTEN', 'APOE', 'CFTR']

suggestion_buttons = []

def clear_suggestions():
    global suggestion_buttons
    for b in suggestion_buttons:
        destroy(b)
    suggestion_buttons = []

def show_suggestions(state, matches):
    clear_suggestions()
    for i, gene in enumerate(matches[:5]):  # cap how many show at once
        btn = Button(
            parent=left_panel,
            text=gene,
            scale=(0.3, 0.04),
            position=(.5, -0.05 - i * 0.045),
            color=color_palette.get('fg'),
            highlight_color=color_palette.get('highlight2'),
            text_color=color.white,
        )
        btn.text_entity.font = '../fonts/IBMPlexSans-Regular.ttf'
        btn.on_click = Func(jump_to_gene, state, gene)
        suggestion_buttons.append(btn)


def update():
    # if search_field.text != state.last_search_text:
    #     state.last_search_text = search_field.text
    #     query = search_field.text.strip().upper()
    #     if query:
    #         matches = [g for g in SUGGESTED_GENES if g.startswith(query)]
    #         show_suggestions(state, matches)
    #     else:
    #         clear_suggestions()

    if not state.data_ready:
        state.loading_dot_timer += time.dt
        if state.loading_dot_timer >= 0.4:
            state.loading_dot_timer = 0
            state.loading_dot_count = (state.loading_dot_count + 1) % 4
            txt.text = "Indexing DNA sequence" + "." * state.loading_dot_count + "\n(this can take up to a minute)"

    search_placeholder.enabled = search_field.text == ''

    
    if state.pending_selection is not None and state.data_ready:
        load_data(state, state.pending_selection)
        state.pending_selection = None

    if state.helix is None:
        return

    if state.twisting:
        state.twist_amount = min(state.twist_amount + 2 * time.dt, 1)
    else:
        state.twist_amount = max(state.twist_amount - 2 * time.dt, 0)

    if camera.fov != state.target_camera_fov:
        if camera.fov < state.target_camera_fov:
            camera.fov += 20 * time.dt
        else:
            camera.fov -= 20 * time.dt
        if abs(camera.fov - state.target_camera_fov) < 1:
            camera.fov = state.target_camera_fov

    state.scroll_offset += (0 - state.scroll_offset) * min(scroll_decay_speed * time.dt, 1)
    state.helix.set_shader_input('scroll_offset', state.scroll_offset)
    state.helix.set_shader_input('twist_amount', state.twist_amount * state.helix.max_twist)
    state.twisting = camera.fov > 13

def input(key):
    if state.helix is None:
        return

    if key == 'scroll up' and state.offset > 0:
        step = min(SCROLL_STEP, state.offset)
        state.offset -= step
        state.scroll_offset -= step * SPACING
        refersh_visible_helix(state)
    if key == 'scroll down':
        step = min(SCROLL_STEP, state.helix_length - state.window_size - state.offset)
        state.offset += step
        state.scroll_offset += step * SPACING
        refersh_visible_helix(state)
    if key == 'enter':
        if search_field.text == "":
            print("Search field is empty")
            return
        if not search_field.text.isdigit():
            jump_to_gene(state)
        else:
            print("Jumping to line:", search_field.text)
            jump_to_line(state)
def get_complement_sequence(x):
    return ''.join([complement(char) for char in x])
def complement(x):
    if(x == 'a'):
        return 't'
    if(x == 't'):
        return 'a'
    if(x == 'c'):
        return 'g'
    if(x == 'g'):
        return 'c'
    if(x == 'n'):
        return 'n'
    return x
        

def find_gene(gene_name):
    # Placeholder implementation for finding a gene by name
    # This function should search the current chromosome data for the given gene name
    return parser.find_gene_by_name(gene_name)
chromosome_menu = Entity(parent=camera.ui)
menu_cols = 6
menu_spacing_x = 0.09
menu_spacing_y = 0.12
menu_rows = math.ceil(len(images) / menu_cols)
grid_width = (menu_cols - 1) * menu_spacing_x
grid_height = (menu_rows- 1) * menu_spacing_y

menu_backing = make_accent_panel(
    parent= chromosome_menu,
    scale=(grid_width+0.4, grid_height+0.5),
    position=(0, grid_height/2 - 0.02),
    panel_color=color_palette.get('fg'),
    radius=0.08,
)

menu_title = Text(
    parent=chromosome_menu,
    text='Select a Chromosome',
    origin=(0, 0),
    y=grid_height / 2 + 0.15,
    scale=1.5,
    font='../fonts/IBMPlexSans-Regular.ttf',
)

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
        text=f"{i+1}",
        font='../fonts/IBMPlexSans-Regular.ttf',
    )
    button.text_entity.font = '../fonts/IBMPlexSans-Regular.ttf'
    button.on_click = lambda i=i: handle_back_btn(state, i)

left_panel = Entity(
    parent=camera.ui,
    origin=(0, 0),
    x=window.left.x - 5
)

left_background = make_accent_panel(
    parent=left_panel,
    scale=(.7, 1),
    position=(0,0),
    panel_color=color_palette.get('fg'),
    radius=0.04,
    # z=-.001,
)

right_panel = Entity(
    parent=camera.ui,
    origin=(0, 0),
    x=window.right.x + 5
)

right_background = make_accent_panel(
    parent=right_panel,
    scale=(.7, 1),
    position=(0,0),
    panel_color=color_palette.get('fg'),
    radius=0.04,
)

panel_width = right_background.scale_x

info_title = Text(
    parent=right_panel,
    text='Information:',
    origin=(0, 0),
    x=-.172,
    y=0.4,
    scale=1.2,
    font='../fonts/IBMPlexSans-Regular.ttf',
)

info_text = Text(
    parent=right_panel,
    text='',
    origin=(0, 0),
    x=-.172,
    y=0.2,
    scale=0.8,
    font='../fonts/IBMPlexSans-Regular.ttf',
)
info_text.font = '../fonts/IBMPlexSans-Regular.ttf'

txt = Text(
    parent=camera.ui, 
    text='Loading...', 
    color=color.light_gray,
    origin=(0, 0),
    y=grid_height / 2 + 0.25,
    scale=1,
    font='../fonts/IBMPlexSans-Regular.ttf',
)
target_camera_fov = camera.fov

zoom_in_btn = Button(
    text='+', parent=left_panel, scale=(0.1, 0.06),
    position=(.1, -.425), highlight_color=color_palette.get("highlight2"),
    text_color=color.white, color=color_palette.get("highlight"), disabled=False,
    font='../fonts/IBMPlexSans-Regular.ttf',
)

zoom_out_btn = Button(
    text='-', parent=left_panel, scale=(0.1, 0.06),
    position=(.25, -.425), highlight_color=color_palette.get("highlight2"),
    text_color=color.white, color=color_palette.get("highlight"), disabled=False,
    font='../fonts/IBMPlexSans-Regular.ttf',
)

zoom_in_btn.on_click = zoom_in
zoom_out_btn.on_click = zoom_out

back_btn = Button(
    text='Back', parent=left_panel, scale=(0.15, 0.06),
    position=(.172, .45), highlight_color=color_palette.get("highlight2"),
    color=color_palette.get("highlight"), text_color=color.white, enabled=False,
    font='../fonts/IBMPlexSans-Regular.ttf',
)

search_field = InputField(
    parent=left_panel,
    default_value='',
    scale=(0.3, 0.05),
    position=(.175, 0.3),
    color=color_palette.get('highlight'),
)
search_field.text_field.font = '../fonts/IBMPlexSans-Regular.ttf'
search_field.text_field.highlight_color = color_palette.get('bg')
search_field.highlight_color = color_palette.get('highlight2')

search_placeholder = Text(
    parent=search_field,
    text='Enter Gene or Line',
    origin=(-0.5, 0),
    color=color.rgba(255, 255, 255, 120),
    scale=5,
    z=-0.1,
    font='../fonts/IBMPlexSans-Regular.ttf',
)


zoom_in_btn.text_entity.font = '../fonts/IBMPlexSans-Regular.ttf'
zoom_out_btn.text_entity.font = '../fonts/IBMPlexSans-Regular.ttf'
back_btn.text_entity.font = '../fonts/IBMPlexSans-Regular.ttf'
back_btn.on_click = go_back_to_menu

SPACING = 0.5
SCROLL_STEP = 1
scroll_decay_speed = 12.0

PANEL_SLIDE_DURATION = 0.3

app.run()