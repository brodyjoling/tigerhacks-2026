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

# ** COLORS **
pair_colors = {
    "A": color.hex('#7ffbaa'),
    "T": color.hex('#d87ef7'),
    "C": color.hex('#f67da2'),
    "G": color.hex('#7f8afa'),
    " ": color.red
    }

def create_backbone(length, position, spacing, scale_y):
    backbone = Entity(model=None)

    for i in range(length + 1):
        position_y = -scale_y / 2 + i * spacing

        left_pos  = [-spacing, position_y, 0]
        right_pos = [spacing, position_y, 0]

        Entity(parent=backbone, model='cube', scale=(0.15, spacing, 0.15), position=left_pos, color='#8790aa')
        Entity(parent=backbone, model='cube', scale=(0.15, spacing, 0.15), position=right_pos, color='#8790aa')
    
    backbone.combine()
    return backbone

def create_base_pairs(sequence, complement, length, position, spacing, scale_y):
    base_pairs = Entity(model=None)

    # loop through one sequence and complement at a time
    for i, (s, c) in enumerate(zip_longest(sequence, complement, fillvalue="")):
        offset_y = (i + 0.5) / length
        position_y = scale_y / 2 - offset_y * scale_y

        left_pos  = [-spacing / 2, position_y, 0]
        right_pos = [spacing / 2, position_y, 0]
        
        Entity(parent=base_pairs, model='cube', scale=(spacing, 0.1, 0.1), position=left_pos, color=pair_colors.get(s))
        Entity(parent=base_pairs, model='cube', scale=(spacing, 0.1, 0.1), position=right_pos, color=pair_colors.get(c))

    base_pairs.combine()
    return base_pairs

def create_dna_helix(sequence, complement, position):
    start_time = perf_counter()

    length = len(sequence)
    spacing = 0.5

    dna = Entity(model=None, position=position)

    scale_y = spacing * length
    backbone = create_backbone(length, position, spacing, scale_y)
    base_pairs = create_base_pairs(sequence, complement, length, position, spacing, scale_y)

    backbone.parent = dna
    base_pairs.parent = dna

    dna.combine()

    ys = [v[1] for v in dna.model.vertices]
    min_height, max_height = min(ys), max(ys)

    dna.shader = twist_shader
    dna.set_shader_input('min_height', min_height)
    dna.set_shader_input('span', max_height - min_height)
    dna.set_shader_input('twist_amount', 1.0)

    dna.max_twist = math.radians(36) * length

    elapsed_time = perf_counter() - start_time
    print(f"Time to render DNA helix: {elapsed_time} seconds")

    return dna

helices = []
pending_data = []
spawn_x = 0
spacing = 5
helices_per_frame = 3

def queue_helix(sequence, complement):
    global spawn_x
    pending_data.append((sequence, complement, [spawn_x, 0, 0]))
    spawn_x += spacing

def load_data():
    p = ProtoDNASequencer("data/gene.fna")
    p.load_sequence()
    # p.printSequenceAndComplement()
    queue_helix(p.getSeq(), p.getComplement())
    destroy(txt)

invoke(load_data, delay=1)

twisting = True
twist_amount = 1 # a number from 0 - 1 of a helix's twist that is applied

def update():
    global twist_amount, twisting

    if mouse.left:
        camera.x -= mouse.velocity[0] * camera.fov
        camera.y -= mouse.velocity[1] * (camera.fov + 0)

    # camera.x = clamp(camera.x, -30, 30)
    # camera.y = clamp(camera.y, -25, 25)

    if twisting:
        twist_amount = min(twist_amount + 2 * time.dt, 1)
    else:
        twist_amount = max(twist_amount - 2 * time.dt, 0)

    for h in helices:
        h.set_shader_input('twist_amount', twist_amount * h.max_twist)

    if camera.fov < 15:
        twisting = False
    else:
        twisting = True

    # Build Helices
    for _ in range(helices_per_frame):
        if not pending_data:
            break
        sequence, complement, position = pending_data.pop(0)
        helices.append(create_dna_helix(sequence, complement, position))


def input(key):
    global twisting

    if key == 'scroll up':
        camera.fov = max(5, camera.fov - 2)
    if key == 'scroll down':
        camera.fov = min(200, camera.fov + 2)

app.run()