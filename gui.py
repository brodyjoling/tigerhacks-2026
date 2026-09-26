from ursina import *
from ursina_config import *
from shader import twist_shader

app = Ursina(title = TITLE, size=SIZE)

apply_window_config()
apply_camera_config()

# ** TEST DICTIONARY ARRAYS **
arr = [
    {"base": "A", "complement": "T"},
    {"base": "C", "`complement": "G"},
    {"base": "A", "complement": "T"},
    {"base": "C", "`complement": "G"},
    {"base": "A", "complement": "T"},
    {"base": "C", "`complement": "G"},
    {"base": "A", "complement": "T"},
    {"base": "A", "complement": "T"},
    ]

arr2 = [
    {"base": "A", "complement": "T"},
    {"base": "C", "`complement": "G"},
    {"base": "A", "complement": "T"},
    {"base": "C", "`complement": "G"},
    ]

def create_backbone(length, position, spacing):
    scale_y = spacing * length
    backbone = Entity(model=None)

    for i in range(length + 1):
        position_y = -scale_y / 2 + i * spacing

        left_pos  = [-spacing, position_y, 0]
        right_pos = [spacing, position_y, 0]

        Entity(parent=backbone, model='cube', scale=(0.15, spacing, 0.15), position=left_pos, color='#8790aa')
        Entity(parent=backbone, model='cube', scale=(0.15, spacing, 0.15), position=right_pos, color='#8790aa')
    
    backbone.combine()
    return backbone

def create_base_pairs(arr, position, spacing):
    length = len(arr)
    total_height = spacing * length

    base_pairs = Entity(model=None)

    for i, pair in enumerate(arr):
        t = (i + 0.5) / length
        position_y = -total_height / 2 + t * total_height

        left_pos  = [-spacing / 2, position_y, 0]
        right_pos = [spacing / 2, position_y, 0]

        Entity(parent=base_pairs, model='cube', scale=(spacing, 0.1, 0.1), position=left_pos, color=color.blue)
        Entity(parent=base_pairs, model='cube', scale=(spacing, 0.1, 0.1), position=right_pos, color=color.green)

    base_pairs.combine()
    return base_pairs

def create_dna_helix(arr, position):
    length = len(arr)

    dna = Entity(model=None, position=position)

    backbone = create_backbone(length, position, 0.5)
    base_pairs = create_base_pairs(arr, position, 0.5)

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

    return dna

helices = [
    create_dna_helix(arr, [0, 0, 0]),
    create_dna_helix(arr2, [-4, .5, 0]),
]

twisting = True
twist_amount = 1 # a number from 0 - 1 of a helix's twist that is applied

def update():
    global twist_amount, twisting

    if mouse.left:
        camera.x -= mouse.velocity[0] * camera.fov
        camera.y -= mouse.velocity[1] * (camera.fov + 0)

    camera.x = clamp(camera.x, -30, 30)
    camera.y = clamp(camera.y, -25, 25)

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


def input(key):
    global twisting

    if key == 'scroll up':
        camera.fov = max(5, camera.fov - 2)
    if key == 'scroll down':
        camera.fov = min(60, camera.fov + 2)

app.run()