from ursina import *

app = Ursina()

camera.orthographic = True
camera.fov = 20
camera.position=(0, 0, -30)

window.title = 'DNA Visualizer'
window.borderless = False
window.fullscreen = False
window.exit_button.visible = False

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
    left_pos  = [position[0] - spacing, position[1], position[2]]
    right_pos = [position[0] + spacing, position[1], position[2]]

    backbone = Entity(model=None)

    left = Entity(parent=backbone, model='cube', scale=(0.15, spacing * length, 0.15), position=left_pos, color='#8790aa')
    right = Entity(parent=backbone, model='cube', scale=(0.15, spacing * length, 0.15), position=right_pos, color='#8790aa')

    backbone.combine()
    return backbone

def create_base_pairs(arr, position, spacing):
    length = len(arr)
    total_height = spacing * length

    base_pairs = Entity(model=None)
    pairs = []

    for i, pair in enumerate(arr):
        t = (i + 0.5) / length
        y = position[1] - total_height / 2 + t * total_height

        left_pos  = [position[0] - spacing /2, y, position[2]]
        right_pos = [position[0] + spacing /2, y, position[2]]

        pairs.append(Entity(parent=base_pairs, model='cube', scale=(spacing, 0.1, 0.1), position=left_pos, color=color.blue))
        pairs.append(Entity(parent=base_pairs, model='cube', scale=(spacing, 0.1, 0.1), position=right_pos, color=color.green))

    base_pairs.combine()
    return base_pairs

def create_dna_helix(arr, position):
    length = len(arr)

    dna = Entity(model=None)

    backbone = create_backbone(length, position, 0.5)
    base_pairs = create_base_pairs(arr, position, 0.5)

    backbone.parent = dna
    base_pairs.parent = dna
    return dna

position = [0, 0, 0]
helix = create_dna_helix(arr, position)

position = [-4, .5, 0]
helix = create_dna_helix(arr2, position)

def update():
    if mouse.left:
        camera.x -= mouse.velocity[0] * camera.fov
        camera.y -= mouse.velocity[1] * (camera.fov + 35)

    camera.x = clamp(camera.x, -30, 30)
    camera.y = clamp(camera.y, -30, 30)


def input(key):
    if key == 'scroll up':
        camera.fov = max(5, camera.fov - 2)
    if key == 'scroll down':
        camera.fov = min(60, camera.fov + 2)

app.run()