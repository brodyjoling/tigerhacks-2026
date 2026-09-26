from ursina import window, camera, Vec2

TITLE = 'DNA Visualizer 4135'
SIZE = Vec2(1500, 1200)

def apply_window_config():
    window.title = 'DNA Visualizer 4135'
    window.borderless = False
    window.fullscreen = False
    window.exit_button.visible = False
    window.entity_counter.enabled = False
    window.cog_button.enabled = False
    window.collider_counter.enabled = False
    window.fps_counter.enabled = False

def apply_camera_config():
    camera.orthographic = True
    camera.fov = 20
    camera.position=(0, 0, -30)