import asyncio
import os
import threading
import time
from typing import Dict, List, Optional, Tuple

import glfw
import numpy as np
import toga
from OpenGL.GL import *
from OpenGL.GLU import *
import OpenGL.GLUT as ogl_glut
from OpenGL.GLUT import *
from toga.style import Pack
from toga.style.pack import CENTER, COLUMN, LEFT, ROW

Atom = Tuple[str, Tuple[float, float, float]]
Frame = List[Atom]


class MolecularViewer:
    """3D molecular viewer with trajectory support using GLFW + OpenGL."""

    def __init__(self):
        # GLFW-related attributes
        self.glfw_initialized = False
        self.window = None
        self._render_thread: Optional[threading.Thread] = None
        self._state_lock = threading.RLock()
        self._status_message = ""
        self._slider_update_in_progress = False
        self._gl_resources_ready = False
        self._sphere_lists: Dict[float, int] = {}
        self._glut_ready = False
        self._projection_update_requested = False
        self._scene_radius_cache: Dict[int, float] = {}
        self._fixed_scene_radius: Optional[float] = None

        # Sample molecule data (used until a trajectory is loaded)
        self.molecule_data: Frame = [
            ("C", (0.0, 0.0, 0.0)),
            ("H", (1.0, 0.0, 0.0)),
            ("H", (-1.0, 0.0, 0.0)),
            ("O", (0.0, 1.0, 0.0)),
            ("N", (0.0, -1.0, 0.0)),
        ]

        # Annotation / measurement options
        self.show_atom_numbers = False
        self.show_atom_symbols = False
        self.label_font = getattr(ogl_glut, "GLUT_BITMAP_HELVETICA_18")
        self.measurement_result = ""
        self.active_measurement_type: Optional[str] = None
        self.active_measurement_indices: List[int] = []

        # Atom colors and approximate covalent/visual radii (Å)
        self.atom_colors = {
            "H": (1.0, 1.0, 1.0),
            "C": (0.0, 1.0, 1.0),
            "N": (0.0, 0.0, 1.0),
            "O": (1.0, 0.0, 0.0),
            "F": (0.6, 0.9, 0.6),
            "P": (1.0, 0.5, 0.0),
            "S": (1.0, 1.0, 0.0),
            "Cl": (0.0, 1.0, 0.0),
            "Br": (0.6, 0.2, 0.2),
            "I": (0.5, 0.0, 0.7),
        }
        self.covalent_radii = {
            "H": 0.31,
            "B": 0.85,
            "C": 0.76,
            "N": 0.71,
            "O": 0.66,
            "F": 0.57,
            "P": 1.07,
            "S": 1.05,
            "Cl": 1.02,
            "Br": 1.20,
            "I": 1.39,
        }
        self.visual_radii = {
            "H": 0.25,
            "C": 0.42,
            "N": 0.40,
            "O": 0.40,
            "F": 0.38,
            "P": 0.50,
            "S": 0.50,
            "Cl": 0.47,
            "Br": 0.52,
            "I": 0.58,
        }
        self.vdw_radii = {
            "H": 1.20,
            "He": 1.40,
            "Li": 1.82,
            "Be": 1.53,
            "B": 1.92,
            "C": 1.70,
            "N": 1.55,
            "O": 1.52,
            "F": 1.47,
            "Ne": 1.54,
            "Na": 2.27,
            "Mg": 1.73,
            "Al": 1.84,
            "Si": 2.10,
            "P": 1.80,
            "S": 1.80,
            "Cl": 1.75,
            "Ar": 1.88,
            "K": 2.75,
            "Ca": 2.31,
            "Br": 1.85,
            "I": 1.98,
        }
        self.default_visual_radius = 0.45
        self.default_vdw_radius = 1.70
        self.bond_tolerance = 0.45
        self.measurement_highlight_color = (1.0, 1.0, 0.0)
        self.measurement_line_color = (1.0, 1.0, 1.0)

        self.rotation_angle_x = 0.0
        self.rotation_angle_y = 0.0
        self.molecule_rotation_x = 0.0
        self.molecule_rotation_y = 0.0
        self.molecule_rotation_z = 0.0
        self.continuous_rotation_axis: Optional[str] = None
        self.continuous_rotation_direction = 0
        self.continuous_rotation_speed = 60.0
        self._last_continuous_rotation_time: Optional[float] = None
        self.translation_vector = [0.0, 0.0, 0.0]
        self.base_camera_distance = 5.0
        self.mouse_down = False
        self.right_mouse_down = False
        self.last_mouse_pos = None
        self.connection_distance = 1.5
        self.atom_scale_factor = 1.0
        self.atom_display_style = "Line style"
        self.bond_thickness_scale_factor = 1.0
        self.base_bond_radius = 0.05
        self.bonds: List[Tuple[int, int]] = []
        self.display = [800, 600]

        # Visualization style and zoom
        self.visualization_style = "Orthographic"
        self.zoom_factor = 1.0
        self.min_zoom = 0.1
        self.max_zoom = 10.0

        # Trajectory data
        self.frames: List[Frame] = []
        self.current_frame = 0
        self.playing = False
        self.frame_elements: List[str] = []
        self.frame_offsets: List[int] = []

        # Cache for bond lists by frame and connectivity settings
        self._bond_cache: Dict[Tuple[int, float, float], List[Tuple[int, int]]] = {}
        self._bond_cache_limit = 250
        self.bond_rendering_mode = "Static first frame"
        self._reference_bonds: Optional[List[Tuple[int, int]]] = None

        # Playback timing
        self.frame_skip = 1
        self.update_delay = 0.05
        self.last_frame_time = time.time()
        self.fast_playback_mode = False
        self.loop_mode = "Loop"
        self.play_direction = 1
        self.auto_frame_zoom = False

        # Orthorhombic box properties
        self.box_sizes = (4.0, 4.0, 4.0)
        self.box_centering_mode = "Geometric center"
        self.show_box = False

        # UI-loop pointer (set by subclass)
        self.loop = None

    # ------------------------------------------------------------------
    # State helpers
    # ------------------------------------------------------------------
    def set_status_message(self, message: str):
        self._status_message = message
        if self.loop and not self.loop.is_closed() and getattr(self, "loading_label", None):
            self.loop.call_soon_threadsafe(self._update_loading_label, message)

    def _update_loading_label(self, message: str):
        if getattr(self, "loading_label", None):
            self.loading_label.text = message

    def get_current_frame_data(self) -> Frame:
        with self._state_lock:
            if self.frames and 0 <= self.current_frame < len(self.frames):
                return list(self.frames[self.current_frame])
            return list(self.molecule_data)

    def invalidate_bond_cache(self):
        with self._state_lock:
            self._bond_cache.clear()
            self._reference_bonds = None

    def invalidate_scene_cache(self):
        with self._state_lock:
            self._scene_radius_cache.clear()
            self._fixed_scene_radius = None

    def request_projection_update(self):
        with self._state_lock:
            self._projection_update_requested = True

    def consume_projection_update_request(self) -> bool:
        with self._state_lock:
            update_requested = self._projection_update_requested
            self._projection_update_requested = False
            return update_requested

    def _bond_cache_key(self, frame_index: int) -> Tuple[int, float, float]:
        return (frame_index, round(self.connection_distance, 4), round(self.bond_tolerance, 4))

    def _trim_bond_cache_if_needed(self):
        if len(self._bond_cache) > self._bond_cache_limit:
            # Remove the oldest inserted item in a simple, deterministic way.
            first_key = next(iter(self._bond_cache))
            del self._bond_cache[first_key]

    # ------------------------------------------------------------------
    # OpenGL / GLFW
    # ------------------------------------------------------------------
    def init_glfw(self):
        if self.glfw_initialized:
            return

        if not glfw.init():
            raise RuntimeError("GLFW initialization failed")

        self.window = glfw.create_window(
            self.display[0], self.display[1], "3D Molecular Viewer", None, None
        )
        if not self.window:
            glfw.terminate()
            raise RuntimeError("GLFW window creation failed")

        glfw.make_context_current(self.window)
        self._ensure_glut_ready()

        glViewport(0, 0, self.display[0], self.display[1])
        self.set_projection()

        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_NORMALIZE)
        glShadeModel(GL_SMOOTH)

        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, [0.2, 0.2, 0.2, 1.0])
        light_ambient = [0.0, 0.0, 0.0, 1.0]
        light_diffuse = [1.0, 1.0, 1.0, 1.0]
        light_specular = [1.0, 1.0, 1.0, 1.0]
        light_position = [1.0, 1.0, 1.0, 0.0]
        glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient)
        glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse)
        glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular)
        glLightfv(GL_LIGHT0, GL_POSITION, light_position)

        # Set GLFW callbacks
        glfw.set_window_size_callback(self.window, self.on_window_resize)
        glfw.set_mouse_button_callback(self.window, self.on_mouse_button)
        glfw.set_cursor_pos_callback(self.window, self.on_cursor_pos)
        glfw.set_scroll_callback(self.window, self.on_scroll)

        self._build_gl_resources()
        self.glfw_initialized = True

    def _build_gl_resources(self):
        if self._gl_resources_ready:
            return

        sphere_radii = (
            set(self.visual_radii.values())
            | set(self.vdw_radii.values())
            | {self.default_visual_radius, self.default_vdw_radius}
        )
        for radius in sorted(sphere_radii):
            list_id = glGenLists(1)
            glNewList(list_id, GL_COMPILE)
            quad = gluNewQuadric()
            try:
                gluQuadricNormals(quad, GLU_SMOOTH)
                gluSphere(quad, radius, 28, 28)
            finally:
                gluDeleteQuadric(quad)
            glEndList()
            self._sphere_lists[radius] = list_id

        self._gl_resources_ready = True

    def _destroy_gl_resources(self):
        if not self._gl_resources_ready:
            return

        for list_id in self._sphere_lists.values():
            try:
                glDeleteLists(list_id, 1)
            except Exception:
                pass
        self._sphere_lists.clear()
        self._gl_resources_ready = False

    def _ensure_glut_ready(self):
        if self._glut_ready:
            return
        try:
            glutInit()
            self._glut_ready = True
        except Exception:
            self._glut_ready = False

    def _draw_text_3d(self, position, text_value: str):
        if not text_value or not self._glut_ready:
            return

        glDisable(GL_LIGHTING)
        glColor3f(1.0, 1.0, 0.0)
        glRasterPos3f(*position)
        for ch in text_value:
            glutBitmapCharacter(self.label_font, ord(ch))
        glEnable(GL_LIGHTING)


    def _draw_text_2d(self, x: float, y: float, text_value: str, color=(1.0, 1.0, 1.0)):
        if not text_value or not self._glut_ready:
            return

        glColor3f(*color)
        glRasterPos2f(x, y)
        for ch in text_value:
            glutBitmapCharacter(self.label_font, ord(ch))

    def draw_cartesian_axes_overlay(self):
        """Draw a small Cartesian axis triad at the bottom-left corner.

        The triad uses the same scene rotations as the molecule, so it rotates
        together with the current orientation while remaining anchored on screen.
        """
        inset_size = max(110, min(self.display[0], self.display[1]) // 5)
        margin = 18
        vx = margin
        vy = margin
        vw = inset_size
        vh = inset_size
        axis_len = 0.9

        prev_viewport = glGetIntegerv(GL_VIEWPORT)

        glViewport(vx, vy, vw, vh)

        glMatrixMode(GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity()
        gluPerspective(35.0, 1.0, 0.1, 10.0)

        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glLoadIdentity()
        glTranslatef(0.0, 0.0, -3.2)
        glRotatef(self.rotation_angle_x, 1, 0, 0)
        glRotatef(self.rotation_angle_y, 0, 1, 0)

        glDisable(GL_LIGHTING)
        glClear(GL_DEPTH_BUFFER_BIT)
        glLineWidth(3.0)

        glBegin(GL_LINES)
        glColor3f(0.0, 0.0, 1.0)
        glVertex3f(0.0, 0.0, 0.0)
        glVertex3f(axis_len, 0.0, 0.0)

        glColor3f(0.0, 1.0, 0.0)
        glVertex3f(0.0, 0.0, 0.0)
        glVertex3f(0.0, axis_len, 0.0)

        glColor3f(1.0, 0.0, 0.0)
        glVertex3f(0.0, 0.0, 0.0)
        glVertex3f(0.0, 0.0, axis_len)
        glEnd()

        if self._glut_ready:
            glColor3f(0.0, 0.0, 1.0)
            glRasterPos3f(axis_len + 0.12, 0.0, 0.0)
            glutBitmapCharacter(self.label_font, ord('x'))

            glColor3f(0.0, 1.0, 0.0)
            glRasterPos3f(0.0, axis_len + 0.12, 0.0)
            glutBitmapCharacter(self.label_font, ord('y'))

            glColor3f(1.0, 0.0, 0.0)
            glRasterPos3f(0.0, 0.0, axis_len + 0.12)
            glutBitmapCharacter(self.label_font, ord('z'))

        glEnable(GL_LIGHTING)

        glPopMatrix()
        glMatrixMode(GL_PROJECTION)
        glPopMatrix()
        glMatrixMode(GL_MODELVIEW)

        glViewport(int(prev_viewport[0]), int(prev_viewport[1]), int(prev_viewport[2]), int(prev_viewport[3]))

    def on_window_resize(self, window, width, height):
        if width <= 0:
            width = 1
        if height <= 0:
            height = 1
        self.display = [width, height]
        glViewport(0, 0, width, height)
        self.set_projection()

    def on_mouse_button(self, window, button, action, mods):
        if button == glfw.MOUSE_BUTTON_LEFT:
            if action == glfw.PRESS:
                self.mouse_down = True
                self.last_mouse_pos = glfw.get_cursor_pos(window)
            elif action == glfw.RELEASE:
                self.mouse_down = False
                self.last_mouse_pos = None
        elif button == glfw.MOUSE_BUTTON_RIGHT:
            if action == glfw.PRESS:
                self.right_mouse_down = True
                self.last_mouse_pos = glfw.get_cursor_pos(window)
            elif action == glfw.RELEASE:
                self.right_mouse_down = False
                self.last_mouse_pos = None

    def on_cursor_pos(self, window, xpos, ypos):
        if self.last_mouse_pos is None:
            self.last_mouse_pos = (xpos, ypos)
            return

        last_x, last_y = self.last_mouse_pos
        dx = xpos - last_x
        dy = ypos - last_y

        if self.mouse_down:
            self.rotation_angle_x += dy * 0.2
            self.rotation_angle_y += dx * 0.2
        elif self.right_mouse_down:
            self.translation_vector[0] += dx * 0.01
            self.translation_vector[1] -= dy * 0.01

        self.last_mouse_pos = (xpos, ypos)

    def on_scroll(self, window, xoffset, yoffset):
        if self.visualization_style == "Perspective":
            self.translation_vector[2] += yoffset * 0.1
            self._clamp_perspective_zoom()
            if self.glfw_initialized:
                self.set_projection()
            return

        if yoffset > 0:
            self.zoom_factor *= 0.9
        elif yoffset < 0:
            self.zoom_factor *= 1.1

        self.zoom_factor = max(self.min_zoom, min(self.zoom_factor, self.max_zoom))
        self.set_projection()

    def set_projection(self):
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()

        aspect_ratio = max(1, self.display[0]) / max(1, self.display[1])
        radius = self.get_scene_radius()
        camera_distance = self.get_camera_distance(radius)
        near_plane, far_plane = self.get_clipping_planes(radius, camera_distance)
        ortho_half_height = self.get_ortho_half_height(radius)
        ortho_half_width = ortho_half_height * aspect_ratio

        if self.visualization_style == "Perspective":
            gluPerspective(45, aspect_ratio, near_plane, far_plane)
        else:
            glOrtho(
                -ortho_half_width * self.zoom_factor,
                ortho_half_width * self.zoom_factor,
                -ortho_half_height * self.zoom_factor,
                ortho_half_height * self.zoom_factor,
                near_plane,
                far_plane,
            )

        glMatrixMode(GL_MODELVIEW)

    def draw_sphere(self, radius: float):
        if radius in self._sphere_lists:
            glCallList(self._sphere_lists[radius])
            return

        quad = gluNewQuadric()
        try:
            gluQuadricNormals(quad, GLU_SMOOTH)
            gluSphere(quad, radius, 24, 24)
        finally:
            gluDeleteQuadric(quad)

    def create_cylinder(self, start, end, radius=None, slices=24):
        if radius is None:
            radius = self.base_bond_radius * self.bond_thickness_scale_factor
        start = np.array(start, dtype=float)
        end = np.array(end, dtype=float)

        direction = end - start
        length = np.linalg.norm(direction)
        if length == 0:
            return

        direction /= length
        up = np.array([0.0, 0.0, 1.0])
        axis = np.cross(up, direction)
        dot = float(np.clip(np.dot(up, direction), -1.0, 1.0))
        angle = np.degrees(np.arccos(dot))

        axis_norm = np.linalg.norm(axis)
        if axis_norm > 1e-12:
            axis /= axis_norm
        else:
            axis = [1.0, 0.0, 0.0] if dot < 0.0 else [0.0, 0.0, 1.0]

        glPushMatrix()
        try:
            glTranslatef(*start)
            if abs(angle) > 1e-12:
                glRotatef(angle, *axis)

            quad = gluNewQuadric()
            try:
                gluQuadricNormals(quad, GLU_SMOOTH)
                gluCylinder(quad, radius, radius, length, slices, 10)
            finally:
                gluDeleteQuadric(quad)
        finally:
            glPopMatrix()

    # ------------------------------------------------------------------
    # Chemistry helpers
    # ------------------------------------------------------------------
    def get_atom_radius(self, element: str) -> float:
        if self.atom_display_style == "vdW style":
            return self.vdw_radii.get(element, self.default_vdw_radius)
        if self.atom_display_style == "CPK style":
            return self.visual_radii.get(element, self.default_visual_radius)
        return max(self.visual_radii.get(element, self.default_visual_radius) * 0.18, 0.06)

    def show_atoms_as_spheres(self) -> bool:
        return self.atom_display_style in {"CPK style", "vdW style"}

    def show_bonds_as_cylinders(self) -> bool:
        return self.atom_display_style == "CPK style"

    def show_bonds_as_lines(self) -> bool:
        return self.atom_display_style == "Line style"

    def clear_measurement_overlay(self):
        with self._state_lock:
            self.active_measurement_type = None
            self.active_measurement_indices = []
        self.measurement_result = "Measurement: "
        if self.loop and not self.loop.is_closed() and getattr(self, "measurement_label", None):
            self.loop.call_soon_threadsafe(self._update_measurement_label, self.measurement_result)

    def set_measurement_overlay(self, measure_type: str, indices: List[int]):
        with self._state_lock:
            self.active_measurement_type = measure_type
            self.active_measurement_indices = list(indices)

    def get_measurement_overlay(self) -> Tuple[Optional[str], List[int]]:
        with self._state_lock:
            return self.active_measurement_type, list(self.active_measurement_indices)

    def get_bond_distance_threshold(self, elem1: str, elem2: str) -> float:
        r1 = self.covalent_radii.get(elem1)
        r2 = self.covalent_radii.get(elem2)
        if r1 is not None and r2 is not None:
            return min(self.connection_distance, r1 + r2 + self.bond_tolerance)
        return self.connection_distance

    def _compute_bonds_for_frame(self, frame_index: int) -> List[Tuple[int, int]]:
        frame_data = self.frames[frame_index]
        coords = np.array([pos for _, pos in frame_data], dtype=float)
        elements = [element for element, _ in frame_data]
        bonds: List[Tuple[int, int]] = []

        natoms = len(frame_data)
        for i in range(natoms):
            elem_i = elements[i]
            for j in range(i + 1, natoms):
                threshold = self.get_bond_distance_threshold(elem_i, elements[j])
                diff = coords[i] - coords[j]
                distance_sq = float(np.dot(diff, diff))
                if distance_sq <= threshold * threshold:
                    bonds.append((i, j))

        return bonds

    def calculate_bonds(self, frame_index: Optional[int] = None):
        with self._state_lock:
            if not self.frames:
                self.bonds = []
                return

            if frame_index is None:
                frame_index = self.current_frame
            if not (0 <= frame_index < len(self.frames)):
                self.bonds = []
                return

            if self.bond_rendering_mode == "Static first frame":
                if self._reference_bonds is None:
                    self._reference_bonds = self._compute_bonds_for_frame(0)
                self.bonds = list(self._reference_bonds)
                return

            if self.bond_rendering_mode == "Dynamic live":
                self.bonds = self._compute_bonds_for_frame(frame_index)
                return

            key = self._bond_cache_key(frame_index)
            cached = self._bond_cache.get(key)
            if cached is not None:
                self.bonds = list(cached)
                return

            bonds = self._compute_bonds_for_frame(frame_index)
            self._bond_cache[key] = list(bonds)
            self._trim_bond_cache_if_needed()
            self.bonds = bonds

    def advance_frame(self):
        with self._state_lock:
            if not self.frames:
                return

            next_frame = self.current_frame + (self.frame_skip * self.play_direction)
            last_index = len(self.frames) - 1

            if self.loop_mode == "Rock":
                if next_frame > last_index or next_frame < 0:
                    self.play_direction *= -1
                    next_frame = self.current_frame + (self.frame_skip * self.play_direction)
                next_frame = max(0, min(next_frame, last_index))
            elif self.loop_mode == "Once":
                if next_frame > last_index:
                    next_frame = last_index
                    self.playing = False
                elif next_frame < 0:
                    next_frame = 0
                    self.playing = False
            else:
                next_frame %= len(self.frames)

            self.current_frame = next_frame
            self.molecule_data = self.frames[self.current_frame]

        self.calculate_bonds(self.current_frame)
        if not self.playing and getattr(self, "play_pause_button", None):
            self.play_pause_button.text = "Play"

    def update_frame_label(self):
        if self.loop and not self.loop.is_closed():
            self.loop.call_soon_threadsafe(self._update_slider_and_label)

    def _update_slider_and_label(self):
        if getattr(self, "frame_slider", None) and self.frames:
            self._slider_update_in_progress = True
            try:
                self.frame_slider.value = self.current_frame
            finally:
                self._slider_update_in_progress = False
            if getattr(self, "frame_input", None):
                self.frame_input.value = str(self.current_frame)
            if getattr(self, "frame_count_label", None):
                self.frame_count_label.text = f"{self.current_frame} / {len(self.frames) - 1}"
        self.refresh_measurement_display()

    def set_frame_index(self, frame_index: int):
        if not self.frames:
            return

        with self._state_lock:
            self.current_frame = max(0, min(frame_index, len(self.frames) - 1))
            self.molecule_data = self.frames[self.current_frame]

        self.calculate_bonds(self.current_frame)
        self.update_frame_label()

    def refresh_measurement_display(self):
        measure_type, measure_indices = self.get_measurement_overlay()
        value, _ = self.compute_measurement_value(measure_type, measure_indices)
        text_value = self.format_measurement_label(measure_type, measure_indices, value)
        self.measurement_result = text_value
        if self.loop and not self.loop.is_closed() and getattr(self, "measurement_label", None):
            self.loop.call_soon_threadsafe(self._update_measurement_label, text_value)

    def _update_measurement_label(self, text_value: str):
        if getattr(self, "measurement_label", None):
            self.measurement_label.value = text_value

    def get_system_center(self, frame_data: Optional[Frame] = None):
        """Compute the center of the molecular system in the current frame."""
        if frame_data is None:
            frame_data = self.get_current_frame_data()
        if not frame_data:
            return (0.0, 0.0, 0.0)

        positions = np.array([pos for (_, pos) in frame_data], dtype=float)
        center = positions.mean(axis=0)
        return tuple(center)


    def get_scene_radius(self, frame_data: Optional[Frame] = None) -> float:
        """Return a conservative radius that encloses the current molecular scene."""
        if frame_data is None:
            with self._state_lock:
                if self.frames and 0 <= self.current_frame < len(self.frames):
                    cached_radius = self._scene_radius_cache.get(self.current_frame)
                    if cached_radius is not None:
                        return cached_radius
            frame_data = self.get_current_frame_data()
        if not frame_data:
            return 2.0

        center = np.array(self.get_system_center(frame_data), dtype=float)
        max_distance = 0.0
        for element, position in frame_data:
            radius = self.get_atom_radius(element) * self.atom_scale_factor
            distance = float(np.linalg.norm(np.array(position, dtype=float) - center)) + radius
            if distance > max_distance:
                max_distance = distance

        if self.show_box:
            a, b, c = self.box_sizes
            hx, hy, hz = a / 2.0, b / 2.0, c / 2.0
            box_center = np.array(self.get_box_center(frame_data), dtype=float)
            for dx in (-hx, hx):
                for dy in (-hy, hy):
                    for dz in (-hz, hz):
                        corner = box_center + np.array((dx, dy, dz), dtype=float)
                        distance = float(np.linalg.norm(corner - center))
                        if distance > max_distance:
                            max_distance = distance

        final_radius = max(max_distance, 2.0)
        with self._state_lock:
            if self.frames and frame_data is self.frames[self.current_frame]:
                self._scene_radius_cache[self.current_frame] = final_radius
        return final_radius

    def get_active_scene_radius(self, frame_data: Optional[Frame] = None) -> float:
        if self.auto_frame_zoom:
            return self.get_scene_radius(frame_data)

        with self._state_lock:
            if self._fixed_scene_radius is not None:
                return self._fixed_scene_radius

        if frame_data is None:
            frame_data = self.get_current_frame_data()
        fixed_radius = self.get_scene_radius(frame_data)
        with self._state_lock:
            self._fixed_scene_radius = fixed_radius
        return fixed_radius

    def get_camera_distance(self, scene_radius: Optional[float] = None) -> float:
        if scene_radius is None:
            scene_radius = self.get_scene_radius()
        self.base_camera_distance = max(5.0, scene_radius * 2.6)
        min_distance = max(scene_radius + 0.8, self.base_camera_distance * 0.35)
        distance = self.base_camera_distance - self.translation_vector[2]
        return max(distance, min_distance)

    def _clamp_perspective_zoom(self, scene_radius: Optional[float] = None):
        if scene_radius is None:
            scene_radius = self.get_scene_radius()
        self.base_camera_distance = max(5.0, scene_radius * 2.6)
        max_zoom_in = self.base_camera_distance - max(scene_radius + 0.8, self.base_camera_distance * 0.35)
        self.translation_vector[2] = min(self.translation_vector[2], max_zoom_in)

    def get_clipping_planes(self, scene_radius: Optional[float] = None, camera_distance: Optional[float] = None):
        if scene_radius is None:
            scene_radius = self.get_scene_radius()
        if camera_distance is None:
            camera_distance = self.get_camera_distance(scene_radius)

        near_plane = max(0.01, camera_distance - scene_radius - 1.0)
        far_plane = camera_distance + scene_radius + 6.0
        if far_plane <= near_plane:
            far_plane = near_plane + 10.0
        return near_plane, far_plane

    def get_ortho_half_height(self, scene_radius: Optional[float] = None) -> float:
        if scene_radius is None:
            scene_radius = self.get_scene_radius()
        return max(scene_radius * 1.3, 4.0)

    def get_box_center(self, frame_data: Optional[Frame] = None):
        """Return the center used for drawing the displayed simulation box."""
        if frame_data is None:
            frame_data = self.get_current_frame_data()

        system_center = self.get_system_center(frame_data)
        if self.box_centering_mode == "Bottom at z=0":
            _, _, c = self.box_sizes
            return (system_center[0], system_center[1], c / 2.0)
        return system_center

    def draw_box(self, frame_data: Optional[Frame] = None):
        """Draw an orthorhombic box using the selected box centering mode."""
        if frame_data is None:
            frame_data = self.get_current_frame_data()

        center = self.get_box_center(frame_data)
        a, b, c = self.box_sizes
        hx, hy, hz = a / 2.0, b / 2.0, c / 2.0
        cx, cy, cz = center

        vertices = [
            (cx - hx, cy - hy, cz - hz),
            (cx + hx, cy - hy, cz - hz),
            (cx + hx, cy + hy, cz - hz),
            (cx - hx, cy + hy, cz - hz),
            (cx - hx, cy - hy, cz + hz),
            (cx + hx, cy - hy, cz + hz),
            (cx + hx, cy + hy, cz + hz),
            (cx - hx, cy + hy, cz + hz),
        ]

        edges = [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7),
        ]

        glDisable(GL_LIGHTING)
        glColor3f(1.0, 1.0, 0.0)
        glLineWidth(2.0)
        glBegin(GL_LINES)
        for edge in edges:
            for vertex in edge:
                glVertex3f(*vertices[vertex])
        glEnd()
        glEnable(GL_LIGHTING)

    def draw_bond_line(self, start, end, start_color, end_color):
        glDisable(GL_LIGHTING)
        glLineWidth(max(1.0, self.bond_thickness_scale_factor * 2.0))
        midpoint = (
            (start[0] + end[0]) * 0.5,
            (start[1] + end[1]) * 0.5,
            (start[2] + end[2]) * 0.5,
        )
        glBegin(GL_LINES)
        glColor3f(*start_color)
        glVertex3f(*start)
        glVertex3f(*midpoint)
        glColor3f(*end_color)
        glVertex3f(*midpoint)
        glVertex3f(*end)
        glEnd()
        glEnable(GL_LIGHTING)

    def update_continuous_molecule_rotation(self):
        with self._state_lock:
            axis = self.continuous_rotation_axis
            direction = self.continuous_rotation_direction
            if axis is None or direction == 0:
                self._last_continuous_rotation_time = None
                return

            now = time.time()
            if self._last_continuous_rotation_time is None:
                self._last_continuous_rotation_time = now
                return

            elapsed = min(now - self._last_continuous_rotation_time, 0.1)
            self._last_continuous_rotation_time = now
            delta = elapsed * self.continuous_rotation_speed * direction
            if axis == "x":
                self.molecule_rotation_x = (self.molecule_rotation_x + delta) % 360.0
            elif axis == "y":
                self.molecule_rotation_y = (self.molecule_rotation_y + delta) % 360.0
            elif axis == "z":
                self.molecule_rotation_z = (self.molecule_rotation_z + delta) % 360.0

    def draw_dotted_line(self, start, end, color):
        glDisable(GL_LIGHTING)
        glColor3f(*color)
        glEnable(GL_LINE_STIPPLE)
        glLineStipple(1, 0x0F0F)
        glLineWidth(max(1.0, self.bond_thickness_scale_factor * 1.5))
        glBegin(GL_LINES)
        glVertex3f(*start)
        glVertex3f(*end)
        glEnd()
        glDisable(GL_LINE_STIPPLE)
        glEnable(GL_LIGHTING)

    def get_measurement_segments(self, measure_type: Optional[str], indices: List[int]) -> List[Tuple[int, int]]:
        if measure_type == "Bond length" and len(indices) == 2:
            return [(indices[0], indices[1])]
        if measure_type == "Bond angle" and len(indices) == 3:
            return [(indices[0], indices[1]), (indices[1], indices[2])]
        if measure_type == "Dihedral angle" and len(indices) == 4:
            return [(indices[0], indices[1]), (indices[1], indices[2]), (indices[2], indices[3])]
        return []

    def get_atom_label(self, atom_index: int, frame_data: Optional[Frame] = None) -> str:
        if frame_data is None:
            frame_data = self.get_current_frame_data()
        if atom_index < 0 or atom_index >= len(frame_data):
            return str(atom_index + 1)
        element, _ = frame_data[atom_index]
        return f"{element}{atom_index + 1}"

    def compute_measurement_value(
        self,
        measure_type: Optional[str],
        indices: List[int],
        frame_data: Optional[Frame] = None,
    ) -> Tuple[Optional[object], str]:
        if frame_data is None:
            frame_data = self.get_current_frame_data()

        if not measure_type or not indices:
            return None, ""
        if any(i < 0 or i >= len(frame_data) for i in indices):
            return None, ""

        coords = [np.array(frame_data[i][1], dtype=float) for i in indices]
        try:
            if measure_type == "Bond length" and len(coords) == 2:
                value = float(np.linalg.norm(coords[1] - coords[0]))
                return value, f"d({indices[0]+1},{indices[1]+1}) = {value:.4f} A"
            if measure_type == "Bond angle" and len(coords) == 3:
                v1 = coords[0] - coords[1]
                v2 = coords[2] - coords[1]
                n1 = np.linalg.norm(v1)
                n2 = np.linalg.norm(v2)
                if n1 == 0.0 or n2 == 0.0:
                    return None, ""
                cosang = float(np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0))
                value = float(np.degrees(np.arccos(cosang)))
                return value, f"angle = {value:.3f} deg"
            if measure_type == "Dihedral angle" and len(coords) == 4:
                b0 = -(coords[1] - coords[0])
                b1 = coords[2] - coords[1]
                b2 = coords[3] - coords[2]
                b1_norm = np.linalg.norm(b1)
                if b1_norm == 0.0:
                    return None, ""
                b1 /= b1_norm
                v = b0 - np.dot(b0, b1) * b1
                w = b2 - np.dot(b2, b1) * b1
                value = float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))
                return value, f"tau = {value:.3f} deg"
            if measure_type == "Atom coordinates" and len(coords) == 1:
                x, y, z = (float(coord) for coord in coords[0])
                return (x, y, z), ""
        except Exception:
            return None, ""

        return None, ""

    def format_measurement_label(self, measure_type: Optional[str], indices: List[int], value: Optional[object]) -> str:
        if value is None or not measure_type:
            return "Measurement: "
        if measure_type == "Bond length" and len(indices) == 2:
            return f"Measurement: d({indices[0]+1},{indices[1]+1}) = {value:.4f} Å"
        if measure_type == "Bond angle" and len(indices) == 3:
            return f"Measurement: ∠({indices[0]+1},{indices[1]+1},{indices[2]+1}) = {value:.3f}°"
        if measure_type == "Dihedral angle" and len(indices) == 4:
            return f"Measurement: τ({indices[0]+1},{indices[1]+1},{indices[2]+1},{indices[3]+1}) = {value:.3f}°"
        if measure_type == "Atom coordinates" and len(indices) == 1:
            x, y, z = value
            atom_label = self.get_atom_label(indices[0])
            return f"Atom coordinates: {atom_label} {x:.4f} {y:.4f} {z:.4f}"
        return "Measurement: "

    def get_measurement_text_position(
        self,
        measure_type: Optional[str],
        indices: List[int],
        frame_data: Frame,
    ) -> Optional[Tuple[float, float, float]]:
        if not measure_type or not indices or any(i < 0 or i >= len(frame_data) for i in indices):
            return None

        coords = [np.array(frame_data[i][1], dtype=float) for i in indices]
        if measure_type == "Bond length" and len(coords) == 2:
            pos = (coords[0] + coords[1]) * 0.5
        elif measure_type == "Bond angle" and len(coords) == 3:
            pos = coords[1]
        elif measure_type == "Dihedral angle" and len(coords) == 4:
            pos = (coords[1] + coords[2]) * 0.5
        else:
            return None

        offset = np.array([0.12, 0.12, 0.12], dtype=float)
        final_pos = pos + offset
        return (float(final_pos[0]), float(final_pos[1]), float(final_pos[2]))

    def _vector_between(self, p1, p2):
        return np.array(p2, dtype=float) - np.array(p1, dtype=float)

    def measure_bond_length(self, atom1: int, atom2: int) -> float:
        frame_data = self.get_current_frame_data()
        p1 = np.array(frame_data[atom1][1], dtype=float)
        p2 = np.array(frame_data[atom2][1], dtype=float)
        return float(np.linalg.norm(p2 - p1))

    def measure_bond_angle(self, atom1: int, atom2: int, atom3: int) -> float:
        frame_data = self.get_current_frame_data()
        p1 = np.array(frame_data[atom1][1], dtype=float)
        p2 = np.array(frame_data[atom2][1], dtype=float)
        p3 = np.array(frame_data[atom3][1], dtype=float)
        v1 = p1 - p2
        v2 = p3 - p2
        n1 = np.linalg.norm(v1)
        n2 = np.linalg.norm(v2)
        if n1 == 0.0 or n2 == 0.0:
            raise ValueError('Cannot compute angle with zero-length bond.')
        cosang = float(np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0))
        return float(np.degrees(np.arccos(cosang)))

    def measure_dihedral_angle(self, atom1: int, atom2: int, atom3: int, atom4: int) -> float:
        frame_data = self.get_current_frame_data()
        p0 = np.array(frame_data[atom1][1], dtype=float)
        p1 = np.array(frame_data[atom2][1], dtype=float)
        p2 = np.array(frame_data[atom3][1], dtype=float)
        p3 = np.array(frame_data[atom4][1], dtype=float)

        b0 = -(p1 - p0)
        b1 = p2 - p1
        b2 = p3 - p2

        b1_norm = np.linalg.norm(b1)
        if b1_norm == 0.0:
            raise ValueError('Cannot compute dihedral with zero-length central bond.')
        b1 /= b1_norm

        v = b0 - np.dot(b0, b1) * b1
        w = b2 - np.dot(b2, b1) * b1
        x = np.dot(v, w)
        y = np.dot(np.cross(b1, v), w)
        return float(np.degrees(np.arctan2(y, x)))

    def render_frame(self, frame_data: Frame, bonds: List[Tuple[int, int]]):
        measure_type, measure_indices = self.get_measurement_overlay()
        highlighted_atoms = set(i for i in measure_indices if 0 <= i < len(frame_data))
        _, measurement_overlay_text = self.compute_measurement_value(
            measure_type, measure_indices, frame_data
        )

        if self.show_bonds_as_lines():
            for i, j in bonds:
                start_element, start = frame_data[i]
                end_element, end = frame_data[j]
                start_color = self.atom_colors.get(start_element, (0.6, 0.6, 0.6))
                end_color = self.atom_colors.get(end_element, (0.6, 0.6, 0.6))
                self.draw_bond_line(start, end, start_color, end_color)

        for idx, (element, position) in enumerate(frame_data, start=1):
            atom_index = idx - 1
            color = (
                self.measurement_highlight_color
                if atom_index in highlighted_atoms
                else self.atom_colors.get(element, (0.6, 0.6, 0.6))
            )
            radius = self.get_atom_radius(element) * self.atom_scale_factor
            if self.show_atoms_as_spheres() or atom_index in highlighted_atoms:
                glPushMatrix()
                glTranslatef(*position)
                glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, (*color, 1.0))
                highlight_radius = radius if self.show_atoms_as_spheres() else max(radius, 0.18)
                if atom_index in highlighted_atoms:
                    highlight_radius = max(highlight_radius, radius * 1.15)
                self.draw_sphere(highlight_radius)
                glPopMatrix()

            if not self.fast_playback_mode:
                label_parts = []
                if self.show_atom_symbols:
                    label_parts.append(element)
                if self.show_atom_numbers:
                    label_parts.append(str(idx))
                if label_parts:
                    label_pos = (
                        position[0] + radius * 0.45,
                        position[1] + radius * 0.45,
                        position[2] + radius * 0.45,
                    )
                    self._draw_text_3d(label_pos, '-'.join(label_parts))

        if self.show_bonds_as_cylinders():
            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, (0.5, 0.5, 0.5, 1.0))
            cylinder_slices = 12 if self.fast_playback_mode else 24
            for i, j in bonds:
                start = frame_data[i][1]
                end = frame_data[j][1]
                self.create_cylinder(start, end, slices=cylinder_slices)

        for i, j in self.get_measurement_segments(measure_type, measure_indices):
            if 0 <= i < len(frame_data) and 0 <= j < len(frame_data):
                self.draw_dotted_line(
                    frame_data[i][1],
                    frame_data[j][1],
                    self.measurement_line_color,
                )
        measurement_text_position = self.get_measurement_text_position(
            measure_type, measure_indices, frame_data
        )
        if measurement_overlay_text and measurement_text_position is not None:
            self._draw_text_3d(
                measurement_text_position,
                measurement_overlay_text,
            )

    def main_loop(self):
        try:
            self.init_glfw()
            with self._state_lock:
                if self.frames:
                    self.molecule_data = self.frames[self.current_frame]
            self.calculate_bonds(self.current_frame)
            self.update_frame_label()
        except Exception as exc:
            self.set_status_message(f"OpenGL/GLFW error: {exc}")
            return

        try:
            while not glfw.window_should_close(self.window):
                if self.playing and (time.time() - self.last_frame_time) >= self.update_delay:
                    self.advance_frame()
                    self.update_frame_label()
                    self.last_frame_time = time.time()

                if self.consume_projection_update_request():
                    self.set_projection()

                self.update_continuous_molecule_rotation()

                with self._state_lock:
                    frame_data = self.get_current_frame_data()
                    bonds = list(self.bonds)

                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
                glLoadIdentity()

                scene_center = self.get_system_center(frame_data)
                scene_radius = self.get_active_scene_radius(frame_data)
                camera_distance = self.get_camera_distance(scene_radius)

                glTranslatef(0.0, 0.0, -camera_distance)
                glTranslatef(self.translation_vector[0], self.translation_vector[1], 0.0)
                glRotatef(self.rotation_angle_x, 1, 0, 0)
                glRotatef(self.rotation_angle_y, 0, 1, 0)
                glTranslatef(-scene_center[0], -scene_center[1], -scene_center[2])

                if frame_data:
                    glPushMatrix()
                    glTranslatef(scene_center[0], scene_center[1], scene_center[2])
                    glRotatef(self.molecule_rotation_x, 1, 0, 0)
                    glRotatef(self.molecule_rotation_y, 0, 1, 0)
                    glRotatef(self.molecule_rotation_z, 0, 0, 1)
                    glTranslatef(-scene_center[0], -scene_center[1], -scene_center[2])
                    self.render_frame(frame_data, bonds)
                    if self.show_box and not self.fast_playback_mode:
                        self.draw_box(frame_data)
                    glPopMatrix()

                if not self.fast_playback_mode:
                    self.draw_cartesian_axes_overlay()

                glfw.swap_buffers(self.window)
                glfw.poll_events()
                time.sleep(0.01)
        finally:
            try:
                self._destroy_gl_resources()
            finally:
                if self.window is not None:
                    glfw.destroy_window(self.window)
                    self.window = None
                glfw.terminate()
                self.glfw_initialized = False
                self._render_thread = None

    # ------------------------------------------------------------------
    # XYZ reader
    # ------------------------------------------------------------------
    async def load_xyz_trajectory(self, file_path: str):
        """Load all frames from an XYZ trajectory file using streaming parsing."""
        with self._state_lock:
            self.frames = []
            self.frame_offsets = []
            self.frame_elements = []
            self.current_frame = 0
            self.molecule_data = []
            self.bonds = []
            self.active_measurement_type = None
            self.active_measurement_indices = []
        self.invalidate_bond_cache()
        self.invalidate_scene_cache()

        count_frames = 0
        valid_atoms_first_frame: Optional[List[str]] = None

        if self.loop is None:
            self.loop = asyncio.get_event_loop()

        self.set_status_message("Loading trajectory...")

        try:
            with open(file_path, "r", encoding="utf-8") as f:
                while True:
                    frame_offset = f.tell()
                    header = f.readline()
                    if not header:
                        break

                    header = header.strip()
                    if not header:
                        continue

                    try:
                        num_atoms = int(header)
                    except ValueError:
                        continue

                    comment = f.readline()
                    if comment == "":
                        break

                    frame_data: Frame = []
                    malformed_frame = False
                    for _ in range(num_atoms):
                        atom_line = f.readline()
                        if atom_line == "":
                            malformed_frame = True
                            break

                        parts = atom_line.split()
                        if len(parts) < 4:
                            malformed_frame = True
                            break

                        element = parts[0]
                        try:
                            x, y, z = map(float, parts[1:4])
                        except ValueError:
                            malformed_frame = True
                            break

                        frame_data.append((element, (x, y, z)))

                    if malformed_frame:
                        break

                    if frame_data:
                        with self._state_lock:
                            self.frames.append(frame_data)
                            self.frame_offsets.append(frame_offset)
                        if valid_atoms_first_frame is None:
                            valid_atoms_first_frame = [element for element, _ in frame_data]
                        count_frames += 1

                    if (count_frames % 100 == 0):
                        self.set_status_message(f"Loading trajectory... {count_frames} frames")
                        await asyncio.sleep(0)
        except OSError as exc:
            self.set_status_message(f"I/O error while reading file: {exc}")
            raise

        with self._state_lock:
            if self.frames:
                self.current_frame = 0
                self.molecule_data = self.frames[0]
                self.frame_elements = valid_atoms_first_frame or []
            else:
                self.current_frame = 0
                self.molecule_data = []
                self.frame_elements = []
                self.bonds = []

        if self.frames:
            self.calculate_bonds(0)
            self.update_frame_label()
            self.set_status_message(f"{count_frames} frames loaded.")
        else:
            self.set_status_message("No valid frames found in XYZ file.")


class MolecularViewerUI(MolecularViewer):
    def __init__(self, *args, **kwargs):
        # Accept and ignore any args Toga might pass
        super().__init__()
        self.loop = asyncio.get_event_loop()

        # UI elements
        self.frame_slider = None
        self.frame_skip_input = None
        self.loop_mode_selection = None
        self.zoom_playback_switch = None
        self.frame_input = None
        self.playback_speed_input = None
        self.play_pause_button = None
        self.frame_count_label = None
        self.textInput_file = None
        self.textInput_length = None
        self.visualization_selection = None
        self.atom_display_style_selection = None
        self.atom_scale_input = None
        self.bond_thickness_input = None
        self.rotation_x_input = None
        self.rotation_y_input = None
        self.rotation_z_input = None
        self.atom_numbers_switch = None
        self.atom_symbols_switch = None
        self.measure_indices_input = None
        self.measure_type_selection = None
        self.measurement_label = None

        # Box sizes UI
        self.box_a_input = None
        self.box_b_input = None
        self.box_c_input = None
        self.box_centering_selection = None
        self.box_visibility_switch = None
        self.fast_playback_switch = None
        self.bond_rendering_selection = None

        self.layout_main_window()

    def layout_main_window(self):
        self.main_window = toga.Window(
            title="Molecular Viewer",
            size=(560, 540),
        )

        main_box = toga.Box(style=Pack(direction=COLUMN, margin=12))
        title_label = toga.Label(
            "3D Molecular Viewer",
            style=Pack(font_size=16, font_weight="bold", text_align=LEFT, margin=(0, 0, 6, 0)),
        )

        file_box = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 5, 0)))
        file_label = toga.Label(
            "Select XYZ file:",
            style=Pack(margin=(0, 0, 5, 5), text_align=LEFT, width=110),
        )
        self.textInput_file = toga.TextInput(
            placeholder="Click Browse to select XYZ file",
            style=Pack(flex=1, margin=(5, 5)),
        )
        browse_button = toga.Button(
            "Browse",
            on_press=self.browse_file,
            style=Pack(margin=5, width=72),
        )
        file_box.add(file_label)
        file_box.add(self.textInput_file)
        file_box.add(browse_button)

        length_box = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 5, 0)))
        length_label = toga.Label(
            "Set max bond length:",
            style=Pack(margin=(0, 0, 5, 5), text_align=LEFT, width=110),
        )
        self.textInput_length = toga.TextInput(
            placeholder="Upper limit for bond length",
            style=Pack(flex=1, margin=(5, 5)),
        )
        set_length_button = toga.Button(
            "Set bond length",
            on_press=self.set_connection_distance,
            style=Pack(margin=5, width=90),
        )
        length_box.add(length_label)
        length_box.add(self.textInput_length)
        length_box.add(set_length_button)

        style_box = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 5, 0)))
        style_label = toga.Label(
            "Visualization style:",
            style=Pack(margin=(0, 0, 5, 5), text_align=LEFT, width=110),
        )
        self.visualization_selection = toga.Selection(
            items=["Orthographic", "Perspective"],
            style=Pack(flex=1, margin=(5, 5)),
            on_change=self.set_visualization_style,
        )
        self.visualization_selection.value = "Orthographic"
        style_box.add(style_label)
        style_box.add(self.visualization_selection)

        atom_style_box = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 4, 0)))
        atom_style_label = toga.Label(
            "Atom display style:",
            style=Pack(margin=(0, 0, 5, 5), text_align=LEFT, width=110),
        )
        self.atom_display_style_selection = toga.Selection(
            items=["Line style", "CPK style", "vdW style"],
            style=Pack(width=120, margin=(5, 5)),
            on_change=self.set_atom_display_style,
        )
        self.atom_display_style_selection.value = self.atom_display_style
        atom_style_box.add(atom_style_label)
        atom_style_box.add(self.atom_display_style_selection)

        scale_box = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 4, 0)))
        atom_scale_label = toga.Label(
            "Atom size scale:",
            style=Pack(margin=(0, 0, 5, 5), text_align=LEFT, width=110),
        )
        self.atom_scale_input = toga.TextInput(
            placeholder="e.g. 1.0",
            style=Pack(width=72, margin=(5, 5)),
        )
        self.atom_scale_input.value = "1.0"
        set_atom_scale_button = toga.Button(
            "Set atom scale",
            on_press=self.set_atom_scale_factor,
            style=Pack(margin=5, width=90),
        )
        bond_thickness_label = toga.Label(
            "Bond thickness:",
            style=Pack(margin=(0, 0, 5, 5), text_align=LEFT, width=110),
        )
        self.bond_thickness_input = toga.TextInput(
            placeholder="e.g. 1.0",
            style=Pack(width=72, margin=(5, 5)),
        )
        self.bond_thickness_input.value = "1.0"
        set_bond_thickness_button = toga.Button(
            "Set bond scale",
            on_press=self.set_bond_thickness_scale_factor,
            style=Pack(margin=5, width=90),
        )
        scale_box.add(atom_scale_label)
        scale_box.add(self.atom_scale_input)
        scale_box.add(set_atom_scale_button)
        scale_box.add(bond_thickness_label)
        scale_box.add(self.bond_thickness_input)
        scale_box.add(set_bond_thickness_button)

        rotation_box = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 4, 0)))
        rotation_label = toga.Label(
            "Rotate:",
            style=Pack(margin=(0, 0, 5, 5), text_align=LEFT, width=58),
        )
        rotation_x_label = toga.Label(
            "X:",
            style=Pack(margin=(0, 0, 5, 2), text_align=LEFT, width=16),
        )
        self.rotation_x_input = toga.TextInput(
            value="0.0",
            placeholder="deg",
            style=Pack(width=46, margin=(2, 2)),
        )
        rotate_x_negative_button = toga.Button(
            "<",
            on_press=lambda widget: self.toggle_continuous_molecule_rotation("x", -1),
            style=Pack(margin=2, width=26),
        )
        rotate_x_positive_button = toga.Button(
            ">",
            on_press=lambda widget: self.toggle_continuous_molecule_rotation("x", 1),
            style=Pack(margin=2, width=26),
        )
        rotation_y_label = toga.Label(
            "Y:",
            style=Pack(margin=(0, 0, 5, 8), text_align=LEFT, width=16),
        )
        self.rotation_y_input = toga.TextInput(
            value="0.0",
            placeholder="deg",
            style=Pack(width=46, margin=(2, 2)),
        )
        rotate_y_negative_button = toga.Button(
            "<",
            on_press=lambda widget: self.toggle_continuous_molecule_rotation("y", -1),
            style=Pack(margin=2, width=26),
        )
        rotate_y_positive_button = toga.Button(
            ">",
            on_press=lambda widget: self.toggle_continuous_molecule_rotation("y", 1),
            style=Pack(margin=2, width=26),
        )
        rotation_z_label = toga.Label(
            "Z:",
            style=Pack(margin=(0, 0, 5, 8), text_align=LEFT, width=16),
        )
        self.rotation_z_input = toga.TextInput(
            value="0.0",
            placeholder="deg",
            style=Pack(width=46, margin=(2, 2)),
        )
        rotate_z_negative_button = toga.Button(
            "<",
            on_press=lambda widget: self.toggle_continuous_molecule_rotation("z", -1),
            style=Pack(margin=2, width=26),
        )
        rotate_z_positive_button = toga.Button(
            ">",
            on_press=lambda widget: self.toggle_continuous_molecule_rotation("z", 1),
            style=Pack(margin=2, width=26),
        )
        apply_rotation_button = toga.Button(
            "Apply",
            on_press=self.apply_molecule_rotation,
            style=Pack(margin=(2, 2, 2, 8), width=52),
        )
        reset_rotation_button = toga.Button(
            "Reset",
            on_press=self.reset_molecule_rotation,
            style=Pack(margin=2, width=52),
        )
        rotation_box.add(rotation_label)
        rotation_box.add(rotation_x_label)
        rotation_box.add(self.rotation_x_input)
        rotation_box.add(rotate_x_negative_button)
        rotation_box.add(rotate_x_positive_button)
        rotation_box.add(rotation_y_label)
        rotation_box.add(self.rotation_y_input)
        rotation_box.add(rotate_y_negative_button)
        rotation_box.add(rotate_y_positive_button)
        rotation_box.add(rotation_z_label)
        rotation_box.add(self.rotation_z_input)
        rotation_box.add(rotate_z_negative_button)
        rotation_box.add(rotate_z_positive_button)
        rotation_box.add(apply_rotation_button)
        rotation_box.add(reset_rotation_button)

        labels_box = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 4, 0)))
        self.atom_numbers_switch = toga.Switch("Atom numbers", on_change=self.toggle_atom_numbers)
        self.atom_numbers_switch.value = False
        self.atom_symbols_switch = toga.Switch("Atomic symbols", on_change=self.toggle_atom_symbols)
        self.atom_symbols_switch.value = False
        labels_box.add(self.atom_numbers_switch)
        labels_box.add(self.atom_symbols_switch)

        measure_box = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 4, 0)))
        measure_label = toga.Label(
            "Measure:",
            style=Pack(margin=(0, 0, 4, 4), text_align=LEFT, width=70),
        )
        self.measure_type_selection = toga.Selection(
            items=["Bond length", "Bond angle", "Dihedral angle", "Atom coordinates"],
            style=Pack(width=150, margin=(4, 4)),
            on_change=self.update_measurement_input_hint,
        )
        self.measure_type_selection.value = "Bond length"
        self.measure_indices_input = toga.TextInput(
            placeholder="1,2 or 1,2,3 or 1,2,3,4",
            style=Pack(flex=1, margin=(4, 4)),
        )
        measure_button = toga.Button(
            "Measure",
            on_press=self.run_measurement,
            style=Pack(margin=4, width=78),
        )
        measure_box.add(measure_label)
        measure_box.add(self.measure_type_selection)
        measure_box.add(self.measure_indices_input)
        measure_box.add(measure_button)

        self.measurement_label = toga.TextInput(
            value="Measurement: ",
            readonly=True,
            style=Pack(margin=(2, 0, 6, 0), font_size=10),
        )

        box_box = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 5, 0)))
        box_a_label = toga.Label(
            "Box size a:",
            style=Pack(margin=(0, 0, 5, 5), text_align=LEFT, width=72),
        )
        self.box_a_input = toga.TextInput(placeholder="a", style=Pack(flex=1, margin=(5, 5)))
        box_b_label = toga.Label(
            "Box size b:",
            style=Pack(margin=(0, 0, 5, 5), text_align=LEFT, width=72),
        )
        self.box_b_input = toga.TextInput(placeholder="b", style=Pack(flex=1, margin=(5, 5)))
        box_c_label = toga.Label(
            "Box size c:",
            style=Pack(margin=(0, 0, 5, 5), text_align=LEFT, width=72),
        )
        self.box_c_input = toga.TextInput(placeholder="c", style=Pack(flex=1, margin=(5, 5)))
        set_box_button = toga.Button(
            "Set Box Sizes",
            on_press=self.set_box_sizes,
            style=Pack(margin=5, width=90),
        )
        box_box.add(box_a_label)
        box_box.add(self.box_a_input)
        box_box.add(box_b_label)
        box_box.add(self.box_b_input)
        box_box.add(box_c_label)
        box_box.add(self.box_c_input)
        box_box.add(set_box_button)

        box_centering_box = toga.Box(
            style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 5, 0))
        )
        box_centering_label = toga.Label(
            "Box center:",
            style=Pack(margin=(0, 0, 5, 5), text_align=LEFT, width=110),
        )
        self.box_centering_selection = toga.Selection(
            items=["Geometric center", "Bottom at z=0"],
            style=Pack(flex=1, margin=(5, 5)),
            on_change=self.set_box_centering_mode,
        )
        self.box_centering_selection.value = self.box_centering_mode
        box_centering_box.add(box_centering_label)
        box_centering_box.add(self.box_centering_selection)

        box_visibility_box = toga.Box(
            style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 5, 0))
        )
        box_visibility_label = toga.Label(
            "Show Box:",
            style=Pack(margin=(0, 0, 5, 5), text_align=LEFT, width=110),
        )
        self.box_visibility_switch = toga.Switch("Show Box", on_change=self.set_box_visibility)
        self.box_visibility_switch.value = False
        box_visibility_box.add(box_visibility_label)
        box_visibility_box.add(self.box_visibility_switch)

        nav_box = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 5, 0)))
        first_button = toga.Button("|<", on_press=self.go_to_first_frame, style=Pack(margin=3, width=34))
        prev_button = toga.Button("<<", on_press=self.step_backward, style=Pack(margin=3, width=34))
        self.frame_input = toga.TextInput(value="0", style=Pack(width=64, margin=(3, 5)))
        go_button = toga.Button("Go", on_press=self.go_to_frame_from_input, style=Pack(margin=3, width=34))
        self.frame_slider = toga.Slider(
            min=0,
            max=1,
            value=0,
            style=Pack(flex=1, margin=5),
            on_change=self.set_current_frame,
        )
        next_button = toga.Button(">>", on_press=self.step_forward, style=Pack(margin=3, width=34))
        last_button = toga.Button(">|", on_press=self.go_to_last_frame, style=Pack(margin=3, width=34))
        self.frame_count_label = toga.Label(
            "0 / 0",
            style=Pack(margin=(0, 0, 5, 5), text_align=LEFT, width=72),
        )
        nav_box.add(first_button)
        nav_box.add(prev_button)
        nav_box.add(self.frame_input)
        nav_box.add(go_button)
        nav_box.add(self.frame_slider)
        nav_box.add(next_button)
        nav_box.add(last_button)
        nav_box.add(self.frame_count_label)

        playback_box = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 5, 0)))
        loop_label = toga.Label(
            "Loop:",
            style=Pack(margin=(0, 0, 5, 4), text_align=LEFT, width=42),
        )
        self.loop_mode_selection = toga.Selection(
            items=["Loop", "Once", "Rock"],
            style=Pack(width=88, margin=(0, 6, 0, 0)),
            on_change=self.set_loop_mode,
        )
        self.loop_mode_selection.value = self.loop_mode
        zoom_label = toga.Label(
            "zoom",
            style=Pack(margin=(0, 0, 5, 12), text_align=LEFT, width=38),
        )
        self.zoom_playback_switch = toga.Switch("", on_change=self.toggle_auto_frame_zoom)
        self.zoom_playback_switch.value = self.auto_frame_zoom
        self.play_pause_button = toga.Button(
            "Play",
            on_press=self.toggle_play_pause,
            style=Pack(margin=4, width=54),
        )
        step_label = toga.Label(
            "step",
            style=Pack(margin=(0, 0, 5, 10), text_align=LEFT, width=34),
        )
        step_back_button = toga.Button("<", on_press=self.decrease_frame_step, style=Pack(margin=3, width=28))
        self.frame_skip_input = toga.TextInput(
            value=str(self.frame_skip),
            placeholder="Step",
            style=Pack(width=50, margin=(3, 5)),
        )
        step_forward_button = toga.Button(">", on_press=self.increase_frame_step, style=Pack(margin=3, width=28))
        playback_speed_label = toga.Label(
            "Delay (s):",
            style=Pack(margin=(0, 0, 5, 10), text_align=LEFT, width=62),
        )
        self.playback_speed_input = toga.TextInput(
            value=f"{self.update_delay:.2f}",
            placeholder="Seconds between frames",
            style=Pack(width=70, margin=(5, 5)),
        )
        set_speed_button = toga.Button(
            "Set Speed",
            on_press=self.set_playback_delay,
            style=Pack(margin=5, width=78),
        )
        playback_box.add(loop_label)
        playback_box.add(self.loop_mode_selection)
        playback_box.add(zoom_label)
        playback_box.add(self.zoom_playback_switch)
        playback_box.add(step_label)
        playback_box.add(step_back_button)
        playback_box.add(self.frame_skip_input)
        playback_box.add(step_forward_button)
        playback_box.add(playback_speed_label)
        playback_box.add(self.playback_speed_input)
        playback_box.add(set_speed_button)
        playback_box.add(self.play_pause_button)

        skip_box = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 5, 0)))
        skip_label = toga.Label(
            "Frame Step:",
            style=Pack(margin=(0, 0, 5, 5), text_align=LEFT, width=72),
        )
        set_skip_button = toga.Button(
            "Set Step",
            on_press=self.set_frame_skip,
            style=Pack(margin=5, width=74),
        )
        skip_box.add(skip_label)
        skip_box.add(self.frame_skip_input)
        skip_box.add(set_skip_button)

        performance_box = toga.Box(
            style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 5, 0))
        )
        fast_playback_label = toga.Label(
            "Fast Playback:",
            style=Pack(margin=(0, 0, 5, 5), text_align=LEFT, width=110),
        )
        self.fast_playback_switch = toga.Switch(
            "Reduce overlays during playback",
            on_change=self.toggle_fast_playback_mode,
        )
        self.fast_playback_switch.value = self.fast_playback_mode
        bond_rendering_label = toga.Label(
            "Bond mode:",
            style=Pack(margin=(0, 0, 5, 12), text_align=LEFT, width=80),
        )
        self.bond_rendering_selection = toga.Selection(
            items=["Static first frame", "Dynamic cached", "Dynamic live"],
            style=Pack(width=135, margin=(5, 5)),
            on_change=self.set_bond_rendering_mode,
        )
        self.bond_rendering_selection.value = self.bond_rendering_mode
        performance_box.add(fast_playback_label)
        performance_box.add(self.fast_playback_switch)
        performance_box.add(bond_rendering_label)
        performance_box.add(self.bond_rendering_selection)

        action_box = toga.Box(style=Pack(direction=ROW, align_items=CENTER, margin=(0, 0, 5, 0)))
        self.display_button = toga.Button(
            "Display Molecule/Trajectory",
            on_press=self.open_opengl_window,
            style=Pack(margin=8, align_items=CENTER),
        )
        save_frame_button = toga.Button(
            "Save Current Frame XYZ",
            on_press=self.save_current_frame_xyz,
            style=Pack(margin=8, align_items=CENTER),
        )
        action_box.add(self.display_button)
        action_box.add(save_frame_button)

        self.loading_label = toga.Label(
            " ",
            style=Pack(margin=10, font_size=18, font_weight="bold", text_align=CENTER),
        )

        main_box.add(title_label)
        main_box.add(file_box)
        main_box.add(length_box)
        main_box.add(style_box)
        main_box.add(atom_style_box)
        main_box.add(scale_box)
        main_box.add(rotation_box)
        main_box.add(labels_box)
        main_box.add(measure_box)
        main_box.add(self.measurement_label)
        main_box.add(box_box)
        main_box.add(box_centering_box)
        main_box.add(box_visibility_box)
        main_box.add(nav_box)
        main_box.add(playback_box)
        main_box.add(skip_box)
        main_box.add(performance_box)
        main_box.add(action_box)
        main_box.add(self.loading_label)

        self.main_window.content = main_box
        self.main_window.show()

    async def _show_error(self, title: str, message: str):
        await self.main_window.dialog(toga.ErrorDialog(title, message))

    async def browse_file(self, widget):
        try:
            selected = None

            selected = await self.main_window.dialog(
                toga.OpenFileDialog("Open XYZ trajectory", multiple_select=False, file_types=["*.*"])
            )

            if not selected:
                return

            if isinstance(selected, (list, tuple)):
                selected = selected[0] if selected else None
            if not selected:
                return

            file_path = str(selected)
            self.textInput_file.value = file_path
            await self.load_xyz_trajectory(file_path)

            if self.frames:
                self.frame_slider.max = len(self.frames) - 1
                self.frame_slider.value = 0
                self.frame_input.value = str(self.current_frame)
                self.frame_count_label.text = f"{self.current_frame} / {len(self.frames) - 1}"
                self.play_pause_button.text = "Play"
                self.playing = False
        except Exception as exc:
            await self._show_error("File Error", f"Could not open XYZ file.\n{exc}")

    async def set_connection_distance(self, widget):
        try:
            distance = float(self.textInput_length.value.strip())
            if distance <= 0:
                raise ValueError

            self.connection_distance = distance
            self.invalidate_bond_cache()
            self.calculate_bonds()
            self.set_status_message(f"Maximum bond length set to {distance:.3f} Å.")
        except (TypeError, ValueError, AttributeError):
            await self._show_error(
                "Invalid Input", "Please enter a valid positive number for connection distance."
            )

    def set_visualization_style(self, widget):
        self.visualization_style = widget.value
        self.zoom_factor = 1.0
        self._clamp_perspective_zoom()
        if self.glfw_initialized:
            self.request_projection_update()

    def set_current_frame(self, widget):
        if not self.frames or self._slider_update_in_progress:
            return

        self.set_frame_index(int(round(widget.value)))

    def toggle_play_pause(self, widget):
        if not self.frames:
            return

        self.playing = not self.playing
        if self.loop_mode == "Rock" and self.current_frame in (0, len(self.frames) - 1):
            self.play_direction = -1 if self.current_frame == len(self.frames) - 1 else 1
        self.play_pause_button.text = "Pause" if self.playing else "Play"
        self.last_frame_time = time.time()

    def set_loop_mode(self, widget):
        self.loop_mode = self.loop_mode_selection.value or "Loop"
        self.set_status_message(f"Loop mode set to {self.loop_mode}.")

    def toggle_auto_frame_zoom(self, widget):
        self.auto_frame_zoom = bool(self.zoom_playback_switch.value)
        if not self.auto_frame_zoom:
            self.invalidate_scene_cache()
        if self.glfw_initialized:
            self.request_projection_update()
        self.set_status_message(
            "Auto zoom enabled for frame playback."
            if self.auto_frame_zoom
            else "Auto zoom disabled for steadier playback."
        )

    async def set_playback_delay(self, widget):
        try:
            delay = float(self.playback_speed_input.value.strip())
            if delay <= 0:
                raise ValueError

            self.update_delay = delay
            self.set_status_message(f"Playback delay set to {delay:.3f} s per frame.")
        except (TypeError, ValueError, AttributeError):
            await self._show_error("Invalid Input", "Playback delay must be a positive number.")

    async def set_atom_scale_factor(self, widget):
        try:
            scale = float(self.atom_scale_input.value.strip())
            if scale <= 0:
                raise ValueError
            self.atom_scale_factor = scale
            self.invalidate_scene_cache()
            self.set_status_message(f"Atom size scale set to {scale:.3f}.")
        except (TypeError, ValueError, AttributeError):
            await self._show_error(
                "Invalid Input", "Please enter a valid positive number for atom size scale."
            )

    def set_atom_display_style(self, widget):
        self.atom_display_style = self.atom_display_style_selection.value or "Line style"
        self.invalidate_scene_cache()
        self._clamp_perspective_zoom()
        if self.glfw_initialized:
            self.request_projection_update()
        self.set_status_message(
            f"Atom display style set to {self.atom_display_style}."
        )

    async def set_bond_thickness_scale_factor(self, widget):
        try:
            scale = float(self.bond_thickness_input.value.strip())
            if scale <= 0:
                raise ValueError
            self.bond_thickness_scale_factor = scale
            self.set_status_message(f"Bond thickness scale set to {scale:.3f}.")
        except (TypeError, ValueError, AttributeError):
            await self._show_error(
                "Invalid Input", "Please enter a valid positive number for bond thickness scale."
            )

    async def apply_molecule_rotation(self, widget):
        try:
            delta_x = float(self.rotation_x_input.value.strip())
            delta_y = float(self.rotation_y_input.value.strip())
            delta_z = float(self.rotation_z_input.value.strip())
        except (TypeError, ValueError, AttributeError):
            await self._show_error(
                "Invalid Input", "Please enter valid numeric rotation angles in degrees."
            )
            return

        with self._state_lock:
            self.molecule_rotation_x = (self.molecule_rotation_x + delta_x) % 360.0
            self.molecule_rotation_y = (self.molecule_rotation_y + delta_y) % 360.0
            self.molecule_rotation_z = (self.molecule_rotation_z + delta_z) % 360.0
            rotation_x = self.molecule_rotation_x
            rotation_y = self.molecule_rotation_y
            rotation_z = self.molecule_rotation_z

        self.set_status_message(
            "Molecular rotation set to "
            f"x={rotation_x:.3f}, "
            f"y={rotation_y:.3f}, "
            f"z={rotation_z:.3f} degrees."
        )

    def toggle_continuous_molecule_rotation(self, axis: str, direction: int):
        if axis not in {"x", "y", "z"} or direction not in {-1, 1}:
            return

        with self._state_lock:
            is_same_rotation = (
                self.continuous_rotation_axis == axis
                and self.continuous_rotation_direction == direction
            )
            if is_same_rotation:
                self.continuous_rotation_axis = None
                self.continuous_rotation_direction = 0
                self._last_continuous_rotation_time = None
                message = f"Continuous {axis.upper()} rotation stopped."
            else:
                self.continuous_rotation_axis = axis
                self.continuous_rotation_direction = direction
                self._last_continuous_rotation_time = None
                direction_label = "positive" if direction > 0 else "negative"
                message = f"Continuous {axis.upper()} rotation started in the {direction_label} direction."

        self.set_status_message(message)

    def reset_molecule_rotation(self, widget):
        with self._state_lock:
            self.molecule_rotation_x = 0.0
            self.molecule_rotation_y = 0.0
            self.molecule_rotation_z = 0.0
            self.continuous_rotation_axis = None
            self.continuous_rotation_direction = 0
            self._last_continuous_rotation_time = None

        self.rotation_x_input.value = "0.0"
        self.rotation_y_input.value = "0.0"
        self.rotation_z_input.value = "0.0"
        self.set_status_message("Molecular rotation reset.")

    async def set_frame_skip(self, widget):
        try:
            skip = int(self.frame_skip_input.value.strip())
            if skip < 1:
                raise ValueError

            self.frame_skip = skip
            self.frame_skip_input.value = str(skip)
            self.set_status_message(f"Frame step set to {skip}.")
        except (TypeError, ValueError, AttributeError):
            await self._show_error("Invalid Input", "Frame step must be an integer >= 1.")

    def go_to_first_frame(self, widget):
        self.play_direction = 1
        self.set_frame_index(0)

    def go_to_last_frame(self, widget):
        self.play_direction = -1
        if self.frames:
            self.set_frame_index(len(self.frames) - 1)

    def step_backward(self, widget):
        self.play_direction = -1
        if self.frames:
            self.set_frame_index(self.current_frame - self.frame_skip)

    def step_forward(self, widget):
        self.play_direction = 1
        if self.frames:
            self.set_frame_index(self.current_frame + self.frame_skip)

    def decrease_frame_step(self, widget):
        self.frame_skip = max(1, self.frame_skip - 1)
        self.frame_skip_input.value = str(self.frame_skip)
        self.set_status_message(f"Frame step set to {self.frame_skip}.")

    def increase_frame_step(self, widget):
        self.frame_skip += 1
        self.frame_skip_input.value = str(self.frame_skip)
        self.set_status_message(f"Frame step set to {self.frame_skip}.")

    async def go_to_frame_from_input(self, widget):
        try:
            if not self.frames:
                raise ValueError("No frames loaded.")
            frame_value = int(self.frame_input.value.strip())
            self.set_frame_index(frame_value)
        except Exception:
            await self._show_error("Invalid Input", "Enter a valid frame index within range.")

    def toggle_fast_playback_mode(self, widget):
        self.fast_playback_mode = bool(self.fast_playback_switch.value)
        self.set_status_message(
            "Fast playback mode enabled."
            if self.fast_playback_mode
            else "Fast playback mode disabled."
        )

    def set_bond_rendering_mode(self, widget):
        self.bond_rendering_mode = self.bond_rendering_selection.value or "Static first frame"
        self.invalidate_bond_cache()
        self.calculate_bonds(self.current_frame if self.frames else None)
        messages = {
            "Static first frame": "Using first-frame connectivity for playback.",
            "Dynamic cached": "Dynamic bond recalculation enabled with per-frame caching.",
            "Dynamic live": "Dynamic live bond recalculation enabled for every frame.",
        }
        self.set_status_message(messages.get(self.bond_rendering_mode, "Bond rendering mode updated."))

    async def set_box_sizes(self, widget):
        try:
            a = float(self.box_a_input.value.strip())
            b = float(self.box_b_input.value.strip())
            c = float(self.box_c_input.value.strip())
            if a <= 0 or b <= 0 or c <= 0:
                raise ValueError

            self.box_sizes = (a, b, c)
            self.invalidate_scene_cache()
            self.set_status_message(f"Box sizes set to a={a:.3f}, b={b:.3f}, c={c:.3f} Å.")
        except (TypeError, ValueError, AttributeError):
            await self._show_error(
                "Invalid Input", "Please enter valid positive numbers for the box sizes."
            )

    def set_box_centering_mode(self, widget):
        selected_mode = getattr(self.box_centering_selection, "value", None)
        if selected_mode not in {"Geometric center", "Bottom at z=0"}:
            selected_mode = "Geometric center"

        self.box_centering_mode = selected_mode
        self.invalidate_scene_cache()
        if selected_mode == "Bottom at z=0":
            self.set_status_message("Box centered in a-b plane with bottom face at z=0.")
        else:
            self.set_status_message("Box centered on the molecular geometric center.")

    def set_box_visibility(self, widget):
        self.show_box = bool(self.box_visibility_switch.value)
        self.invalidate_scene_cache()
        self.set_status_message("Box display enabled." if self.show_box else "Box display disabled.")

    def toggle_atom_numbers(self, widget):
        self.show_atom_numbers = bool(self.atom_numbers_switch.value)
        self.set_status_message(
            "Atom numbering enabled." if self.show_atom_numbers else "Atom numbering disabled."
        )

    def toggle_atom_symbols(self, widget):
        self.show_atom_symbols = bool(self.atom_symbols_switch.value)
        self.set_status_message(
            "Atomic symbols enabled." if self.show_atom_symbols else "Atomic symbols disabled."
        )

    def update_measurement_input_hint(self, widget):
        if not getattr(self, "measure_indices_input", None):
            return
        measure_type = self.measure_type_selection.value
        if measure_type == "Atom coordinates":
            self.measure_indices_input.placeholder = "Atom label, e.g. 5 or C-5"
        else:
            self.measure_indices_input.placeholder = "1,2 or 1,2,3 or 1,2,3,4"

    def _parse_measurement_indices(self, expected_count: int) -> List[int]:
        raw = (self.measure_indices_input.value or "").replace(';', ',')
        parts = [p.strip() for p in raw.split(',') if p.strip()]
        if len(parts) != expected_count:
            raise ValueError(f"Enter exactly {expected_count} atom indices.")
        indices = [int(p) - 1 for p in parts]
        if not self.frames and not self.molecule_data:
            raise ValueError("No molecular data loaded.")
        frame_data = self.get_current_frame_data()
        natoms = len(frame_data)
        if any(i < 0 or i >= natoms for i in indices):
            raise ValueError(f"Atom indices must be between 1 and {natoms}.")
        return indices

    def _parse_atom_label(self) -> int:
        raw = (self.measure_indices_input.value or "").strip()
        if not raw:
            raise ValueError("Enter one valid atom label, for example 5 or C-5.")

        frame_data = self.get_current_frame_data()
        natoms = len(frame_data)
        if natoms == 0:
            raise ValueError("No molecular data loaded.")

        compact = raw.replace(" ", "").replace("_", "-")
        if compact.isdigit():
            atom_index = int(compact) - 1
            if 0 <= atom_index < natoms:
                return atom_index
            raise ValueError(f"Atom label must identify an atom between 1 and {natoms}.")

        label_parts = compact.split("-", 1)
        if len(label_parts) == 2 and label_parts[1].isdigit():
            symbol_part, number_part = label_parts
        else:
            symbol_part = "".join(ch for ch in compact if ch.isalpha())
            number_part = "".join(ch for ch in compact if ch.isdigit())

        if symbol_part and number_part:
            atom_index = int(number_part) - 1
            if 0 <= atom_index < natoms:
                expected_symbol = symbol_part[0].upper() + symbol_part[1:].lower()
                actual_symbol = frame_data[atom_index][0]
                if actual_symbol == expected_symbol:
                    return atom_index
                raise ValueError(
                    f"Atom {atom_index + 1} is {actual_symbol}, not {expected_symbol}."
                )

        if symbol_part and not number_part:
            expected_symbol = symbol_part[0].upper() + symbol_part[1:].lower()
            matches = [
                index for index, (element, _) in enumerate(frame_data)
                if element == expected_symbol
            ]
            if len(matches) == 1:
                return matches[0]
            if len(matches) > 1:
                raise ValueError(
                    f"Multiple {expected_symbol} atoms found. Use a numbered label such as {expected_symbol}-{matches[0] + 1}."
                )

        raise ValueError("Enter one valid atom label, for example 5 or C-5.")

    async def run_measurement(self, widget):
        try:
            measure_type = self.measure_type_selection.value
            if measure_type == "Bond length":
                i, j = self._parse_measurement_indices(2)
                value = self.measure_bond_length(i, j)
                self.set_measurement_overlay(measure_type, [i, j])
                text_value = self.format_measurement_label(measure_type, [i, j], value)
            elif measure_type == "Bond angle":
                i, j, k = self._parse_measurement_indices(3)
                value = self.measure_bond_angle(i, j, k)
                self.set_measurement_overlay(measure_type, [i, j, k])
                text_value = self.format_measurement_label(measure_type, [i, j, k], value)
            elif measure_type == "Dihedral angle":
                i, j, k, l = self._parse_measurement_indices(4)
                value = self.measure_dihedral_angle(i, j, k, l)
                self.set_measurement_overlay(measure_type, [i, j, k, l])
                text_value = self.format_measurement_label(measure_type, [i, j, k, l], value)
            elif measure_type == "Atom coordinates":
                i = self._parse_atom_label()
                value, _ = self.compute_measurement_value(measure_type, [i])
                self.set_measurement_overlay(measure_type, [i])
                text_value = self.format_measurement_label(measure_type, [i], value)
            else:
                raise ValueError("Select a valid measurement type.")

            self.measurement_result = text_value
            self.measurement_label.value = text_value
            self.set_status_message(text_value)
        except Exception as exc:
            self.clear_measurement_overlay()
            await self._show_error("Measurement Error", str(exc))

    async def save_current_frame_xyz(self, widget):
        try:
            frame_data = self.get_current_frame_data()
            if not frame_data:
                await self._show_error("No Data", "No frame is available to save.")
                return

            base_name = "frame"
            input_path = (self.textInput_file.value or "").strip()
            if input_path:
                base_name = os.path.splitext(os.path.basename(input_path))[0] or "frame"
            suggested = f"{base_name}_frame_{self.current_frame:04d}.xyz"

            selected = await self.main_window.dialog(
                toga.SaveFileDialog(title="Save Current Frame as XYZ", suggested_filename=suggested)
            )
            if not selected:
                return

            output_path = str(selected)
            with open(output_path, "w", encoding="utf-8") as f:
                f.write(f"{len(frame_data)}\n")
                f.write(f"Frame {self.current_frame} exported from molecularViewer\n")
                for element, (x, y, z) in frame_data:
                    f.write(f"{element:<3s} {x: .8f} {y: .8f} {z: .8f}\n")

            await self.main_window.dialog(
                toga.InfoDialog("Frame Saved", f"Current frame saved to:\n{output_path}")
            )
        except Exception as exc:
            await self._show_error("Save Error", str(exc))

    async def open_opengl_window(self, widget):
        if self._render_thread and self._render_thread.is_alive():
            self.set_status_message("Viewer window is already open.")
            return

        if not self.frames and not self.molecule_data:
            await self._show_error("No Data", "Load an XYZ file before opening the viewer.")
            return

        self.set_status_message("Opening OpenGL viewer...")
        self._render_thread = threading.Thread(target=self.main_loop, daemon=True)
        self._render_thread.start()

def main():
    return MolecularViewerUI()
