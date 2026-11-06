# -*- encoding: utf-8 -*-
# import wx

import numpy as np
from threading import Timer
# from time import clock
# import copy

import cairo
import wx.lib.wxcairo

import gi
gi.require_version("Gtk", "3.0")
from gi.repository import GLib, Gio, Gtk, Gdk, GObject
GObject.threads_init()

# import ase
from ase.data import covalent_radii
from ase.data.colors import jmol_colors


class ClusterViewer(Gtk.Window):
    """
    :param mx: last mouse position in x
    :param my: last mouse position in y
    """
    MODE_NONE = 0b0000000
    MODE_SELECTION = 0b0000001
    MODE_SELECTION_BOX = 0b0000010
    MODE_SELECTION_APPEND = 0b0000100
    MODE_SELECTION_TOGGLE = 0b0001000
    MODE_TRANSLATION = 0b0010000
    MODE_ROTATION = 0b0100000

    def __init__(self, *args, **kwargs):
        #GTK# kwargs['style'] = wx.NO_FULL_REPAINT_ON_RESIZE | wx.CLIP_CHILDREN
        Gtk.Window.__init__(self, *args, **kwargs)
        self.drawing_area = Gtk.DrawingArea()
        self.add(self.drawing_area)

        self.ox = self.oy = 0  # offset in x and y
        self.im_ox = self.im_oy = 0  # image offset in x and y
        self.last_mouse_move_x = self.last_mouse_move_y = 0  # last mouse move
        self.mx = self.my = 0  # last mouse position
        self.theta = self.phi = 0
        self.scale = self.scale0 = 100
        self.im_factor = self.im_scale = 1
        self.atoms = None
        # self.do_rescale = False
        # self.do_center = False
        self.atoms_center_of_mass = np.zeros(3)
        self.atoms_largest_dimension = 1.0  # float, in angstrom
        self.selection = []
        self.selection_box = None
        self.__outer_margin = 0
        self.surface = None
        self.busy = False
        self.refresh_delay = 200
        self.back_buffer = None
        self.screenshot = None
        self.atom_numbers = None
        self.atom_surfaces = None
        self.atoms_sprite = None
        self.background_sprite = None

        self.mode = self.MODE_NONE

        self.colors = {
            'selection_box': (0.0, 0.4, 1.0),
            'boulding_box_line': (0.0, 0.4, 1.0, 1.0),
            'boulding_box_fill': (0.0, 0.4, 1.0, 0.3),
        }
        self.sprites_opts = {'alpha': 1, 'glow': True}

        self.light_mode = False
        self.light_mode_threshold = 2000

        self.rotation_matrix = np.identity(4)
        self.scale_matrix = np.identity(4)
        self.translation_matrix = np.identity(4)
        self.model_matrix = np.identity(4)
        self.projection_matrix = np.identity(4)
        # model to world matrix
        self.m2w_matrix = np.identity(4)
        # world to view matrix
        self.w2v_matrix = np.identity(4)
        # view to projection matrix
        viewport = (-1., 1., -1., 1., -1., 1.)
        self.v2p_matrix = self.create_v2p_matrix(*viewport)

        self.projections = None

        #GTK#
        #self.timer = wx.Timer(self)
        #self.Bind(wx.EVT_PAINT, self.__evt_paint_cb)
        #self.Bind(wx.EVT_SIZE, self.__evt_size_cb)
        #self.Bind(wx.EVT_MOUSEWHEEL, self.__evt_mousewheel_cb)
        #self.Bind(wx.EVT_MOTION, self.__evt_motion_cb)
        #self.Bind(wx.EVT_LEFT_DOWN, self.__evt_left_down_cb)
        #self.Bind(wx.EVT_LEFT_UP, self.__evt_left_up_cb)
        #self.Bind(wx.EVT_RIGHT_UP, self.__evt_right_up_cb)
        #self.Bind(wx.EVT_TIMER, self.__evt_timer_cb, self.timer)

        self.drawing_area.add_events(Gdk.EventMask.SCROLL_MASK |
                                     #Gdk.EventMask.POINTER_MOTION_MASK |
                                     Gdk.EventMask.BUTTON_PRESS_MASK |
                                     Gdk.EventMask.BUTTON_RELEASE_MASK |
                                     Gdk.EventMask.BUTTON1_MOTION_MASK |
                                     Gdk.EventMask.BUTTON3_MOTION_MASK)

        self.connect("size-allocate", self.__evt_size_cb)
        self.drawing_area.connect("draw", self.__evt_draw_cb)
        self.drawing_area.connect("scroll-event", self.__evt_mousewheel_cb)
        self.drawing_area.connect("motion-notify-event", self.__evt_motion_cb)
        self.drawing_area.connect("button-press-event", self.__evt_press_cb)
        self.drawing_area.connect("button-release-event", self.__evt_release_cb)
        #self.drawing_area.connect("button-release-event", self.__evt_release_cb)

        self.timer_id = None

    def show_emitter(self, show=True, alpha=0.25):
        _opts = self.sprites_opts.copy()
        if show:
            self.sprites_opts['alpha'] = alpha
            self.sprites_opts['glow'] = False
        else:
            self.sprites_opts = _opts.copy()

    def set_atoms(self, atoms, rescale=False, center=True):
        """
        Attach an Atoms object to the view.

        This will translate the model to the center of mass, move the model
        center to the center of screen and adjust the scale to the largest
        dimension of the model

        :param rescale: if True, the zoom is computed to view the atoms; if
                        False, a fixed zoom value is used
        """
        if atoms is None:
            self.light_mode = False
            self.atoms_center_of_mass = np.zeros(3)
            self.atoms_largest_dimension = 1.0
            self.atom_numbers = None
            self.atom_surfaces = None
            self.atoms_sprite = None
            self.projections = None
        else:
            # Set the light mode according to the number of atoms
            if len(atoms) > self.light_mode_threshold:
                self.light_mode = True
            else:
                self.light_mode = False

            # get the center of mass
            self.atoms_center_of_mass = atoms.get_center_of_mass()
            # get the largest dimension
            p = atoms.get_positions()
            self.atoms_largest_dimension = np.max(np.amax(p, axis=0) -
                                                  np.amin(p, axis=0))
            if self.atoms_largest_dimension == 0:
                self.atoms_largest_dimension = 1.0

        # make atoms a class attribute
        self.atoms = atoms
        # self.do_rescale = rescale
        # self.do_center = center
        self.update_camera(center=center, rescale=rescale)
        # create the textures
        self.create_atom_sprites()
        # finally update the view
        self.update_drawing()

    def rotate_atoms(self, dtheta, dphi):
        self.theta += dtheta
        self.phi += dphi

        tx, ty = (self.theta, self.phi)
        m_mat = np.zeros((4, 4))
        m_mat[0, 0] = m_mat[3, 3] = 1
        m_mat[1, 1] = m_mat[2, 2] = np.cos(np.radians(tx))
        m_mat[2, 1] = -np.sin(np.radians(tx))
        m_mat[1, 2] = np.sin(np.radians(tx))

        n_mat = np.zeros((4, 4))
        n_mat[1, 1] = n_mat[3, 3] = 1
        n_mat[0, 0] = n_mat[2, 2] = np.cos(np.radians(ty))
        n_mat[0, 2] = -np.sin(np.radians(ty))
        n_mat[2, 0] = np.sin(np.radians(ty))

        self.rotation_matrix = np.dot(m_mat, n_mat)
        self.update_model_matrix()
        self.scale_atoms(self.scale)

    def scale_atoms(self, factor):
        self.scale = factor
        self.scale_matrix[(0, 1, 2), (0, 1, 2)] = factor
        self.create_atom_sprites()
        self.update_projection_matrix()

    def translate_atoms(self, x, y):
        """
        sets the translation of the atoms
        """
        self.ox = x
        self.oy = y
        self.im_ox += self.last_mouse_move_x
        self.im_oy += self.last_mouse_move_y
        self.last_mouse_move_x = self.last_mouse_move_y = 0
        self.translation_matrix[-1, (0, 1)] = (x, y)
        self.update_projection_matrix()

    def select_atoms(self, x, y, w=None, h=None, append=False, toggle=False):
        print("Append=", append, "Toggle=", toggle)
        selection = np.array([])
        if w is None and h is None:
            # get the projections
            p = self.projections.copy()
            # translate to the event point
            p[:, :2] -= (x, y)
            # compute the norm and the radius for each projected atom
            norm = np.linalg.norm(p[:, :2], axis=1)
            radii = covalent_radii[p[:, 4].astype(int)] * self.scale
            # search where the norm is inside an atom
            i = np.where(norm < radii)
            # pick up the atom index of the one with the z min
            try:
                selection = np.array([int(p[i][np.argmin(p[i, 2]), 5])])
                # self.selection = np.array([selection])
            except:
                pass
        else:
            if w < 0:
                x += w
                w = abs(w)
            if h < 0:
                y += h
                h = abs(h)
            p = self.projections.copy()
            p = p[np.where(p[:, 0] > x)]
            p = p[np.where(p[:, 0] < x + w)]
            p = p[np.where(p[:, 1] > y)]
            p = p[np.where(p[:, 1] < y + h)]
            selection = p[:, -1].astype(int)

        if toggle:
            print(self.selection)
            # whether atoms in the current selection were previously selected
            i = np.in1d(self.selection, selection)
            print(i)
            self.selection = self.selection[np.invert(i)]

        if append:
            self.selection = np.append(self.selection, selection)
            self.selection = np.unique(self.selection)
        else:
            self.selection = selection

    def __evt_paint_cb(self, event):
        self.swap_buffers()

    def OLD__evt_size_cb(self, event):
        self.timer.Stop()
        self.timer.Start(self.refresh_delay)
        size = self.GetClientSize()
        self.back_buffer = cairo.ImageSurface(cairo.FORMAT_RGB24, *size)
        self.create_background_sprite(*size)
        self.update_drawing()

    def __evt_size_cb(self, widget, data):
        #self.timer.cancel()
        #self.timer.start()
        self._postpone_drawing()
        size = self.get_size()
        self.back_buffer = cairo.ImageSurface(cairo.FORMAT_RGB24, *size)
        self.create_background_sprite(*size)
        self.update_drawing()

    def __evt_timer_cb(self, event):
        self.update_drawing(light=False)
        self.timer.Stop()

    def __evt_left_down_cb(self, event):
        self.mx = event.GetX()
        self.my = event.GetY()
        self.capture_screen()
        if event.ControlDown():
            self.mode |= self.MODE_SELECTION
        if event.ShiftDown():
            self.mode |= self.MODE_SELECTION_APPEND
        if event.AltDown():
            self.mode |= self.MODE_SELECTION_TOGGLE

    def __evt_press_cb(self, widget, event):
        if event.button == 1:
            print("press_cb", event.state)
            self.mx = event.x
            self.my = event.y
            self.capture_screen()
            if event.state & Gdk.ModifierType.CONTROL_MASK:
                self.mode |= self.MODE_SELECTION
            if event.state & Gdk.ModifierType.SHIFT_MASK:
                self.mode |= self.MODE_SELECTION_APPEND
            if event.state & Gdk.ModifierType.MOD1_MASK:
                self.mode |= self.MODE_SELECTION_TOGGLE

    def __evt_release_cb(self, widget, event):
        if event.button not in (1,3):
            return

        if self.mode & self.MODE_SELECTION:
            self.mode ^= self.MODE_SELECTION
            # search for atoms in the selection box
            x, y = event.x, event.y
            w = h = None
            if self.mode & self.MODE_SELECTION_BOX:
                self.mode ^= self.MODE_SELECTION_BOX
                x, y, w, h = self.selection_box

            append = False
            if self.mode & self.MODE_SELECTION_APPEND:
                self.mode ^= self.MODE_SELECTION_APPEND
                append = True

            toggle = False
            if self.mode & self.MODE_SELECTION_TOGGLE:
                self.mode ^= self.MODE_SELECTION_TOGGLE
                toggle = True

            self.select_atoms(x, y, w, h, append=append, toggle=toggle)

        if self.mode == self.MODE_TRANSLATION:
            self.mode ^= self.MODE_TRANSLATION

        if self.mode & self.MODE_ROTATION:
            self.mode ^= self.MODE_ROTATION

        self.update_drawing(light=False)


    def __evt_left_up_cb(self, event):
        if self.mode & self.MODE_SELECTION:
            self.mode ^= self.MODE_SELECTION
            # search for atoms in the selection box
            x, y = event.GetPosition()
            w = h = None
            if self.mode & self.MODE_SELECTION_BOX:
                self.mode ^= self.MODE_SELECTION_BOX
                x, y, w, h = self.selection_box

            append = False
            if self.mode & self.MODE_SELECTION_APPEND:
                self.mode ^= self.MODE_SELECTION_APPEND
                append = True

            toggle = False
            if self.mode & self.MODE_SELECTION_TOGGLE:
                self.mode ^= self.MODE_SELECTION_TOGGLE
                toggle = True

            self.select_atoms(x, y, w, h, append=append, toggle=toggle)

        if self.mode == self.MODE_TRANSLATION:
            self.mode ^= self.MODE_TRANSLATION

        self.update_drawing(light=False)

    def __evt_right_up_cb(self, event):
        if self.mode & self.MODE_ROTATION:
            self.mode ^= self.MODE_ROTATION
        self.update_drawing(light=False)

    def OLD__evt_motion_cb(self, event):
        self.timer.Stop()
        self.timer.Start(self.refresh_delay)
        if event.LeftIsDown():
            mx, my = event.GetPosition()
            dx, dy = (mx - self.mx, my - self.my)
            # if event.ControlDown():
            if self.mode & self.MODE_SELECTION:
                self.mode |= self.MODE_SELECTION_BOX
                # if event.ShiftDown():
                #    self.mode |= self.MODE_SELECTION_APPEND
                self.selection_box = [self.mx, self.my, dx, dy]
            else:
                self.mode = self.MODE_TRANSLATION
                self.mx, self.my = (mx, my)
                self.last_mouse_move_x = int(dx)
                self.last_mouse_move_y = int(dy)
                self.ox = int(self.ox + dx)
                self.oy = int(self.oy + dy)
                self.translate_atoms(self.ox, self.oy)
            self.update_drawing()
        elif event.RightIsDown():
            self.mode = self.MODE_ROTATION
            theta = 2. * (float(self.scale0) / self.scale)
            theta = max(1., theta)
            mx, my = event.GetPosition()
            dx, dy = (mx - self.mx, my - self.my)
            self.mx, self.my = (mx, my)

            tx = theta * np.sign(dy)
            ty = theta * np.sign(dx)
            self.rotate_atoms(tx, ty)

            self.update_drawing()

    def __evt_motion_cb(self, widget, event):
        #self.timer.cancel()
        #self.timer.start()
        self._postpone_drawing()
        self._update_flag = True
        GLib.timeout_add(self.refresh_delay, self.update_drawing, False)
        if event.state & Gdk.ModifierType.BUTTON1_MASK:
            print("motion cb...", bool(event.state &
                  Gdk.ModifierType.BUTTON1_MASK))
            mx, my = event.x, event.y
            dx, dy = (mx - self.mx, my - self.my)
            # if event.ControlDown():
            if self.mode & self.MODE_SELECTION:
                self.mode |= self.MODE_SELECTION_BOX
                # if event.ShiftDown():
                #    self.mode |= self.MODE_SELECTION_APPEND
                self.selection_box = [self.mx, self.my, dx, dy]
            else:
                self.mode = self.MODE_TRANSLATION
                self.mx, self.my = (mx, my)
                self.last_mouse_move_x = int(dx)
                self.last_mouse_move_y = int(dy)
                self.ox = int(self.ox + dx)
                self.oy = int(self.oy + dy)
                self.translate_atoms(self.ox, self.oy)
            self.update_drawing()
        #elif event.RightIsDown():
        elif event.state & Gdk.ModifierType.BUTTON3_MASK:
            self.mode = self.MODE_ROTATION
            theta = 2. * (float(self.scale0) / self.scale)
            theta = max(1., theta)
            mx, my = event.x, event.y
            dx, dy = (mx - self.mx, my - self.my)
            self.mx, self.my = (mx, my)

            tx = theta * np.sign(dy)
            ty = theta * np.sign(dx)
            self.rotate_atoms(tx, ty)

            self.update_drawing()

    def OLD__evt_mousewheel_cb(self, event):
        if wx.GetKeyState(wx.WXK_CONTROL):
          alpha = self.sprites_opts['alpha']
          rot = event.GetWheelRotation()
          if rot > 0:
            alpha *= 1.2
            alpha = min(1, alpha)
          elif rot < 0:
            alpha /= 1.2
            alpha = max(0, alpha)
          self.sprites_opts['alpha'] = alpha
          self.create_atom_sprites()
          self.update_drawing()
        else:
          rot = event.GetWheelRotation()
          self.timer.Stop()
          self.timer.Start(self.refresh_delay)
          if rot > 0:
              factor = self.scale * 1.1
              im_factor = 1 * 1.1
          elif rot < 0:
              factor = self.scale / 1.1
              im_factor = 1 / 1.1
          self.im_factor = im_factor
          self.scale_atoms(factor)
          self.update_drawing()

    def __evt_mousewheel_cb(self, widget, event):
        #if wx.GetKeyState(wx.WXK_CONTROL):
        if event.state & Gdk.ModifierType.CONTROL_MASK:
          alpha = self.sprites_opts['alpha']
          #rot = event.GetWheelRotation()
          rot = event.direction
          if rot == Gdk.ScrollDirection.UP:
            alpha *= 1.2
            alpha = min(1, alpha)
          elif rot == Gdk.ScrollDirection.DOWN:
            alpha /= 1.2
            alpha = max(0, alpha)
          self.sprites_opts['alpha'] = alpha
          self.create_atom_sprites()
          self.update_drawing()
        else:
          #rot = event.GetWheelRotation()
          rot = event.direction
          #self.timer.cancel()
          #self.timer.start()
          self._postpone_drawing()
          if rot == Gdk.ScrollDirection.UP:
              factor = self.scale * 1.1
              im_factor = 1 * 1.1
          elif rot == Gdk.ScrollDirection.DOWN:
              factor = self.scale / 1.1
              im_factor = 1 / 1.1
          self.im_factor = im_factor
          self.scale_atoms(factor)
          self.update_drawing()


    def target_fn(self):
        print("Timer function!!")
        return True

    def _postpone_drawing(self):
        if self.timer_id:
            GLib.source_remove(self.timer_id)
        self.timer_id = GLib.timeout_add(self.refresh_delay,
                                         self.update_drawing, False)


    def __postpone_drawing(self):
        return
        try:
            self.timer.cancel()
        except:
            #self.timer = Timer(self.refresh_delay, self.update_drawing,
            #                   kwargs={'light_mode': False})
            self.timer = Timer(self.refresh_delay, self.target_fn)
            self.timer.start()
            print("Timer created")

    def capture_screen(self):
        # get size of screen
        w, h = self.get_size()
        # create a cairo surface and context
        surface = cairo.ImageSurface(cairo.FORMAT_RGB24, w, h)
        ctx = cairo.Context(surface)
        # trick here: blit the last back_buffer onto the newly created surface
        ctx.set_source_surface(self.back_buffer)
        ctx.paint()
        # store it as an attribute
        self.screenshot = surface

    def create_atom_sprites(self):
        """
        This function creates a list of cairo surfaces for each kind
        of atoms
        """

        # Get out if there are no atoms
        if not self.atoms:
            return

        # First get an array of all atoms numbers
        atom_numbers = np.unique(self.atoms.numbers)

        # Now, for each kind of atoms create a surface in memory
        atom_surfaces = np.empty((2, len(atom_numbers)), dtype=object)
        self.__outer_margin = 0
        def create_surface(atom_number, alpha=1, glow=True):
            #global margin
            # get the radius, and the color
            radius = int(covalent_radii[atom_number] * 1. * self.scale)
            r, g, b = jmol_colors[atom_number]
            # actually create the surface
            size = 2 * radius #+ 4
            self.__outer_margin = np.maximum(self.__outer_margin, size / 2.)
            surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, size, size)
            # draw the ball
            ctx = cairo.Context(surface)
            # ctx.set_antialias(cairo.ANTIALIAS_NONE)
            ctx.set_line_width(1.)
            ctx.set_source_rgba(r, g, b, alpha)
            ctx.arc(radius, radius, radius - 0.5, 0, 2 * np.pi)
            ctx.fill_preserve()
            if glow:
                gradient = cairo.RadialGradient(radius, radius, radius / 2,
                                                radius, radius, radius)
                gradient.add_color_stop_rgba(0., 1., 1., 1., .5*alpha)
                gradient.add_color_stop_rgba(0.5, 1., 1., 1., 0)
                gradient.add_color_stop_rgba(1., 1., 1., 1., 0.)
                ctx.set_source(gradient)
                ctx.fill_preserve()
            ctx.set_source_rgba(0., 0., 0., alpha)
            ctx.stroke()

            # Create the overlay for selection
            overlay = cairo.ImageSurface(cairo.FORMAT_ARGB32, size, size)
            # draw the circle
            ctx = cairo.Context(overlay)
            ctx.set_source_surface(surface)
            ctx.paint()
            ctx.set_line_width(2.)
            ctx.set_source_rgb(1 - r, 1 - g, 1 - b)
            ctx.arc(radius, radius, radius - 2., 0, 2 * np.pi)
            ctx.stroke()

            return surface, overlay

        for i, a in enumerate(atom_numbers):
            surface, overlay = create_surface(a,
                                              alpha=self.sprites_opts['alpha'],
                                              glow=self.sprites_opts['glow'])
            atom_surfaces[0, i] = surface
            atom_surfaces[1, i] = overlay
            """
            # get the radius, and the color
            radius = int(covalent_radii[a] * 1. * self.scale)
            # b, g, r = jmol_colors[a]
            r, g, b = jmol_colors[a]
            # actually create the surface
            size = 2 * radius + 4
            margin = np.maximum(margin, size / 2.)
            surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, size, size)
            # draw the ball
            ctx = cairo.Context(surface)
            # ctx.set_antialias(cairo.ANTIALIAS_NONE)
            ctx.set_line_width(1.)
            ctx.set_source_rgba(r, g, b, self.sprites_opts['alpha'])
            ctx.arc(radius, radius, radius - 0.5, 0, 2 * np.pi)
            ctx.fill_preserve()
            if self.sprites_opts['glow']:
                gradient = cairo.RadialGradient(radius, radius, radius / 2,
                                                radius, radius, radius)
                gradient.add_color_stop_rgba(0., 1., 1., 1., .5)
                gradient.add_color_stop_rgba(0.5, 1., 1., 1., 0)
                gradient.add_color_stop_rgba(1., 1., 1., 1., 0.)
                ctx.set_source(gradient)
                ctx.fill_preserve()
            ctx.set_source_rgba(0., 0., 0., self.sprites_opts['alpha'])
            ctx.stroke()
            # store it
            atom_surfaces[0, i] = surface

            # Create the overlay for selection
            overlay = cairo.ImageSurface(cairo.FORMAT_ARGB32, size, size)
            # draw the circle
            ctx = cairo.Context(overlay)
            ctx.set_source_surface(surface)
            ctx.paint()
            ctx.set_line_width(2.)
            ctx.set_source_rgb(1 - r, 1 - g, 1 - b)
            ctx.arc(radius, radius, radius - 2., 0, 2 * np.pi)
            ctx.stroke()
            atom_surfaces[1, i] = overlay
        """

        self.atom_numbers = atom_numbers
        self.atom_surfaces = atom_surfaces
        try:
            absorber_number = self.atoms[self.atoms.info['absorber']].number
            self.absorber_surface = create_surface(absorber_number, alpha=1,
                                                   glow=True)
        except:
            self.atoms.info['absorber'] = -1
        self.__outer_margin *= 1.1

    def create_background_sprite(self, w, h):
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, w, h)
        ctx = cairo.Context(surface)

        if True:  # pylint: disable=using-constant-test
            g = cairo.LinearGradient(0, 0, 0, h)
            g.add_color_stop_rgba(0.0, 1.0, 1.0, 1.0, 1.0)
            g.add_color_stop_rgba(0.7, 1.0, 1.0, 1.0, 1.0)
            g.add_color_stop_rgba(1.0, 0.5, 0.5, 0.5, 1.0)
            ctx.set_source(g)
            ctx.rectangle(0, 0, w, h)
            ctx.fill()

        g = cairo.LinearGradient(0, 0, 0, h)
        #g.add_color_stop_rgba(0., 1., 1., 1., 1.)
        #g.add_color_stop_rgba(2 / 3., 0.5, 0.5, 0.5, 1)
        #g.add_color_stop_rgba(0.0, 1.0, 1.0, 1.0, 1.0)
        #g.add_color_stop_rgba(0.9, 0.8, 0.8, 0.8, 1.0)
        #g.add_color_stop_rgba(1.0, 0.2, 0.2, 0.2, 1.0)
        #ctx.set_source(g)
        ctx.set_source_rgb(1, 1, 1)
        ctx.rectangle(0, 0, w, h)
        ctx.fill()

        ctx.save()

        if False:
            ctx.set_source_rgb(0.8, 0.8, 0.8)
            rect = (0, 2 * h / 3, w, h / 3)
            ctx.rectangle(*rect)
            ctx.clip()
            ctx.paint()
            ctx.set_line_width(1.)
            for i in np.arange(0, 2 * np.pi, np.pi / 30):
                ctx.move_to(w / 2, 2 * h / 3)
                x1 = w * np.cos(i)
                y1 = w * np.sin(i)
                ctx.rel_line_to(x1, y1)
            for i in np.arange(2 * h / 3, h, 10):
                ctx.move_to(0, i)
                ctx.line_to(w, i)

            ctx.set_source_rgb(0.7, 0.7, 0.7)
            ctx.stroke()
            ctx.restore()

        self.background_sprite = surface

    @classmethod
    def create_v2p_matrix(cls, left, right, bottom, top, near, far):
        """
        creates the matrix that transforms coordinates from view space (space
        defined by the bounding box passed as argument) to projection space

        this transformation is a scale and offset that maps [left; right],
        [bottom; top], [near; far] to [-1;1], [-1;1], [0;1]
        """
        v2p_matrix = np.eye(4) * -1
        v2p_matrix[0, 0] = 2. / (right - left)
        v2p_matrix[1, 1] = 2. / (right - left)
        v2p_matrix[2, 2] = 1. / (near - far)
        v2p_matrix[3, 0] = (left + right) / (left - right)
        v2p_matrix[3, 1] = (top + bottom) / (bottom - top)
        v2p_matrix[3, 2] = near / (near - far)
        return v2p_matrix

    def update_projection_matrix(self):
        m_matrix = np.dot(self.v2p_matrix, self.scale_matrix)
        m_matrix = np.dot(m_matrix, self.translation_matrix)
        self.projection_matrix = m_matrix

    def update_model_matrix(self):
        m_matrix = np.dot(self.m2w_matrix, self.rotation_matrix)
        self.model_matrix = m_matrix

    def get_projections(self, points, save=False):
        m_matrix = np.dot(self.model_matrix, self.projection_matrix)

        p = points[:, :4]
        v = np.dot(p, m_matrix)
        v = v / v[:, -1, None]

        # add the other columns
        v = np.c_[v, points[:, 4:]]
        # and sort by Z
        # v = v[v[:,2].argsort()[::-1]]

        if save:
            self.projections = v
        return v

    def filter_projections(self, projections, w, h):
        try:
            # filtering
            margin = self.__outer_margin
            projections = projections[projections[:, 0] >= -1 * margin, :]
            projections = projections[projections[:, 0] <= w + margin, :]
            projections = projections[projections[:, 1] >= -1 * margin, :]
            projections = projections[projections[:, 1] <= h + margin, :]
            return projections

        except:
            pass

    def render_background(self, ctx):
        surface = self.background_sprite
        ctx.set_source_surface(surface, 0, 0)
        ctx.paint()

    def render_scalebar(self, ctx):
        x, y, w, h = ctx.clip_extents()  # @UnusedVariable
        scalebar_bb_width = 200
        scalebar_bb_height = 20
        ctx.set_source_rgba(0., 0., 0., 0.7)
        ctx.rectangle(x + w - scalebar_bb_width - 6, h - scalebar_bb_height - 6,
                      scalebar_bb_width, scalebar_bb_height)
        ctx.fill()

        ctx.set_source_rgb(1, 1, 1)
        ctx.rectangle(x + w - scalebar_bb_width, h - scalebar_bb_height, 100,
                      scalebar_bb_height - 12)
        ctx.fill()

        ctx.move_to(x + w - scalebar_bb_width / 2 + 6, h - 9)
        ctx.set_source_rgb(1, 1, 1)
        ctx.set_font_size(16)
        ctx.show_text("%.2f \xc5" % (100. / self.scale))

    def render_axes(self, ctx):
        _, _, w, h = ctx.clip_extents()  # @UnusedVariable
        m_matrix = np.dot(self.rotation_matrix, self.v2p_matrix)

        d = 20
        offset = 12

        points = np.array([[d, 0, 0, 1, 0],
                           [0, d, 0, 1, 1],
                           [0, 0, d, 1, 2]])

        # project onto viewport
        projections = np.dot(points[:, :4], m_matrix)
        projections /= projections[:, -1, None]
        # add the other columns
        projections = np.c_[projections, points[:, 4:]]
        # and sort by Z
        projections = projections[projections[:,2].argsort()[::-1]]

        red = (1, 0, 0)
        green = (0, 0.7, 0)
        blue = (0, 0, 1)
        colors = (red, green, blue)
        labels = ('X', 'Y', 'Z')

        # draw a white circle
        so = np.array([d + offset, h - d - offset, 0])
        ctx.move_to(d + offset, h - d - offset)
        ctx.arc(d + offset, h - d - offset, d + offset - 2, 0, 2 * np.pi)
        #ctx.set_source_rgb(1, 1, 1)
        ctx.set_source_rgba(0., 0., 0., 0.7)
        ctx.set_line_width(1)
        #ctx.stroke_preserve()
        ctx.set_source_rgba(0.95, 0.95, 0.95, 1)
        ctx.fill()

        for x, y, z, w, n in projections:
            n = int(n)
            ctx.move_to(*so[:2])
            ctx.rel_line_to(x, y)
            ctx.set_source_rgb(*colors[n])
            ctx.set_line_width(2)
            ctx.set_font_size(10)
            ctx.show_text(labels[n])
            ctx.stroke()


    def render_atoms(self, ctx):
        try:
            atoms = self.atoms
        except:
            return

        # create a points matrix with homogeneous coordinates
        # x,y,z,w and the atom number and index
        points = atoms.get_positions()
        points = np.c_[points, np.ones(len(points))]
        points = np.c_[points, atoms.numbers]
        points = np.c_[points, np.arange(len(atoms))]

        # Get the points array projected to the screen display
        projections = self.get_projections(points, save=True)
        # print('render_atoms : len(projections) = %d' % len(projections))
        # print('render_atoms : projections = %s' % str(projections))

        # Reduce the number of atoms to be drawn if outside the viewport
        w, h = self.get_size()
        # print('render_atoms : w = %f, h = %f' % (w, h))
        projections = self.filter_projections(projections, w, h)
        # self.projections = projections

        if len(projections) == 0:
            print('render_atoms : no atom is visible')
            return  # no atom is visible from the camera

        # sort by z
        projections = projections[projections[:, 2].argsort()[::-1]]

        margin = self.__outer_margin
        xmin, ymin = np.min(projections[:, :2], axis=0)
        xmax, ymax = np.max(projections[:, :2], axis=0)
        sw = xmax - xmin + 2 * margin
        sh = ymax - ymin + 2 * margin
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, int(sw), int(sh))
        surface_ctx = cairo.Context(surface)

        # set the local references
        set_source_surface = surface_ctx.set_source_surface
        paint = surface_ctx.paint
        selection = self.selection
        atom_numbers = self.atom_numbers
        atom_surfaces = self.atom_surfaces

        self.im_ox = int(xmin - margin)
        self.im_oy = int(ymin - margin)
        surface_ctx.translate(-self.im_ox, -self.im_oy)
        # surface_ctx.set_source_rgb(1,1,0)
        # surface_ctx.set_line_width(2.)
        # surface_ctx.rectangle(self.im_ox, self.im_oy, int(sw), int(sh))
        # surface_ctx.stroke()
        for i, p in enumerate(projections):  # @UnusedVariable
            x, y, z, w, n, index = p  # @UnusedVariable
            # load the surface
            if index in selection:
                if index == self.atoms.info['absorber']:
                    sprite = self.absorber_surface[1]
                else:
                    sprite = atom_surfaces[1, np.where(atom_numbers == n)[0][0]]
            else:
                if index == self.atoms.info['absorber']:
                    sprite = self.absorber_surface[0]
                else:
                    sprite = atom_surfaces[0, np.where(atom_numbers == n)[0][0]]

            sx = x - sprite.get_width() / 2.
            sy = y - sprite.get_height() / 2.
            set_source_surface(sprite, int(sx), int(sy))
            paint()
            if False:  # pylint: disable=using-constant-test
                ctx.set_source_rgb(0., 0., 0.)
                r = sprite.get_width() / 2
                ctx.arc(x, y, r, 0, 2 * np.pi)
                ctx.stroke_preserve()
                ctx.set_source_rgb(0.8, 0.8, 0.8)
                ctx.fill()
                ctx.move_to(x, y)
                ctx.set_source_rgb(0, 0, 0)
                ctx.set_font_size(14)
                ctx.show_text("%d" % index)

        # save the rendering
        self.atoms_sprite = surface
        ctx.set_source_surface(surface, self.im_ox, self.im_oy)
        ctx.paint()

    def render_selection_box(self, ctx):
        r, g, b = self.colors['selection_box']
        ctx.set_source_surface(self.screenshot, 0, 0)
        ctx.paint()
        ctx.set_source_rgba(r, g, b, 0.3)
        ctx.rectangle(*self.selection_box)
        ctx.fill_preserve()
        ctx.set_source_rgb(r, g, b)
        ctx.stroke()

    def render_boundingbox(self, ctx):
        # print('render_boundingbox : start')
        try:
            atoms = self.atoms
        except:
            return

        # create a points matrix with homogeneous coordinates
        # x,y,z,w for atoms extrema
        points = atoms.get_positions()
        margin = self.__outer_margin / float(self.scale)
        xmin, ymin, zmin = np.min(points, axis=0) - margin
        xmax, ymax, zmax = np.max(points, axis=0) + margin
        points = np.array([
            [xmax, ymax, zmax, 1],
            [xmax, ymin, zmax, 1],
            [xmin, ymax, zmax, 1],
            [xmin, ymin, zmax, 1],
            [xmax, ymax, zmin, 1],
            [xmax, ymin, zmin, 1],
            [xmin, ymax, zmin, 1],
            [xmin, ymin, zmin, 1]])

        # Get the points array projected to the screen display
        projections = self.get_projections(points)
        x0, y0 = (np.min(projections[:, :2], axis=0)).astype(int)
        x1, y1 = (np.max(projections[:, :2], axis=0)).astype(int)

        # Declare the 6 faces with their vertex index
        # the order of the numbers define if the normal plane points outward or not
        faces = np.array([
            [6, 7, 5, 4],
            [2, 3, 7, 6],
            [3, 1, 5, 7],
            [2, 6, 4, 0],
            [0, 4, 5, 1],
            [2, 0, 1, 3]])

        # kind of backface culling
        ind = []
        for i, f in enumerate(faces):
            # Get 2 vectors of the plane
            v1 = projections[f[1], :3] - projections[f[0], :3]
            v2 = projections[f[3], :3] - projections[f[0], :3]
            # cross multiply them to get the normal
            n = np.cross(v2, v1)
            # If the normal z coordinate is <0, the plane is not visible, so,
            # draw it first, otherwise draw it last
            if n[-1] > 0:
                ind.append(i)
            else:
                ind.insert(0, i)
        # faces are now re-ordered
        faces = faces[ind]

        # plane transparency and color
        color_plane = self.colors['boulding_box_fill']
        color_line = self.colors['boulding_box_line']

        ctx.save()
        #ctx.set_source_rgb(1, 0, 0)
        #ctx.rectangle(x0, y0, x1 - x0, y1 - y0)
        #ctx.stroke()
        ctx.set_fill_rule(cairo.FILL_RULE_EVEN_ODD)
        ctx.set_line_join(cairo.LINE_JOIN_ROUND)
        ctx.set_line_width(1.)
        ctx.set_dash([8., 8.])
        for i, p in enumerate(faces):
            # remove dash for front faces
            if i > 2:
                ctx.set_dash([])
            if i == 3 and not(self.mode & self.MODE_ROTATION):
                try:
                    f = self.scale / self.scale0
                    # sprite_w = self.atoms_sprite.get_width()
                    # sprite_h = self.atoms_sprite.get_height()
                    # m = self.__outer_margin
                    # sw = float(sprite_w) / (x1 - x0)
                    # sh = float(sprite_h) / (y1 - y0)
                    ctx.save()
                    # ctx.translate((x0 + x1 - sprite_w)/2. ,
                    #                            (y0 + y1 - sprite_h)/2.)
                    ctx.set_source_surface(self.atoms_sprite, self.im_ox, self.im_oy)
                    ctx.paint()
                    ctx.restore()
                except:
                    pass
            ctx.set_source_rgba(*color_plane)
            # get the projected points
            p0, p1, p2, p3 = projections[p, :2]
            # move to the first
            ctx.move_to(*p0)
            # line to the others
            list(map(lambda _: ctx.line_to(*_), (p1, p2, p3)))
            # close the polygon
            ctx.close_path()
            # fill it and stroke the path
            ctx.fill_preserve()
            ctx.set_source_rgba(*color_line)
            ctx.stroke()
        ctx.restore()

    def update_camera(self, rescale=False, center=False):
        # update the scale
        w, h = self.get_size()
        l = min(w, h)
        # print('set_atoms : w=%d h = %d l = %f' % (w, h, l))
        self.scale0 = int(.5 * l / self.atoms_largest_dimension) + 1

        # move the model to the center of mass
        t_matrix = np.eye(4)
        t_matrix[-1, :-1] = -1 * self.atoms_center_of_mass
        self.m2w_matrix = t_matrix  # self.m2w_matrix.dot(t_matrix)
        self.update_model_matrix()

        if rescale:
            assert self.scale0 > 0.0
            self.scale = self.scale0

        #print "scale  = ", self.scale, "scale0 = ", self.scale0
        #self.scale_atoms(1.)

        if center:
            self.translate_atoms(w / 2, h / 2)

    def update_drawing(self, light=None):
        # print('update_drawing : light=%s' % str(light))

        try:
            ctx = cairo.Context(self.back_buffer)
        except:
            #self.scale = self.scale0 = 10
            return

        light_mode = self.light_mode if light is None else light

        if self.mode & self.MODE_SELECTION:
            self.render_selection_box(ctx)
        else:
            self.render_background(ctx)
            if self.atoms:
                if light_mode:
                    self.render_boundingbox(ctx)
                else:
                    self.render_atoms(ctx)
                self.render_scalebar(ctx)
            self.render_axes(ctx)

        #Gtk#
        #self.Refresh(eraseBackground=False)
        allocation = self.drawing_area.get_allocation()
        self.drawing_area.queue_draw_area(allocation.x, allocation.y,
                                          allocation.width, allocation.height)
        return False

    def swap_buffers(self):
        if self.back_buffer:
            back_buffer = self.back_buffer
            # w = back_buffer.get_width()
            # h = back_buffer.get_height()

            bitmap = wx.lib.wxcairo.BitmapFromImageSurface(back_buffer)
            dc = wx.PaintDC(self)
            dc.DrawBitmap(bitmap, 0, 0)

    def __evt_draw_cb(self, da, ctx):
        if self.back_buffer:
            print("Draw called.", ctx)
            back_buffer = self.back_buffer
            ctx.set_source_surface(back_buffer)
            ctx.paint()


if __name__ == "__main__":
  from ase.build import bulk
  from ase.visualize import view
  import wx

  MgO = bulk('MgO', crystalstructure='rocksalt', a=4.21, cubic=True)
  MgO = MgO.repeat((10,10,10))
  #view(MgO)

  #app = wx.App(False)
  #app.SetAppName('MsSpec Data Viewer')
  #win = wx.Frame(None, size=wx.Size(480, 340))
  #cv = ClusterViewer(win, size=wx.Size(480, 340))
  #cv.set_atoms(MgO, rescale=True, center=True)
  #win.Show()
  #app.MainLoop()

  cv = ClusterViewer()
  cv.connect("destroy", Gtk.main_quit)
  cv.set_atoms(MgO, rescale=True, center=True)
  cv.show_all()
  Gtk.main()

