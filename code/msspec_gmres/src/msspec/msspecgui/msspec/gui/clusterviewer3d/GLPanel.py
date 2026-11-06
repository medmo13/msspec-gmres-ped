import wx
from wx import glcanvas
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
import sys

from ase import Atom, Atoms
from ase.data import covalent_radii
from ase.data.colors import jmol_colors

from ase.build import molecule

import numpy as np

class GLPanel(wx.Panel):
    def __init__(self, *args, **kwargs):
        #style = wx.DEFAULT_FRAME_STYLE | wx.NO_FULL_REPAINT_ON_RESIZE
        wx.Panel.__init__(self, *args, **kwargs)

        self.GLinitialized = False
        attribList = (glcanvas.WX_GL_RGBA, # RGBA
                      glcanvas.WX_GL_DOUBLEBUFFER, # Double Buffered
                      glcanvas.WX_GL_DEPTH_SIZE, 24) # 24 bit
        # Create the canvas
        self.canvas = glcanvas.GLCanvas(self, attribList=attribList)
        self.glcontext = glcanvas.GLContext(self.canvas)
        # Set the event handlers.
        self.canvas.Bind(wx.EVT_SIZE, self.onSizeEvent)
        self.canvas.Bind(wx.EVT_PAINT, self.onPaintEvent)

        szr = wx.BoxSizer(wx.HORIZONTAL)
        szr.Add(self.canvas, 1, wx.EXPAND|wx.ALL, 0)
        self.SetSizer(szr)

        self.Bind(wx.EVT_MOUSEWHEEL, self.onMouseWheelEvent)
        self.canvas.Bind(wx.EVT_MOTION, self.onMotionEvent)
        self.canvas.Bind(wx.EVT_LEFT_DOWN, self.onLeftDownEvent)
        #self.Bind(wx.EVT_LEFT_UP, self.__evt_left_up_cb)
        #self.Bind(wx.EVT_RIGHT_UP, self.__evt_right_up_cb)
        self.canvas.Bind(wx.EVT_RIGHT_DOWN, self.onRightDownEvent)

        self.atoms = None
        self.scale = 1.
        self.theta = self.phi = 0.
        self.tx = self.ty = 0
        self.tx0 = self.ty0 = 0.
        self.znear = 1.
        self.zfar = -1.
        self.coords = [0., 0., 0.]
        self.translation = [0., 0., 0]

    def set_atoms(self, atoms):
        self.atoms = atoms
        self.com = atoms.get_center_of_mass()
        self.sorted_indices = np.argsort(atoms.numbers)
        self.theta = 45
        self.phi = -45


    def onMouseWheelEvent(self, event):
        rotation = event.GetWheelRotation()
        if rotation > 0:
            self.scale *= 1.1
        else:
            self.scale /= 1.1
        #self._gl_scale()
        print(self.scale)
        w, h = self.canvas.GetClientSize()
        self._gl_init_view(w, h)
        self._gl_draw()

    def onRightDownEvent(self, event):
        self.mx0, self.my0 = event.GetPosition()

    def onLeftDownEvent(self, event):
        mvmatrix = glGetDoublev(GL_MODELVIEW_MATRIX)
        projmatrix = glGetDoublev(GL_PROJECTION_MATRIX)
        viewport = glGetIntegerv(GL_VIEWPORT)
        winx, winy = event.GetPosition()
        winy = viewport[3] - winy
        winz = glReadPixels(winx, winy, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT)[0][0]
        coords = gluUnProject(winx, winy, winz, mvmatrix, projmatrix, viewport)
        self.coords0 = coords
        print("coords0 = ", coords)


    def onMotionEvent(self, event):
        w, h = self.canvas.GetClientSize()
        if event.RightIsDown():
            mx, my = event.GetPosition()
            dy = my - self.my0
            dx = mx - self.mx0
            self.my0 = my
            self.mx0 = mx
            self.theta += dx
            self.phi += dy
        elif event.LeftIsDown():
            mvmatrix = glGetDoublev(GL_MODELVIEW_MATRIX)
            projmatrix = glGetDoublev(GL_PROJECTION_MATRIX)
            viewport = glGetIntegerv(GL_VIEWPORT)

            # get x0, y0 and z0
            #x0, y0, z0 = gluUnProject(self.tx0, self.ty0, self.tz0, mvmatrix, projmatrix, viewport)

            winx, winy = event.GetPosition()
            winy = viewport[3] - winy
            winz = glReadPixels(winx, winy, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT)[0][0]
            x, y, z = gluUnProject(winx, winy, winz, mvmatrix, projmatrix, viewport)

            #translation = [x - x0, y - y0, 0]
            translation = [x - self.coords0[0], y - self.coords0[1], 0]
            #print x, y, z#translation
            if winx < w and winx > 0 and winy < h and winy > 0:
                print('here')
                glTranslatef(*translation)
                #self.translation = translation
            print(translation)

        self._gl_draw()

    def onSizeEvent(self, event):
        """Process the resize event."""
        width, height = self.canvas.GetClientSize()
        self._gl_resize(width, height)
        self.canvas.Refresh(False)
        event.Skip()

    def onPaintEvent(self, event):
        """Process the drawing event."""
        if not self.GLinitialized:
            w, h = self.canvas.GetClientSize()
            self._gl_init(w, h) 

        self._gl_draw() 
        event.Skip()

    # GLFrame OpenGL Event Handlers
    def _gl_init(self, width, height):
        """Initialize OpenGL for use in the window."""
        self.canvas.SetCurrent(self.glcontext)

        self.quadric = gluNewQuadric()
        gluQuadricNormals(self.quadric, GLU_SMOOTH)

        glClearColor(1., 1., 1., 1.)
        glClearDepth(1.)
        glDepthFunc(GL_LESS)
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_MULTISAMPLE)
        glShadeModel(GL_SMOOTH)

        #glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)



        #glLightfv(GL_LIGHT0, GL_AMBIENT,  (0.5, 0.5, 0.5, 1.0))
        a = 0.4
        d = 0.4
        s = 1.
        glLightfv(GL_LIGHT0, GL_AMBIENT,  (a, a, a, 1.0))
        glLightfv(GL_LIGHT0, GL_DIFFUSE,  (d, d, d, 1.0))
        glLightfv(GL_LIGHT0, GL_SPECULAR, (s, s, s, 1.0))
        glLightfv(GL_LIGHT0, GL_POSITION, (0.0, 100.0, 10., 1.0))
        glEnable(GL_LIGHT0)

        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE)
        #glLightfv(GL_LIGHT1, GL_AMBIENT, (0.3, 0.3, 0.3, 1.0))
        #glLightfv(GL_LIGHT1, GL_DIFFUSE, (.2, .2, .2, 1.0))
        #glLightfv(GL_LIGHT1, GL_POSITION, (0.0, -10.0, 0., 1.0))
        #glEnable(GL_LIGHT1)

        glEnable(GL_LIGHTING)
        glEnable(GL_COLOR_MATERIAL)

        self._gl_init_view(width, height)
        self._gl_resize(width, height)
        self.GLinitialized = True

    def _gl_init_view(self, width, height):
        glViewport(0, 0, width, height)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(45., float(width) / float(height), self.znear, self.zfar)
        # find the largest cluster dimension
        max_d = np.max(np.linalg.norm(self.atoms.get_positions(), axis=1))
        # set the scale accordingly
        gluLookAt(0, 0, 5*max_d*self.scale, 0, 0, 0, 0, 1, 0)
        #gluLookAt(0, 0, self.scale, 0, 0, 0, 0, 1, 0)
        glMatrixMode(GL_MODELVIEW)


    def _gl_resize(self, width, height):
        """Reshape the OpenGL viewport based on the dimensions of the window."""
        if height == 0:
            height = 1
        self._gl_init_view(width, height)

    def _gl_draw(self, *args, **kwargs):
        "Draw the window."
        #glClear(GL_COLOR_BUFFER_BIT)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        # Drawing an example triangle in the middle of the screen

        #glBegin(GL_TRIANGLES)
        #glColor(0, 0, 0)
        #color = [1.0, 0.0, 0.0, 1.0]
        #glMaterialfv(GL_FRONT, GL_DIFFUSE, color)
        #glVertex3f(-.25, -.25, 0.)
        #glVertex3f(.25, -.25, 0.)
        #glVertex3f(0., .25, 0.)
        ##glVertex(-.25, -.25)
        ##glVertex(.25, -.25)
        ##glVertex(0., .25)
        #glEnd()

        def draw_circle(radius, n):
            theta = 2*np.pi/float(n)
            c = np.cos(theta)
            s = np.sin(theta)

            x = radius
            y = 0

            glBegin(GL_LINE_LOOP)
            for i in range(n):
                glVertex2f(x , y )
                t = x
                x = c * x - s * y
                y = s * t + c * y
            glEnd()

        def render_atom(atom):
            glPushMatrix()

            glRotate(self.theta, 0, 1, 0)
            glRotate(self.phi, 1, 0, 0)
            glTranslatef(-self.com[0], -self.com[1], -self.com[2])

            glTranslatef(*atom.position)
            #r, g, b = jmol_colors[atom.number]
            radius = covalent_radii[atom.number]
            #color = [r, g, b, 1.0]
            #glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, (.1, .1, .1, 1))
            #glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, (0, 0, 0, 1))
            #glColor3f(r, g, b)
            resolution = 16
            gluSphere(self.quadric, radius, resolution, resolution)
            #draw_circle(radius, 32)
            glPopMatrix()

        def set_material(number):
            r, g, b = jmol_colors[number]
            color = [r, g, b, 1.0]
            s = 0.7
            glMaterialfv(GL_FRONT, GL_SPECULAR, (s, s, s, 1))
            glMaterialfv(GL_FRONT, GL_EMISSION, (0., 0., 0., 1))
            glMateriali(GL_FRONT, GL_SHININESS, 30)
            glColor4f(r, g, b, 0.5)

        last_number = -1
        #glPushMatrix()
        #glLoadIdentity()
        #glTranslatef(*self.translation)
        #glRotate(self.theta, 0, 1, 0)
        #glRotate(self.phi, 1, 0, 0)
        #glTranslatef(-self.com[0], -self.com[1], -self.com[2])

        #glPopMatrix()
        for i in self.sorted_indices:
            atom = self.atoms[i]
            number = atom.number
            if number != last_number:
                set_material(number)
            render_atom(atom)

        self.canvas.SwapBuffers()

    def _gl_scale(self):
        scale = self.scale
        glScale(scale, scale, scale)
        self._gl_draw()


