# encoding:utf-8
'''
Created on Feb 19, 2016

@author: graffy
'''
import cairo
import math
import re
from ..scenegraph.element._path import _quadric


def scenegraph_matrix_to_cairo_matrix(sg_matrix):
    """
    :type sg_matrix: scenegraph._Transform
    """
    return cairo.Matrix(*sg_matrix.abcdef)

# from logging import log


def render_group(group, cairo_context):
    """ renders the scenegraph group 'group' in the cairo drawing context 'cairo_context'
    
        :param group: the scenegraph group to render
        :type group: scenegraph2d.Group
        :param cairo_context: the cairo drawing context
        :type cairo_context: cairo.Context
    """
    # print('rendering group')
    cairo_context.transform(scenegraph_matrix_to_cairo_matrix(group.local_to_parent_matrix()))
    for child in group.children:
        render_node(child, cairo_context)


def render_rectangle(rect, cairo_context):
    """ renders the scenegraph rectangle 'rect' in the cairo drawing context 'cairo_context'
    
        :param rect: the scenegraph rectangle to render
        :type rect: scenegraph2d.Rectangle
        :param cairo_context: the cairo drawing context
    """
    # cairo_context.set_source_rgb(0.5, 0.5, 0.5)
    cairo_context.rectangle(rect.x, rect.y, rect.width, rect.height)
    cairo_context.fill()

    # print('rendering rectangle')


def render_path(path, cairo_context):
    """ renders the scenegraph path 'path' in the cairo drawing context 'cairo_context'
    
        :param path: the scenegraph path to render
        :type path: scenegraph2d.Path
        :param cairo_context: the cairo drawing context
        :type cairo_context: cairo.Context
    """
    color = path.fill
    cairo_context.set_source_rgb(*color.rgb)
    path_as_svg_code = path.d
    # print('render_path : rendering path %s' % str(path_as_svg_code))
    paths = []
    path = []
    joins = []
    du2 = 1.
    path_data_iter = iter(path_as_svg_code)
    
    def next_d():
        return next(path_data_iter)
    
    def findCircleCenterPassingThroughPoints(p1, p2, r):
        """
        the set of points (x,y) that are at a distance r from p1=(x1,y1) satisfy :
            (x-x1)^2 + (y-y1)^2 = r^2
            
        the set of points (x,y) that are at a distance r from p2=(x2,y2) satisfy :
            (x-x2)^2 + (y-y2)^2 = r^2
        

        So here's another idea for solving the problem. It should lead to the
        same answer and the algebra will be less tedious.  This method,
        however, uses some fancy ideas from "vector theory", which are
        probably strange and unfamiliar to you.
        
        Think of the geometry. You know that for any two points there's a
        "mirror line" that goes halfway between them. Technically, the line
        consists of the locus of all points that are equidistant from your two
        circle points; you can think of the line as a mirror where each of
        your two points appears as a reflection of the other.
        
        Well, this line will help us a lot in constructing our center, because
        we know that the center is on the line AND because we can use
        Pythagoras to tell us where on the line the point is.  Here's how we
        can do all that with algebra.
        
        First, find the distance between points 1 and 2.  We'll call that q,
        and it's given by sqrt((x2-x1)^2 + (y2-y1)^2).
        
        Second, find the point halfway between your two points.  We'll call it
        (x3, y3).  x3 = (x1+x2)/2  and  y3 = (y1+y2)/2.
        
        Now find the direction of the line between point 1 and point 2.  That
        direction is (x1-x2, y1-y2).
        
        What we really want to know is the direction of the mirror line, which
        is perpendicular to the line between point 1 and point 2.
        
        Here's a crucial trick: you can find the direction of a perpendicular
        to a line just by swapping x with y and changing the sign of one.  In
        other words, if the direction of the joining line was (x1-x2, y1-y2)
        then the direction of the mirror line is (y1-y2,x2-x1).
        
        It will be helpful to "normalize" our direction which means to make
        the length of the line equal to 1.  You do this just by dividing both
        the (y1-y2) and the (x2-x1) by q.
        
        The two circle centers will both be on the mirror line, and we can use
        geometry to find how far they are from the point (x3,y3):  First
        notice that the distance from point (x3,y3) to either point 1 or point
        2 is just half of q.  Then the distance to move along the mirror line
        is:
                              sqrt(r^2-(q/2)^2)
          
        Move up the mirror line by adding a multiple of the direction line to
        the point (x3,y3) or move down the mirror line by subtracting the same
        multiple.
        
        One answer will be:
         
        x = x3 + sqrt(r^2-(q/2)^2)*(y1-y2)/q
        y = y3 + sqrt(r^2-(q/2)^2)*(x2-x1)/q
        
        The other will be:
         
        x = x3 - sqrt(r^2-(q/2)^2)*(y1-y2)/q
        y = y3 - sqrt(r^2-(q/2)^2)*(x2-x1)/q
        
        """
        #
        q = math.sqrt((p2[1] - p1[1]) ** 2 + (p2[0] - p1[0]) ** 2)  # @UnusedVariable
        # let p3 be the midpoint between p1 and p2
        p3 = ((p1[0] + p2[0]) * 0.5, (p1[1] + p2[1]) * 0.5)  # @UnusedVariable
        assert(False)  # todo : finish this implementation
        
    pn = p0 = (0., 0.)
    cn = None
    for c in path_data_iter:
        x0, y0 = p0
        xn, yn = pn
        
        if c.islower():  # coordinates are then relative coordinates
            def next_p():
                dx, dy = next_d()
                return (x0 + dx, y0 + dy)
            
            def next_x():
                dx = next_d()
                return x0 + dx
            
            def next_y():
                dy = next_d()
                return y0 + dy
            c = c.upper()
        else:
            next_x = next_y = next_p = next_d
        
        if c == 'M':  # Moveto
            p1 = next_p()
            # if path:
            #    paths.append((path, False, joins))
            # path = [p1]
            # joins = []
            cairo_context.move_to(*p1)
            
            pn, p0 = p0, p1
        
        elif c in "LHV":
            if c == 'L':  # Lineto
                p1 = next_p()
            elif c == 'H':  # Horizontal Lineto
                p1 = (next_x(), y0)
            elif c == 'V':  # Vertical Lineto
                p1 = (x0, next_y())
            cairo_context.line_to(*p1)
            # path.append(p1)
            pn, p0 = p0, p1
        
        elif c in "CS":  # cubic bezier curve
            if c == 'C':
                p1 = next_p()
            else:  # 'S'
                p1 = (2. * x0 - xn, 2 * y0 - yn) if cn in "CS" else p0
            p2, p3 = next_p(), next_p()
            # path += _cubic(p0, p1, p2, p3, du2)
            cairo_context.rel_curve_to(p1, p2, p3)
            
            pn, p0 = p2, p3
        
        elif c in 'QT':  # quadratic bezier curve
            if c == 'Q':
                p1 = next_p()
            else:  # 'T'
                p1 = (2. * x0 - xn, 2 * y0 - yn) if cn in "QT" else p0
            p2 = next_p()
            path += _quadric(p0, p1, p2, du2)
            pn, p0 = p1, p2
        
        elif c == 'A':  # Arcto
            rs, phi, flags = next_d(), next_d(), next_d()  # @UnusedVariable
            # rs = (rx, ry) : radius in each direction
            # phi = rotation of the axis of the ellipse
            # flags = (large-arc-flag, sweep-flag)
            #   large-arc-flag, indique si on doit afficher l’arc dont la mesure fait plus de la moitié du périmètre de l’ellipse (dans ce cas, la valeur est 1), ou l’arc dont la mesure fait moins de la moitié du périmètre (valeur : 0).
            #   sweep-flag, indique quant à lui si l’arc doit être dessiné dans la direction négative des angles (dans lequel cas sa valeur est 0) ou dans la direction positive des angles (valeur : 1)
            p1 = next_p()
            # p1 : end point
            # path += _arc(p0, rs, phi, flags, p1, du2)
            assert(False)  # incomplete implementation
            # cairo_context.rel_curve_to
            pn, p0 = p0, p1
        
        elif c == 'Z':  # Closepath
            x1, y1 = p1 = path[0]
            dx, dy = x1 - x0, y1 - y0
            if (dx * dx + dy * dy) * du2 > 1.:
                path.append(p1)
            paths.append((path, True, joins))
            path = []
            joins = []
            pn, p0 = p0, p1
        
        cn = c
        joins.append(len(path) - 1)
        
    if path:
        paths.append((path, False, joins))
    cairo_context.stroke()
    cairo_context.fill()

    # print('rendering path')


def render_circle(circle_node, cairo_context):
    """ renders the scenegraph circle 'circle_node' in the cairo drawing context 'cairo_context'
    
        :param circle_node: the scenegraph text to render
        :type circle_node: scenegraph2d.Circle
        :param cairo_context: the cairo drawing context
    """
    color = circle_node.fill
    cairo_context.set_source_rgb(*color.rgb)
    cairo_context.arc(circle_node.cx, circle_node.cy, circle_node.r, 0.0, math.pi * 2.0)
    cairo_context.fill()
    
    # print('render_circle : rendering circle (%f, %f, %f)' % (circle_node.cx, circle_node.cy, circle_node.r))


def render_text(text_node, cairo_context):
    """ renders the scenegraph text 'text' in the cairo drawing context 'cairo_context'
    
        :param text_node: the scenegraph text to render
        :type text_node: scenegraph2d.Text
        :param cairo_context: the cairo drawing context
    """
    # face = wx.lib.wxcairo.FontFaceFromFont(wx.FFont(10, wx.SWISS, wx.FONTFLAG_BOLD))
    # ctx.set_font_face(face)
    cairo_context.set_font_size(10)
    (x, y, width, height, dx, dy) = cairo_context.text_extents(text_node.text)  # @UnusedVariable
    # ctx.move_to(*center.Get())
    # cairo_context.move_to(center.x - width/2, center.y + height/2)
    # cairo_context.set_source_rgb(0, 0, 0)

    cairo_context.move_to(text_node._anchor(), 0.0)
    cairo_context.show_text(text_node.text)

    # print('rendering text')


def render_node(node, cairo_context):
    """ renders the scenegraph node 'node' in the cairo drawing context 'cairo_context'
    
        :param node: the scenegraph node to render
        :type node: scenegraph2d.Element
        :param cairo_context: the cairo drawing context
        :type cairo_context: cairo.Context
    """
    # print node.tag
    try:
        handler = {
            # "svg":      render_group,
            "g":        render_group,
            # "symbol":   render_group,
            # "a":        render_group,
            # "defs":     render_group,
            # "clipPath": render_group,
            # "mask":     render_group,
            # "path":     sg.Path,
            "rect":     render_rectangle,
            "circle":   render_circle,
            # "ellipse":  sg.Ellipse,
            "line":     render_path,
            # "polyline": sg.Polyline,
            # "polygon":  sg.Polygon,
            "text":     render_text,
        }[node.tag]
    except KeyError:
        # log.warning("unhandled tag : %s" % node.tag)
        print("unhandled tag : %s" % node.tag)
        return

    # save the parent transform matrix
    # parent_matrix = cairo_context.get_matrix()
    cairo_context.save()

    if(hasattr(node, 'fill')):
        # print("node.fill = %s" % str(node.fill))
        if (isinstance(node.fill, str)):
            assert(False)  # this shouldn't happen now that fill attribute is set to a Color node
            match = re.match('^rgb\((?P<r>[0-9]+),(?P<g>[0-9]+),(?P<b>[0-9]+)\)', node.fill)
            if match:
                r = float(match.group('r')) / 255.0
                g = float(match.group('g')) / 255.0
                b = float(match.group('b')) / 255.0
            else:
                assert(False)
        else:
            r, g, b = node.fill.get_rgb()
        cairo_context.set_source_rgb(r, g, b)

    # render the node
    handler(node, cairo_context)

    cairo_context.restore()

    # restore the parent transform matrix (so that brother nodes car render relative to their parent)
    # cairo_context.set_matrix(parent_matrix)


def render_scene(scene, cairo_context):
    """ renders the scenegraph 'scene' in the cairo drawing context 'cairo_context'
    
        :param scene: the scenegraph to render
        :type scene: scenegraph2d.Group
        :param cairo_context: the cairo drawing context
        :type cairo_context: cairo.Context
    """
    render_node(scene, cairo_context)
