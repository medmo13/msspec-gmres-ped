# -*- coding: utf-8 -*-

"""
scenegraph.element.path
"""


# imports ####################################################################

# from collections import defaultdict
from copy import deepcopy
from math import log, floor, sqrt

# from ...opengl.utils import create_vbo
from . import Element
from ._path import (_cubic, _quadric, _arc, _stroke, _evenodd_hit, _nonzero_hit, _stroke_hit, _bbox)


# flattening #################################################################

def _flatten(path_data, du2=1.):
    """discretize path into straight segments.
    
    :param path_data: path description encoded with svg format, as a list of elements, which are either one-letter keywords or numeric arguments. For example, for a circle, this would be 'M', (cx-r, cy), 'a', (r, r), 0, (0, 0), ( 2*r, 0), 'a', (r, r), 0, (0, 0), (-2*r, 0), 'Z']
    """
    paths = []
    path = []
    joins = []
    
    path_data_iter = iter(path_data)
    
    def next_d():
        return next(path_data_iter)
    
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
            if path:
                paths.append((path, False, joins))
            path = [p1]
            joins = []
            
            pn, p0 = p0, p1
        
        elif c in "LHV":
            if c == 'L':  # Lineto
                p1 = next_p()
            elif c == 'H':  # Horizontal Lineto
                p1 = (next_x(), y0)
            elif c == 'V':  # Vertical Lineto
                p1 = (x0, next_y())
            path.append(p1)
            pn, p0 = p0, p1
        
        elif c in "CS":  # cubic bezier curve
            if c == 'C':
                p1 = next_p()
            else:  # 'S'
                p1 = (2. * x0 - xn, 2 * y0 - yn) if cn in "CS" else p0
            p2, p3 = next_p(), next_p()
            path += _cubic(p0, p1, p2, p3, du2)
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
            rs, phi, flags = next_d(), next_d(), next_d()
            # rs = (rx, ry) : radius in each direction
            # phi = rotation of the axis of the ellipse
            # flags = (large-arc-flag, sweep-flag)
            #   large-arc-flag, indique si on doit afficher l’arc dont la mesure fait plus de la moitié du périmètre de l’ellipse (dans ce cas, la valeur est 1), ou l’arc dont la mesure fait moins de la moitié du périmètre (valeur : 0).
            #   sweep-flag, indique quant à lui si l’arc doit être dessiné dans la direction négative des angles (dans lequel cas sa valeur est 0) ou dans la direction positive des angles (valeur : 1)
            p1 = next_p()
            # p1 : end point
            path += _arc(p0, rs, phi, flags, p1, du2)
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
    
    return paths


# utils ######################################################################

_WIDTH_LIMIT = 1.
_SCALE_STEP = 1.2


def _du2(transform):
    """surface of a pixel in local coordinates."""
    a, b, c, d, _, _ = transform.abcdef
    return abs(a * d - b * c)


def _scale_index(du2, scale_step=_SCALE_STEP):
    """log discretization of the scale suitable as key for hashing cache."""
    try:
        return int(floor(log(du2, scale_step) / 2.))
    except:
        return None


def _strip_range(stop):
    """sort verticies in triangle strip order, i.e. 0 -1 1 -2 2 ..."""
    i = 0
    while i < stop:
        i += 1
        v, s = divmod(i, 2)
        yield v * (s * 2 - 1)


def _join_strips(strips):
    """concatenate strips"""
    strips = iter(strips)
    strip = next(strips, [])
    for s in strips:
        if len(strip) % 2 == 1:
            strip += [strip[-1], s[0], s[0]]
        else:
            strip += [strip[-1], s[0]]
        strip += s
    return strip


# cache ######################################################################

def _fill_state(path):
    return path.d


def _stroke_state(path):
    return (
        path.d,
        path.stroke_width, path.stroke_miterlimit,
        path.stroke_linecap, path.stroke_linejoin,
    )


def _cache(_state):
    """caching decorator
    
    cache is a dict maintained by path element mapping scale index to data
    the cache is cleared if the state characterized by attributes has changed
    """
    def decorator(method):
        def decorated(path, du2=1.):
            state = _state(path)
            if state != path._states.get(method, None):
                path._caches[method] = cache = {}
                path._states[method] = deepcopy(state)
                path._bbox_du2 = 0.
            else:
                cache = path._caches[method]
            scale_index = _scale_index(du2)
            try:
                result = cache[scale_index]
            except KeyError:
                cache[scale_index] = result = method(path, du2)
            return result
        return decorated
    return decorator


# path #######################################################################

class Path(Element):
    tag = "path"
    
    _state_attributes = Element._state_attributes + [
        "d",
    ]
    
    _bbox = (0., 0.), (0., 0.)
    _bbox_du2 = 0.

    def __init__(self, **attributes):
        super(Path, self).__init__(**attributes)
        self._caches = {}
        self._states = {}
    
    def _bbox_is_valid(self):
        return (self._bbox_du2 != 0.)
    
    def _update_bbox(self):
        """ensures that the bounding box of this Path is available
        """
        # print('_update_bbox : self._bbox_du2 = %f' % self._bbox_du2)
        if not self._bbox_is_valid():
            # print("_update_bbox : recomputing self._bbox")
            # self._paths() # this call recomputes the bounding box with the default precision, but it actually doesn't do anything (I suspect the @_cache decorator does something funny)
            du2 = 1.0
            paths = _flatten(self.d, du2)
            self._bbox_du2 = du2
            self._bbox = _bbox(path for (path, _, _) in paths)
    
    @_cache(_fill_state)
    def _paths(self, du2=1.):
        """returns the polygonal approximation of this Path
        
        :param du2: precision of the polygonal approximation
        :returns: a list of polylines
        """
        # print( 'Path._paths : start for %s' % str(self) )
        paths = _flatten(self.d, du2)
        if du2 > self._bbox_du2:
            self._bbox_du2 = du2
            self._bbox = _bbox(path for (path, _, _) in paths)
        return paths
    
    @_cache(_fill_state)
    def _fills(self, du2=1.):
        paths = self._paths(du2)
        return _join_strips([path[i] for i in _strip_range(len(path))]
                            for path, _, _ in paths)
    
    @_cache(_fill_state)
    def _fills_data(self, du2):
        raise NotImplementedError  # only the gl version 
        fills = self._fills(du2)
        return None  # create_vbo(fills)
    
    @_cache(_stroke_state)
    def _strokes(self, du2=1.):
        # print( 'Path._strokes : start for %s' % str(self) )
        paths = self._paths(du2)
        
        # better thin stroke rendering
        du = sqrt(du2)
        adapt_width = self.stroke_width * du
        if adapt_width < _WIDTH_LIMIT:
            width = 1. / du
            opacity_correction = adapt_width
        else:
            width = self.stroke_width
            opacity_correction = 1.
        
        return _join_strips(_stroke(path, closed, joins, width, du,
                                    self.stroke_linecap, self.stroke_linejoin,
                                    self.stroke_miterlimit)
                            for path, closed, joins in paths), opacity_correction
    
    @_cache(_stroke_state)
    def _strokes_data(self, du2):
        strokes, opacity_correction = self._strokes(du2)
        return None, opacity_correction
        # return create_vbo(strokes), opacity_correction
    
    def _aabbox(self, transform, inheriteds):
        du2 = _du2(transform)
        
        points = []
        if self.fill:
            fills = self._fills(du2)
            if fills:
                points.append(transform.project(*p) for p in fills)
        if self.stroke and self.stroke_width > 0.:
            strokes, _ = self._strokes(du2)
            if strokes:
                points.append(transform.project(*p) for p in strokes)
        
        return _bbox(points)
    
    def _render(self, transform, inheriteds, context):
        du2 = _du2(transform)
        origin = self.x, self.y
        
        fill = self._color(self.fill)
        if fill:
            assert False  # doesn't seem to work withoult opengl
            fills = self._fills_data(du2)
            paint = {
                "nonzero": fill.paint_nonzero,
                "evenodd": fill.paint_evenodd,
            }[self.fill_rule]
            paint(self.fill_opacity, fills, transform, context, origin, self._bbox)
        
        stroke = self._color(self.stroke)
        if stroke and self.stroke_width > 0.:
            strokes, correction = self._strokes_data(du2)
            opacity = self.stroke_opacity * correction
            stroke.paint_one(opacity, strokes, transform, context, origin, self._bbox)
    
    def _hit_test(self, x, y, transform):
        """ Tests whether the position (x,y) hits the shape self, when modified with the transformation transform
        """
        # print('Path._hit_test : start : x = %f, y=%f' % (x,y))
        x, y = transform.inverse().project(x, y)
        # print('Path._hit_test : after transform : x = %f, y=%f' % (x,y))
        du2 = _du2(transform)
        hit = False
        
        if not hit and self.fill:
            # print('Path._hit_test : path is filled, self._bbox=%s' % str(self._bbox))
            self._update_bbox()
            (x_min, y_min), (x_max, y_max) = self._bbox
            # print( '_hit_test : (x_min = %f, y_min = %f), (x_max = %f, y_max = %f)' % (x_min, y_min, x_max, y_max) )
            if (x_min <= x <= x_max) and (y_min <= y <= y_max):
                # print('Path._hit_test : in bounding box')
                fills = self._fills(du2)
                if fills:
                    fill_hit = {
                        "nonzero": _nonzero_hit,
                        "evenodd": _evenodd_hit,
                    }[self.fill_rule]
                    hit = fill_hit(x, y, fills)

        if not hit and self.stroke and self.stroke_width > 0.:
            strokes, _ = self._strokes(du2)
            if strokes:
                hit = _stroke_hit(x, y, strokes)
        # print( '_hit_test : hit = %d ' % hit )
        return [([self], (x, y))] if hit else []
