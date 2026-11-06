# -*- coding: utf-8 -*-

"""
scenegraph.element.group
"""


# imports ####################################################################

from . import Element
from ._path import _bbox


# group ######################################################################

_empty_bbox = _bbox([])


class Group(Element):
    tag = "g"

    _state_attributes = Element._state_attributes + [
        "children",
    ]
        
    def __init__(self, children=None, **attributes):
        super(Group, self).__init__(**attributes)
        self.children = children if children is not None else []
    
    def _aabbox(self, transform, inheriteds):
        bboxes = (child.aabbox(transform, inheriteds) for child in self.children)
        return _bbox(bbox for bbox in bboxes if bbox != _empty_bbox)
    
    def _render(self, transform, inheriteds, context):
        for child in self.children:
            child.render(transform, inheriteds, context)
                    
    def _pick_content(self, x, y, transform):
        hits = []
        for child in self.children:
            # print('Group._pick_content : child = %s' % str(child) )
            hits += child.pick(x, y, transform)
            # print('Group._pick_content : len(hits) = %d' % len(hits) )
            
        return hits
        
    def _xml_content(self, defs):
        return "\n".join(child._xml(defs) for child in self.children)

    def add_child(self, child):
        """
        :param child: the new child to append to the list of children of this svg group
        :type child: a scenegraph.element.Element derived object
        """
        self.children.append(child)
        child.parent = self
        # assert(child.parent == self)
        
    def remove_child(self, child_to_remove):
        self.children = [child for child in self.children if child is not child_to_remove]
