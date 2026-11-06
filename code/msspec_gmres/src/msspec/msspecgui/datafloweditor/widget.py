'''
Created on May 20, 2016

@author: graffy
'''

import wx


class Widget(object):
    '''
    The widget representing an interactive object in the dataflow (an operator, a plug, etc.)
    '''
    
    def __init__(self, data_flow_view):
        """
        :param data_flow_view: the dataflowview to which this operator belongs
        :type data_flow_view: datafloweditor.DataflowView
        """
        self._container_group = None  # the 2d scene graph group node representing this widget
        self._data_flow_view = data_flow_view

    @property
    def scenegraph_group(self):
        """
        :rtype: scenegraph.Group
        """
        return self._container_group
    
    @scenegraph_group.setter
    def scenegraph_group(self, scenegraph_group):
        """
        :type scenegraph_group: scenegraph.Group
        """
        print('Widget.scenegraph_group.setter : scenegraph_group=%s' % str(scenegraph_group))
        self._container_group = scenegraph_group
        self._data_flow_view.scenegraph_group_to_widget[scenegraph_group] = self

    def get_bounding_box(self, border=0):
        """
        :param int border: the width of the border surrounding the box
        :return wx.Rect: the smallest box containing this widget
        """
        (x_min, y_min), (x_max, y_max) = self.scenegraph_group.aabbox()  # (scenegraph_group.setter seen as a method, see https://github.com/PyCQA/pylint/issues/870) pylint: disable=no-member
        return wx.Rect(x_min - border, y_min - border, x_max - x_min + border * 2, y_max - y_min + border * 2)

    def on_hover(self, is_entering):
        '''
        :param is_entering: True is the mouse pointer enters this widget, False if it leaves the widget
        '''
        pass
