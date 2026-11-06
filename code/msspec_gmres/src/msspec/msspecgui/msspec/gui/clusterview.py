
import wx
from msspecgui.msspec.gui.clusterviewer import ClusterViewer
from msspecgui.dataflow import DataFlow
from msspecgui.msspec.gui.viewmanager import IViewCreator


class ClusterView(wx.Window):
    """
    a view that allows the user to view a cluster and also to select which cluster is viewed
    """

    class Creator(IViewCreator):
        
        VIEW_TYPE_NAME = '3d cluster viewer'
        
        def __init__(self, cluster_flow):
            """
            :param msspecgui.msspec.cluster.clusterflow.ClusterFlow cluster_flow: the cluster flow that is associated with the cluster viewers created with this creator
            """
            self._cluster_flow = cluster_flow
        
        @property
        def view_type_name(self):
            """
            :return str:
            """
            return self.VIEW_TYPE_NAME
        
        def create_view(self, parent):
            """
            :param wx.Window parent: the wx.Window that owns the view
            :return wx.Panel:
            """
            return ClusterView(self._cluster_flow, parent)
    
    class ClusterFlowEventsHandler(DataFlow.IDataFlowEventsHandler):
        
        def __init__(self, cluster_view):
            """
            :type cluster_view: ClusterView
            """
            super(ClusterView.ClusterFlowEventsHandler, self).__init__()
            self._cluster_view = cluster_view
            # self._selected_cluster_chc = selected_cluster_chc
            # self.__dataflow = dataflow
        
        def on_added_operator(self, operator):
            """
            :type operator: Operator
            """
            super(ClusterView.ClusterFlowEventsHandler, self).on_added_operator(operator)
            self._cluster_view.update_cluster_cbx()

        def on_deleted_operator(self, operator):
            """
            :param Operator operator:
            """
            super(ClusterView.ClusterFlowEventsHandler, self).on_deleted_operator(operator)
            # the available clusters might have changed
            self._cluster_view.update_cluster_cbx()
            # force an update of the displayed atoms
            if self._cluster_view.selected_cluster is not None:
                self._cluster_view.selected_cluster = self._cluster_view.selected_cluster

        def on_modified_operator(self, operator):
            """
            :type operator: Operator
            """
            super(ClusterView.ClusterFlowEventsHandler, self).on_modified_operator(operator)
            self._cluster_view.update_cluster_cbx()
            # force an update of the displayed atoms
            if self._cluster_view.selected_cluster is not None:
                self._cluster_view.selected_cluster = self._cluster_view.selected_cluster
        
        def on_added_wire(self, wire):
            super(ClusterView.ClusterFlowEventsHandler, self).on_added_wire(wire)
            self._cluster_view.update_cluster_cbx()

        def on_deleted_wire(self, wire):
            super(ClusterView.ClusterFlowEventsHandler, self).on_deleted_wire(wire)
            self._cluster_view.update_cluster_cbx()
            # force an update of the displayed atoms
            if self._cluster_view.selected_cluster is not None:
                self._cluster_view.selected_cluster = self._cluster_view.selected_cluster

    def __init__(self, cluster_dataflow, *args, **kwargs):
        """
        :param msspecgui.msspec.cluster.clusterflow.ClusterFlow cluster_dataflow: the cluster flow to which this view is attached
        """
        super(ClusterView, self).__init__(*args, **kwargs)
        self._cluster_dataflow = cluster_dataflow
        self._available_clusters = None
        self._selected_cluster = None

        main_box = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(main_box)
        
        # add the choice widget that allows the user to select the cluster he wants to view
        widgets_border = 1  # number of pixels used as a border between widgets
        widgets_spacing = 3  # size of the spacings separating widgets, in pixels
        selected_cluster_box = wx.BoxSizer(wx.HORIZONTAL)
        selected_cluster_lbl = wx.StaticText(self, id=-1, label=u'Viewed cluster')
        selected_cluster_box.Add(selected_cluster_lbl)
        selected_cluster_box.AddSpacer(widgets_spacing)
        self._selected_cluster_chc = wx.Choice(self, size=wx.Size(400, -1))
        selected_cluster_box.Add(self._selected_cluster_chc)
        main_box.Add(selected_cluster_box, proportion=0, flag=wx.TOP | wx.BOTTOM, border=widgets_border)
        self._selected_cluster_chc.Bind(wx.EVT_CHOICE, self.on_cluster_selection_changed)

        main_box.AddSpacer(widgets_spacing)
        
        self._cluster_viewer = ClusterViewer(self)
        #self._cluster_viewer.light_mode_threshold = 2
        main_box.Add(self._cluster_viewer, proportion=1, flag=wx.EXPAND | wx.TOP | wx.BOTTOM, border=widgets_border)
        
        self._clusterflow_events_handler = ClusterView.ClusterFlowEventsHandler(self)
        self._cluster_dataflow.add_dataflow_events_handler(self._clusterflow_events_handler)
        self._clusterflow_events_handler.on_added_operator(operator=None)
        self.Bind(wx.EVT_CLOSE, self.on_close)

    def on_close(self, event):
        print("ClusterView.on_close")
        self._cluster_dataflow.remove_dataflow_events_handler(self._clusterflow_events_handler)
        self.Close(True)
        
    @property
    def cluster_viewer(self):
        return self._cluster_viewer
    
    @classmethod
    def _get_cluster_id_string(cls, plug):
        """
        :type plug: Plug
        :rtype: str
        """
        op = plug.operator
        cluster_id = '%s(%d).%s' % (op.creator.get_operator_type_id(), op.id, plug.name)
        return cluster_id
    
    def update_cluster_cbx(self):
        """update the widget that allows the user to choose the viewed cluster
        """
        print('update_cluster_cbx')
        available_clusters = []
        
        for op in self._cluster_dataflow.operators:
            # cluster_id = '%s (%d)' % (op.creator.get_operator_type_id(), op.id)
            for plug in op.get_output_plugs():
                # cluster_id += '.' + plug.name
                # available_clusters.append(self._get_cluster_id_string(plug))
                if plug.value_is_available():  # only propose the clusters that can be computed
                    available_clusters.append(plug)
        # if len(available_clusters) == 0:
            # available_clusters.append('no cluster available')
        self._selected_cluster_chc.SetItems([self._get_cluster_id_string(plug) for plug in available_clusters])
        self._available_clusters = dict(enumerate(available_clusters))
        
        # keep the selected cluster in the choices, if possible
        if self._selected_cluster is not None:
            selected_cluster_indices = [k for k, plug in self._available_clusters.items() if plug == self._selected_cluster]
            if len(selected_cluster_indices) != 0:
                selected_cluster_index = selected_cluster_indices[0]
                self._selected_cluster_chc.SetSelection(selected_cluster_index)
            else:
                # the selected cluster is no longer available (for example because its operator has been deleted)
                self._selected_cluster_chc.SetSelection(wx.NOT_FOUND)
                self.selected_cluster = None

    @property
    def selected_cluster(self):
        """
        :rtype: dataflow.Plug
        """
        return self._selected_cluster

    @selected_cluster.setter
    def selected_cluster(self, selected_cluster):
        """
        :param dataflow.Plug selected_cluster:
        """
        self._selected_cluster = selected_cluster
        if selected_cluster is not None:
            # print('ClusterView.selected_cluster : setting self._selected_cluster to %s' % selected_cluster.name)
        
            # print("ClusterView.selected_cluster : updating viewed atoms")
            self._cluster_viewer.set_atoms(selected_cluster.get_value(), rescale=True)
        else:
            self._cluster_viewer.set_atoms(None, rescale=True)

    def on_cluster_selection_changed(self, event):
        """
        callback when the user changes the selected cluster in the choice widget
        
        :type event: wx.CommandEvent
        """
        if len(self._available_clusters) != 0:
            new_selected_cluster = self._available_clusters[event.GetInt()]
        else:
            new_selected_cluster = None
        self.selected_cluster = new_selected_cluster
