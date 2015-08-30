#!env python

import graph_tool.all as gt

G = gt.Graph()
label_prop = G.new_vertex_property("string")
color_prop = G.new_vertex_property("vector<double>")

v1 = G.add_vertex()
label_prop[v1] = "orange"
color_prop[v1] = (1,0.5,0,1)

v2 = G.add_vertex()
label_prop[v2] = "red"
color_prop[v2] = (1,0,0,1)

v3 = G.add_vertex()
label_prop[v3] = "green"
color_prop[v3] = (0,1,0,1)

v4 = G.add_vertex()
label_prop[v4] = "blue"
color_prop[v4] = (0,0,1,1)

v5 = G.add_vertex()
label_prop[v5] = "light red"
color_prop[v5] = (1,0,0,.4)

v6 = G.add_vertex()
label_prop[v6] = "light green"
color_prop[v6] = (0,1,0,.3)

v7 = G.add_vertex()
label_prop[v7] = "light blue"
color_prop[v7] = (0,0,1,.3)

G.vertex_properties["vertex_fill_color"] = color_prop
G.vertex_properties["label"] = label_prop

gt.graph_draw(G,vertex_text=G.vertex_properties['label'],
              vertex_fill_color=G.vertex_properties['vertex_fill_color'],
              bg_color=[1,1,1,1],
              output="/tmp/try.png")

