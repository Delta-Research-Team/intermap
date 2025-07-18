# Created by rglez at 7/18/25
from collections import defaultdict

import networkx as nx
from pyvis.network import Network

from intermap.shiny.app.css import all_interactions_colors as inter_colors


def normalize_to_range(values, target_min=0, target_max=1):
    """
    Normalize a list of values to a predefined range.

    Args:
        values: List of numeric values to normalize
        target_min: Minimum value of target range (default: 0)
        target_max: Maximum value of target range (default: 1)

    Returns:
        List of normalized values in the target range
    """
    if not values:
        return []

    # Find min and max in the original data
    data_min = min(values)
    data_max = max(values)

    # Normalize
    normalized = []
    for value in values:
        normalized_value = (value - data_min) / (data_max - data_min) * (
                target_max - target_min) + target_min
        normalized.append(normalized_value)

    return normalized


class InterNetwork:
    """
    A class to represent the interaction data as a network.
    """

    def __init__(self, master_df, plot_size=(800, 600), node_sizes=(20, 50),
                 edge_widths=(50, 200)):
        """
        Initializes the Network with the provided DataFrame.

        Args:
            master_df (pd.DataFrame): DataFrame containing interaction data.
            plot_size (tuple): Size of the plot (width, height).
            node_sizes (tuple): Minimum and maximum sizes for nodes.
            edge_widths (tuple): Minimum and maximum widths for edges.
        """
        self.master = master_df
        self.min_node_size = node_sizes[0]
        self.max_node_size = node_sizes[1]
        self.min_edge_width = edge_widths[0]
        self.max_edge_width = edge_widths[1]
        self.width = plot_size[0]
        self.height = plot_size[1]

    def get_graph(self):
        """
        Constructs a network graph from the master DataFrame.

        Returns:
            nx.Graph: A NetworkX graph object representing the interactions.
        """
        # Add nodes
        G = nx.Graph()
        sel1_nodes = self.master['sel1'].unique()
        G.add_nodes_from(sel1_nodes, group='sel1')
        receptor_nodes = self.master['sel2'].unique()
        G.add_nodes_from(receptor_nodes, group='sel2')

        # Add edges
        for _, row in self.master.iterrows():
            G.add_edge(row['sel1'], row['sel2'],
                       interaction=row['interaction_name'],
                       color=inter_colors[row['interaction_name']],
                       prevalence=row['prevalence'],
                       length=400 - (row['prevalence'] * 2.5)
                       )

        # Add nodes size based on prevalence
        node_sizes = defaultdict(float)
        for node in G.nodes:
            neighbors = G[node]
            for neighbor in neighbors:
                node_sizes[node] += neighbors[neighbor]['prevalence']
        nodes_scaled = normalize_to_range(list(node_sizes.values()),
                                          self.min_node_size,
                                          self.max_node_size)
        for node, size in zip(G.nodes, nodes_scaled):
            G.nodes[node]['size'] = size

        # Add edges length based on prevalence
        edge_lengths = defaultdict(float)
        for u, v, data in G.edges(data=True):
            edge_lengths[(u, v)] = data['prevalence']
        edges_scaled = normalize_to_range(list(edge_lengths.values()),
                                          self.min_edge_width,
                                          self.max_edge_width)
        for (u, v), width in zip(G.edges, edges_scaled):
            G.edges[u, v]['width'] = width

        return G

    def create_network_plot(self):
        """Create interactive network visualization.

        Returns:
            Network: A pyvis Network object containing the graph.
        """
        G = self.get_graph()
        if not len(G):
            return None

        net = Network(height=f"{self.height}px", font_color='black')
        net.set_options("""
        {
            "nodes": {
                "font": {"size": 12}
            },
            "edges": {
                "font": {"size": 12},
                "smooth": {"type": "continuous"}
            },
            "physics": {
                "enabled": true,
                "forceAtlas2Based": {
                    "gravitationalConstant": -100,
                    "centralGravity": 0.001,
                    "springLength": 200,
                    "springConstant": 0.08,
                    "avoidOverlap": 0.5
                },
                "solver": "forceAtlas2Based",
                "minVelocity": 0.75,
                "stabilization": {
                    "enabled": true,
                    "iterations": 1000
                }
            },
            "interaction": {
                "hover": true,
                "tooltipDelay": 100
            }
        }
        """)

        for node in G.nodes():
            if 'WAT' in node or 'HOH' in node:
                node_color = inter_colors.get('WaterBridge', '#00BFFF')
            else:
                node_color = '#4051b5' if G.nodes[node][
                                              'group'] == 'sel1' else '#ff5555'

            net.add_node(node, label=node, color=node_color,
                         size=G.nodes[node]['size'],
                         )

        for edge in G.edges(data=True):
            net.add_edge(
                edge[0],
                edge[1],
                color=edge[2]['color'],
                label=f"{edge[2]['interaction']} ({edge[2]['prevalence']:.1f}%)",
                length=edge[2]['length'],
                title=(f"Interaction: {edge[2]['interaction']}\n"
                       f"Prevalence: {edge[2]['prevalence']:.1f}%"),
                width=edge[2]['width'],
            )

        return net

# =============================================================================
#
# =============================================================================
# from intermap.shiny.app.icsv import CSVFilter
#
# pickle = '/media/rglez/Roy2TB/Dropbox/RoyData/intermap/tutorial-mayank/outputs/ligs-channel_InterMap.pickle'
# cfg = '/media/rglez/Roy2TB/Dropbox/RoyData/intermap/tutorial-mayank/outputs/ligs-channel_InterMap.cfg'
# csv_obj = CSVFilter(pickle, cfg)
# master_df = csv_obj.master
#
# self = InterNetwork(master_df)
# G = self.get_graph()
