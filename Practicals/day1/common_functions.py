import random
import numpy as np
import networkx as nx
import scipy.sparse as ss
import pylab as plt
import seaborn as sns

# ---
# RANDOM GRAPH GENERATION AND COMPARISON FUNCTIONS
# ---
def random_conf(degree_seq, fun, correct):
    """
    Generate a random graph using the configuration model (preserving the degree sequence),
    and apply a function to it. Optionally correct for self-loops and multiedges.

    Parameters:
        degree_seq (list): Degree sequence to preserve in the random graph.
        fun (function): Function to apply to the generated graph (e.g., a measurement or plot).
        correct (bool): If True, attempts to remove self-loops and multiedges.

    Returns:
        Result of fun(G_r), or np.nan if the function fails.
    """
    G_r = nx.configuration_model(degree_seq)
    len_full = len(G_r.edges())
    G_r = nx.Graph(G_r)  # Convert to simple graph (no parallel edges)

    if correct:
        # Remove multiedges by regenerating until the number of edges matches
        while len_full != len(G_r.edges()):
            # The configuration model can create multiedges; we keep trying until the number of edges matches
            G_r = nx.configuration_model(degree_seq)
            len_full = len(G_r.edges())
            G_r = nx.Graph(G_r)
        # Remove self-loops
        self_edges = list(nx.selfloop_edges(G_r))
        if len(self_edges) == 0:
            pass  # No self-loops, done
        elif len(self_edges) % 2 == 1:
            # If there's an odd number of self-loops, it's not possible to pair them all for rewiring, so try again
            return random_conf(degree_seq, fun, correct)
        else:
            # Remove and rewire self-loops by pairing their endpoints
            in_n, out_n = zip(*self_edges)
            G_r.remove_edges_from(self_edges)
            # Reverse the order of one set to avoid creating new self-loops
            G_r.add_edges_from(zip(in_n, out_n[::-1]))
            G_r = nx.Graph(G_r)
            # If still not correct (self-loops remain or edge count changed), try again
            if (len(list(nx.selfloop_edges(G_r))) > 0) or len_full != len(G_r.edges()):
                return random_conf(degree_seq, fun, correct)
    try:
        return fun(G_r)
    except:
        return np.nan

def random_g(n, m, fun):
    """
    Generate an Erdős–Rényi (ER) random graph with n nodes and m edges, and apply a function to it.

    Parameters:
        n (int): Number of nodes.
        m (int): Number of edges.
        fun (function): Function to apply to the generated graph.

    Returns:
        Result of fun(G_r), or np.nan if the function fails.
    """
    G_r = nx.random_graphs.gnm_random_graph(n, m)
    try:
        return fun(G_r)
    except:
        return np.nan

def random_pl(n, m, fun):
    """
    Generate a Barabási–Albert (BA) scale-free random graph with n nodes and m edges, and apply a function to it.
    Edges are randomly removed to match the exact number of edges m.

    Parameters:
        n (int): Number of nodes.
        m (int): Number of edges.
        fun (function): Function to apply to the generated graph.

    Returns:
        Result of fun(G_r), or np.nan if the function fails.
    """
    G_r = nx.random_graphs.barabasi_albert_graph(n, int(m/n)+1)
    # The BA model may create more edges than needed, so we randomly remove the excess
    to_remove = random.sample(list(G_r.edges()), len(G_r.edges())-m)
    G_r.remove_edges_from(to_remove)
    try:
        return fun(G_r)
    except:
        return np.nan



# ---
# CONFIDENCE INTERVALS FOR NETWORK MEASURES
# ---
def conf_int(G, fun, iterations=100, correct=True):
    """
    Compare a network measure to three random graph models and print 90% confidence intervals for each model.

    Parameters:
        G (networkx.Graph): The original network.
        fun (function): Function to compute a measure on a graph.
        iterations (int): Number of random graphs to generate for each model.
        correct (bool): If True, corrects configuration model for self-loops/multiedges.
    """
    degree_seq = [v for k, v in G.degree()]
    n = len(G)
    m = len(G.edges())
    # Configuration model: preserves degree sequence
    values = np.array([random_conf(degree_seq, fun, correct) for i in range(iterations)])
    values_c = np.percentile(values[np.isfinite(values)], [5, 95])
    # ER random graph: same number of nodes and edges
    values = np.array([random_g(n, m, fun) for i in range(iterations)])
    values_g = np.percentile(values[np.isfinite(values)], [5, 95])
    # BA scale-free graph: same number of nodes and edges
    values = np.array([random_pl(n, m, fun) for i in range(iterations)])
    values_pl = np.percentile(values[np.isfinite(values)], [5, 95])
    print(f"Conf. model {values_c[0]:1.3f} - {values_c[1]:1.3f}")
    print(f"ER graph {values_g[0]:1.3f} - {values_g[1]:1.3f}")
    print(f"BA graph {values_pl[0]:1.3f} - {values_pl[1]:1.3f}")

# ---
# PLOTTING FUNCTIONS
# ---
def plot_cdf(values, scale="log", ax=None, cum=True, compl=False, marker='o-', xlabel="Degree (d)", ylabel="p(Degree < d)"):
    """
    Plot the cumulative distribution function (CDF) or probability distribution of a set of values (e.g., node degrees).

    Parameters:
        values (array-like): Data to plot (e.g., degrees).
        scale (str): Axis scale ('log' or 'linear').
        ax (matplotlib axis): Axis to plot on (optional).
        cum (bool): If True, plot cumulative distribution; else, plot probability distribution.
        compl (bool): If True, plot complementary CDF (1-CDF).
        marker (str): Marker style for the plot.
        xlabel (str): Label for the x-axis.
        ylabel (str): Label for the y-axis.
    """
    from collections import Counter
    # Count occurrences of each value
    C = Counter(values)
    deg, cnt = zip(*sorted(C.items()))
    # Compute cumulative or probability distribution
    if cum:
        cs = np.cumsum(cnt)/np.sum(cnt)
    else:
        cs = cnt/np.sum(cnt)
    if compl:
        cs = 1 - cs
    if ax is None:
        ax = plt.subplot()
    ax.plot(deg, cs, marker)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    plt.tight_layout()
    sns.despine(left=True, bottom=True)
    plt.xscale(scale)
    plt.yscale(scale)

def plot_network_distribution(G, values, mult=1000, dist=True):
    """
    Visualize a network with node color/size based on a value (e.g., centrality),
    and show the value's histogram and CDF.

    Parameters:
        G (networkx.Graph): The network to plot.
        values (array-like): Node values for color/size.
        mult (float): Scaling factor for node size.
        dist (bool): If True, plot histogram and CDF of values.
    """
    import matplotlib as mpl
    norm = mpl.colors.Normalize(vmin=min(values), vmax=max(values), clip=True)
    mapper = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.coolwarm)
    f, (a0, a1, a2) = plt.subplots(1, 3, gridspec_kw={'width_ratios': [2, 1, 1]}, figsize=(12,4))
    mapper._A = []
    cb = plt.colorbar(mapper, ax=a0, location="bottom", shrink=0.8, pad=0.02, label="Value")
    cb.outline.set_visible(False)
    node_size = mult*np.array(list(values))
    if min(node_size) < 0:
        node_size -= min(node_size)
        node_size += 1
    nx.draw(G, pos=nx.spring_layout(G, seed=1), with_labels=True, node_size=node_size, edge_color="gray", 
           node_color=[mapper.to_rgba(i) for i in values], ax=a0)
    if dist:
        sns.histplot(values, ax=a1)
        plot_cdf(values, ax=a2, compl=False, xlabel="Cent c", ylabel="p(Cent > c)")

def plot_network(G, a0=None, values=None, cmap=None, pos=None, with_labels=True, scatter_factor=1000):
    """
    Plot a network with node color and size based on a value (e.g., centrality).

    Parameters:
        G (networkx.Graph): The network to plot.
        a0 (matplotlib axis): Axis to plot on (optional).
        values (array-like): Node values for color/size (default: degree centrality).
        cmap (colormap): Colormap to use.
        pos (dict): Node positions (optional).
        with_labels (bool): Whether to show node labels.
        scatter_factor (float): Scaling factor for node size.
    """
    import matplotlib as mpl
    if cmap is None:
        cmap = mpl.cm.coolwarm
    if values is None:
        values = nx.degree_centrality(G).values()
    if a0 is None:
        a0 = plt.gca()
    norm = mpl.colors.Normalize(vmin=min(values), vmax=max(values), clip=True)
    mapper = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    mapper._A = []
    cb = plt.colorbar(mapper, ax=a0, location="bottom", shrink=0.8, pad=0.02, label="Value")
    cb.outline.set_visible(False)
    # Choose layout and node color
    if nx.is_bipartite(G):
        top = [_ for _ in G.nodes() if _[0] != "S"]
        pos = nx.bipartite_layout(G, top)
        node_color = ["#e6af2e" if node in top else "#e0e2db" for node in G]
    else:
        if pos is None:
            pos = nx.spring_layout(G, seed=1)
        node_color = [mapper.to_rgba(i) for i in values]
    nx.draw(G, pos=pos, with_labels=with_labels, node_size=scatter_factor*np.array(list(values)), edge_color="darkgray", 
           node_color=node_color, ax=a0)

# ---
# RANDOM WALKER AND ASSORTATIVITY FUNCTIONS
# ---
def calculateRWRrange(W, i, alphas, n, maxIter=1000):
    """
    Calculate the personalized TotalRank and personalized PageRank vectors for a node,
    using the power iteration method. Based on Peel et al (2018).

    Parameters:
        W (array): Row-normalized adjacency (transition) matrix.
        i (int): Index of the personalization node.
        alphas (array): Array of (1 - restart probabilities).
        n (int): Number of nodes.
        maxIter (int): Maximum number of iterations.

    Returns:
        pPageRank_all (array): Personalized PageRank for all alpha values.
        pTotalRank (array): Personalized TotalRank (PageRank with alpha integrated out).
        it (int): Number of iterations performed.
    """
    alpha0 = alphas.max()
    WT = alpha0 * W.T  # Transpose for right multiplication
    diff = 1
    it = 1
    # Initialize PageRank vectors
    pPageRank = np.zeros(n)
    pPageRank_all = np.zeros((n, len(alphas)))
    pPageRank[i] = 1  # Start with all probability at node i
    pPageRank_all[i, :] = 1
    pPageRank_old = pPageRank.copy()
    pTotalRank = pPageRank.copy()
    oneminusalpha0 = 1 - alpha0
    while diff > 1e-9:
        # Power iteration: spread probability to neighbors
        pPageRank = WT @ pPageRank
        # Add restart probability (random jump back to node i)
        pPageRank[i] += oneminusalpha0
        delta_pPageRank = pPageRank - pPageRank_old
        # Update TotalRank (integrates over all alpha)
        pTotalRank += (delta_pPageRank) / ((it+1)*(alpha0**it))
        # Update all personalized PageRanks if needed
        if len(alphas) > 1:
            # Outer product distributes delta_pPageRank for each alpha
            pPageRank_all += np.outer((delta_pPageRank), (alphas/alpha0)**it)
        # Check for convergence (L2 norm of change)
        diff = np.sum((delta_pPageRank)**2)/n
        it += 1
        if it > maxIter:
            print(i, "max iterations exceeded")
            diff = 0
        pPageRank_old = pPageRank.copy()
    return pPageRank_all, pTotalRank, it

def calculate_local_assort(G, attribute):
    """
    Calculate the local assortativity for each node in an undirected network.
    Local assortativity measures how similar a node's neighbors are in terms of a given attribute.

    Parameters:
        G (networkx.Graph): Undirected network.
        attribute (array-like): Attribute values for each node (must be in node order).

    Returns:
        loc_ass (np.ndarray): Local assortativity values for each node.
    """
    # Create adjacency matrix
    A = nx.to_scipy_sparse_array(G)
    degree = A.sum(1)
    n = len(G)
    # Normalize attribute to mean 0, std 1
    attribute = (attribute - np.mean(attribute))/np.std(attribute)
    # Row-normalize adjacency matrix (transition matrix)
    D = ss.diags(1./degree, 0, format='csc')
    W = D @ A
    # Calculate personalized PageRank for all nodes
    pr = np.arange(0., 1., 0.1)
    per_pr = []
    for i in range(n):
        # For each node, get the personalized PageRank vector (stationary distribution of a random walker starting at i)
        pis, ti, it = calculateRWRrange(W, i, pr, n)
        per_pr.append(ti)
    per_pr = np.array(per_pr)
    # Compute local assortativity:
    # (per_pr * ((A.T * attribute).T * attribute )).sum(1) / degree
    # This computes, for each node, the weighted sum of attribute similarity with neighbors, normalized by degree
    loc_ass = (per_pr * ((A.T * attribute).T * attribute )).sum(1) / degree
    return loc_ass