import random
import numpy as np
import networkx as nx
import scipy.sparse as ss
import pylab as plt
import seaborn as sns

def random_conf(degree_seq, fun, correct):
    """
    Compare with configuration model (correct self-edges if desired)
    """
    G_r = nx.configuration_model(degree_seq)
    len_full = len(G_r.edges())
    G_r = nx.Graph(G_r)

    if correct:
        # Avoid multiedges (low probability)
        while len_full != len(G_r.edges()):
            G_r = nx.configuration_model(degree_seq)
            len_full = len(G_r.edges())
            G_r = nx.Graph(G_r)
        # Self-edges
        self_edges = list(nx.selfloop_edges(G_r))

        # Start over
        if len(self_edges) == 0:
            pass #great
        elif len(self_edges) % 2 == 1:
            return random_conf(degree_seq, fun, correct)
        else:
            in_n,out_n = zip(*self_edges)
            G_r.remove_edges_from(self_edges)
            G_r.add_edges_from(zip(in_n,out_n[::-1]))
            G_r = nx.Graph(G_r)

            if (len(list(nx.selfloop_edges(G_r))) > 0) | len_full != len(G_r.edges()):
                return random_conf(degree_seq, fun, correct)
        
        
    
    try:
        return fun(G_r)#nx.draw(G_r)
    except:
        return np.nan

def random_g(n, m,fun):
    """
    Compare with ER graph
    """
    G_r = nx.random_graphs.gnm_random_graph(n,m)#nx.random_graphs.erdos_renyi_graph(n, m/((n)*(n-1)/2))
    
    try:
        return fun(G_r)#nx.draw(G_r)
    except:
        return np.nan
    
def random_pl(n, m,fun):
    """
    Compare with BA graph
    """
    G_r = nx.random_graphs.barabasi_albert_graph(n,int(m/n)+1)#nx.random_graphs.erdos_renyi_graph(n, m/((n)*(n-1)/2))
    to_remove = random.sample(list(G_r.edges()),len(G_r.edges())-m)

    G_r.remove_edges_from(to_remove)

    try:
        return fun(G_r)#nx.draw(G_r)
    except:
        return np.nan


def conf_int(G, fun, iterations = 100, correct=True):
    """
    Run all models
    """
    degree_seq = [v for k,v in G.degree()]
    n = len(G)
    m = len(G.edges())
    values = np.array([random_conf(degree_seq, fun, correct) for i in range(iterations)])
    values_c = np.percentile(values[np.isfinite(values)],  [5,95])

    values = np.array([random_g(n, m, fun) for i in range(iterations)])
    values_g = np.percentile(values[np.isfinite(values)],  [5,95])
    
    values = np.array([random_pl(n, m, fun) for i in range(iterations)])
    values_pl = np.percentile(values[np.isfinite(values)],  [5,95])
    
    print(f"Conf. model {values_c[0]:1.3f} - {values_c[1]:1.3f}")
    print(f"ER graph {values_g[0]:1.3f} - {values_g[1]:1.3f}")
    print(f"BA graph {values_pl[0]:1.3f} - {values_pl[1]:1.3f}")
    
# Use the following function to plot the CDF of the degree distributions
def plot_cdf(values, scale = "log", ax = None, cum = True, compl = False, marker = 'o-', xlabel = "Degree (d)", ylabel = "p(Degree < d)"):
    """
    Calculates and plot CDF
    """
    
    from collections import Counter

    # count the number of instance per each degree, sort it
    C = Counter(values)
    deg, cnt = zip(*sorted(C.items()))
    
    # calcualte the cumulative distribution, normalize to be a probability instead of a count
    if cum:
        cs = np.cumsum(cnt)/np.sum(cnt)
    else:
        cs = cnt/np.sum(cnt)
    
    if compl:
        cs = 1 - cs
        
    if ax is None:
        ax = plt.subplot()
    # plot
    ax.plot(deg, cs, marker)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    plt.tight_layout()
    sns.despine(left=True, bottom=True)
    plt.xscale(scale)
    plt.yscale(scale)
    
    
def plot_network_distribution(G, values, mult = 1000, dist=True):
    """
    Plots network (color and node size depends on values) and distributions
    """
    import matplotlib as mpl
    
    norm = mpl.colors.Normalize(vmin=min(values), vmax=max(values), clip=True)
    mapper = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.coolwarm)


    f, (a0, a1, a2) = plt.subplots(1, 3, gridspec_kw={'width_ratios': [2, 1, 1]}, figsize=(12,4))
    mapper._A = []
    cb = plt.colorbar(mapper, ax = a0, location="bottom", shrink=0.8, pad = 0.02, label = "Value")
    cb.outline.set_visible(False)
    
    node_size = mult*np.array(list(values))
    if min(node_size) < 0:
        node_size -= min(node_size)
        node_size += 1
        
    nx.draw(G, pos = nx.spring_layout(G, seed = 1), with_labels = True, node_size = node_size, edge_color = "gray", 
           node_color = [mapper.to_rgba(i) for i in values], ax = a0,)

    if dist:

        sns.histplot(values, ax = a1)
    
    
        plot_cdf(values, ax = a2, compl = True, xlabel = "Cent c", ylabel = "p(Cent > c)")
    
def plot_network(G, a0 = None, values = None, cmap = None, pos = None, with_labels = True, scatter_factor = 1000):
    import matplotlib as mpl
    if cmap is None:
        cmap=mpl.cm.coolwarm
    if values is None:
        values = nx.degree_centrality(G).values()
    
    if a0 is None:
        a0 = plt.gca()
        
    norm = mpl.colors.Normalize(vmin=min(values), vmax=max(values), clip=True)
    mapper = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    mapper._A = []
    cb = plt.colorbar(mapper, ax = a0, location="bottom", shrink=0.8, pad = 0.02, label = "Value")
    cb.outline.set_visible(False)

    # NEtwork
    if nx.is_bipartite(G):
        top = [_ for _ in G.nodes() if _[0] != "S"]
        pos = nx.bipartite_layout(G, top)
        node_color = ["#e6af2e" if node in top else "#e0e2db" for node in G]
    else:
        if pos is None:
            pos = nx.spring_layout(G, seed = 1)
        node_color = [mapper.to_rgba(i) for i in values]

    nx.draw(G, pos = pos, with_labels = with_labels, node_size=scatter_factor*np.array(list(values)), edge_color = "darkgray", 
           node_color = node_color, ax = a0)
    
    
def calculateRWRrange(W, i, alphas, n, maxIter=1000):
    """
    Calculate the personalised TotalRank and personalised PageRank vectors. From Peel et al (2018)
    Parameters
    ----------
    W : array_like
        transition matrix (row normalised adjacency matrix)
    i : int
        index of the personalisation node
    alphas : array_like
        array of (1 - restart probabilties)
    n : int
        number of nodes in the network
    maxIter : int, optional
        maximum number of interations (default: 1000)
    Returns
    -------
    pPageRank_all : array_like
        personalised PageRank for all input alpha values (only calculated if
        more than one alpha given as input, i.e., len(alphas) > 1)
    pTotalRank : array_like
        personalised TotalRank (personalised PageRank with alpha integrated
        out)
    it : int
        number of iterations
    References
    ----------
    See [2]_ and [3]_ for further details.
    .. [2] Boldi, P. (2005). "TotalRank: Ranking without damping." In Special
        interest tracks and posters of the 14th international conference on
        World Wide Web (pp. 898-899).
    .. [3] Boldi, P., Santini, M., & Vigna, S. (2007). "A deeper investigation
        of PageRank as a function of the damping factor." In Dagstuhl Seminar
        Proceedings. Schloss Dagstuhl-Leibniz-Zentrum für Informatik.
    """
    alpha0 = alphas.max()
    WT = alpha0*W.T
    diff = 1
    it = 1

    # initialise PageRank vectors
    pPageRank = np.zeros(n)
    pPageRank_all = np.zeros((n, len(alphas)))
    pPageRank[i] = 1
    pPageRank_all[i, :] = 1
    pPageRank_old = pPageRank.copy()
    pTotalRank = pPageRank.copy()

    oneminusalpha0 = 1-alpha0

    while diff > 1e-9:
        # calculate personalised PageRank via power iteration
        pPageRank = WT @ pPageRank
        pPageRank[i] += oneminusalpha0
        # calculate difference in pPageRank from previous iteration
        delta_pPageRank = pPageRank-pPageRank_old
        # Eq. [S23] Ref. [1]
        pTotalRank += (delta_pPageRank)/((it+1)*(alpha0**it))
        # only calculate personalised pageranks if more than one alpha
        if len(alphas) > 1:
            pPageRank_all += np.outer((delta_pPageRank), (alphas/alpha0)**it)

        # calculate convergence criteria
        diff = np.sum((delta_pPageRank)**2)/n
        it += 1
        if it > maxIter:
            print(i, "max iterations exceeded")
            diff = 0
        pPageRank_old = pPageRank.copy()

    return pPageRank_all, pTotalRank, it

def calculate_local_assort(G, attribute):
    """
    Input:
        G: Undirected graph (networkx)
        attribute: array of values 
    Output:
        loc_ass: array of values representing the local assortativity
    """
    # Create adjacency matrix
    A = nx.to_scipy_sparse_array(G)
    degree = A.sum(1)
    n = len(G)

    # normalize attribute
    attribute = (attribute - np.mean(attribute))/np.std(attribute)

    ## Row normalize matrix
    # construct diagonal inverse degree matrix
    D = ss.diags(1./degree, 0, format='csc')
    # construct transition matrix (row normalised adjacency matrix)
    W = D @ A


    ## Calculate personalized pagerank for all nodes
    pr=np.arange(0., 1., 0.1)
    per_pr = []
    for i in range(n):
        pis, ti, it = calculateRWRrange(W, i, pr, n)
        per_pr.append(ti)
    per_pr = np.array(per_pr)

    # calculate local assortativity (G is undirected, A is symmetric)
    loc_ass = (per_pr * ((A.T * attribute).T * attribute )).sum(1) / degree
    #print(list(sorted(zip(np.round(d,2), G.nodes()))))
    return loc_ass
