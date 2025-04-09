"""Module for graph functions in ligandome calculations."""
from __future__ import annotations

import numpy as np
from tqdm import tqdm
import networkx as nx
from pathlib import Path
from scipy.sparse import lil_matrix
from networkx.classes.graph import Graph
from heapq import heappush, heappop, heapify
from ligandome.utils.sachica import (
    run_sachica,
    peptide_list_to_sachica_input,
    sachica_scores_to_sparse_matrix,
)

def get_isolated_nodes(sparse_edit_distance_matrix: lil_matrix, edit_distance_threshold: int):
    """Get a list of isolated nodes at a given threshold

    Args:
        sparse_edit_distance_matrix (csr_matrix): Sparse edit distance matrix
        peptides (List[str]): List of peptides in correct order for matrix
        
    Returns:
        List[int]: list of peptide indexes which are unconnected at given edit distance
    """
    ones_sparse_matrix = sparse_edit_distance_matrix.copy()
    ones_sparse_matrix[ones_sparse_matrix > edit_distance_threshold] = 0
    return np.argwhere(ones_sparse_matrix.getnnz(0) == 0)

def make_nxgraph_from_sparse_matrix(sparse_edit_distance_matrix: lil_matrix, edit_distance_threshold: int, isolated_nodes: np.array) -> nx.Graph:
    """_summary_

    Args:
        sparse_edit_distance_matrix (lil_matrix): _description_

    Returns:
        Tuple[int, Dict]: _description_
    """
    graph=nx.Graph()
    isolated_nodes_list = isolated_nodes.flatten().tolist()
    isolated_nodes_table = dict([(node, 1) for node in isolated_nodes_list])
    edge_dict = {}

    for i in tqdm(range(sparse_edit_distance_matrix.shape[0])):
        if isolated_nodes_table.get(i) is None:
            graph.add_node(i)
        
    for es, ee in zip(sparse_edit_distance_matrix.nonzero()[0].astype(int), tqdm(sparse_edit_distance_matrix.nonzero()[1].astype(int))):
        if edge_dict.get((es, ee)) is None and edge_dict.get((ee, es)) is None:
            if sparse_edit_distance_matrix[es, ee] <= edit_distance_threshold:
                edge_dict[(es, ee)] = 1
        
    graph.add_edges_from(list(edge_dict.keys()))
    
    return graph

def min_weighted_dominating_set(G: Graph, weight: str | None = None, post_process: bool = True, names_break_ties: bool = False) -> set:
    """Calculate minimum dominating set from a networkX graph.

    This function is an adaptation of the NetworkX implementation with
    a much shorter runtime & produces valid results as verified by
    NetworkX's is_min_dominating_set function. It was sourced from
    https://gist.github.com/bbphd/72543f4d37ac965dd7e4 and minorly 
    adapted. 

    Args:
        G (Graph): NetworkX graph.
        weight (str | None, optional): Name of node weight characteristic. Defaults to None.
        post_process (bool, optional): Whether to further postprocess and minimuse dominating set. Defaults to True.
        names_break_ties (bool, optional): Whether to use names to break tied selections. Defaults to False.

    Returns:
        set: Estimated minimum dominating set of nodes. 
    """    
    pq=[]    # wrapper for heapq from Python 2.7.10 docs section 8.4.2
    entry_finder = {}

    #counter = intertools.count()  don't need the counter really?

    def init_queue(task, priority):  # faster but must use heapify after
        'add tasks to list; must heapify later'
        assert task not in entry_finder
        entry = QueueEntry( task, priority)  
        entry_finder[task] = entry
        pq.append(entry)

    def add_task(task, priority):
        'Add new task or update priority of existing task'
        if task in entry_finder:
            remove_task(task)
        entry = QueueEntry( task, priority)  
        entry_finder[task] = entry
        heappush(pq, entry)
    
    def remove_task(task):
        'Mark existing task as REMOVED.'
        entry = entry_finder.pop(task)
        entry.removed = True

    def pop_task():
        'Remove and return the lowest priority task.'
        while pq:
            entry = heappop(pq)
            task = entry.task
            if not entry.removed :
                del entry_finder[task]
                return task
        raise KeyError('pop from empty priority queue')


    class QueueEntry:
        def __init__(self,task,priority):
            self.task = task
            self.priority = priority
            self.removed = False
        def __lt__(self,other):              #try: not needed if <= only used in heapq
            if not names_break_ties:
                return (self.priority < other.priority)
            else: 
                return ((self.priority, self.task) < (other.priority, other.task)) # this particular choice for historical reasons w/ earlier versions of this program
            


    dom_set = set()      # to be returned

    _coverage_depth = {}   # how many times each node covered so far; 0 if not entered
    
    def _value(u):
        val = len(set(G[u]) | set([u]))    # sometimes self already listed due to self edge
        try: val /= G[u][weight]
        except KeyError: return val
        except ValueError: return float('inf')  # weight must be zero

    def neighbors_incl(u):
        yield u
        for v in set(G[u])-set([u]):
            yield v
    
    def _update_heap_values_now_that(u): # u newly added to dom_set
        delta_val = {}
        for v in neighbors_incl(u):
            if v not in _coverage_depth:  # newly covered
                for w in neighbors_incl(v):
                    delta_val[w] = delta_val.get(w,0) + 1
        for w in set(delta_val.keys()) - dom_set :
            add_task(w, - _find_new_val(w, delta_val[w]))


    def _find_new_val(w, delta_v):
        assert delta_v >= 0
        old_val_tuple = entry_finder[w]
        val = - old_val_tuple.priority
        try: delta_v /= G[w].get(weight)
        except KeyError : pass
        except TypeError : pass
        except ValueError : return float('inf')        # 0/0 = inf seems OK here
        val -= delta_v
        return val
        
    for u in G:           # here think of u as a possible member of dom_set
        init_queue(u, - _value(u))  # value = inverse cost; use - because python heap finds smallest

    heapify(pq)  # faster than doing heappush over and over (says heapq source comments)

    for u in G:           # this time think of u as a node that needs to be covered
        while u not in _coverage_depth:
            #try: 
            v = pop_task()
            #except: 
            #    raise KeyError('no dominating set exists?  directed graphs not supported')

            dom_set.add(v)
            _update_heap_values_now_that(v)
            for w in neighbors_incl(v): 
                _coverage_depth[w] = _coverage_depth.get(w,0) +1

    if post_process:    # Time for Wool and Grossman.
                           # It seems that in the worst case this step can become
                           # the most expensive so we need to use the priority queue all over...
        del pq[:]          # erase old queue
        pq = []

        def _update_heap_redun(u):  # u newly removed
            needs_redun_update = set([])
            for v in neighbors_incl(u): 
                _coverage_depth[v] -= 1
                #for w in (neighbors_incl(v) & dom_set):  oops not legal...do:
                for w in (set(G[v]) | set([v])) & dom_set:
                    needs_redun_update.add(w)
            for v in needs_redun_update:
                old_redun = - (entry_finder[v].priority[0])
                if old_redun != _redun(v):
                    add_task(v,(- _redun(v), - G[u].get(weight,1)))

        def _redun(u):
             retval =  -1 + min((_coverage_depth[v]) for v in neighbors_incl(u))
             if retval < 0: 
                 print("bug?", u, "negative redun")
             return retval
             

        #  tie breaker is weight
        for u in dom_set:  #  set up priority queue
            init_queue(u, (- _redun(u), -G[u].get(weight,1)))

        heapify(pq)

        max_redun = -pq[0].priority[0]
        while max_redun > 0:
            u = pop_task()
            assert _redun(u) == max_redun # just checking
            dom_set.remove(u)
            _update_heap_redun(u)
            max_redun = -pq[0].priority[0]
            
    return dom_set

def calculate_dominating_set_members(peptides: list[str], edit_distance_threshold: int, tmpfile_dir: Path) -> list[bool]:
    """Calculate list of bools for dominating set membership.

    Args:
        peptides (list[str]): Peptides to calculate dominating set for.
        edit_distance_threshold (int): Threshold to draw edges to graph.
        tmpfile_dir (Path): Where to store tempfiles.

    Returns:
        list[bool]: List of dominating set membership booleans of the same size as input peptides list.
    """
    index_table = {peptide: i for i, peptide in enumerate(peptides)}
    assert len(index_table) == len(peptides)

    peptide_tempfile = tmpfile_dir / f'sachica_unique_peptides.fa' 
    
    if not peptide_tempfile.exists():
        peptide_list_to_sachica_input(peptides, peptide_tempfile)

    edges = run_sachica(peptide_tempfile, edit_distance=edit_distance_threshold, temp_dir=tmpfile_dir)
    edit_distance_matrix = sachica_scores_to_sparse_matrix(edges, index_table)

    isolated_nodes = get_isolated_nodes(edit_distance_matrix, 
                                        edit_distance_threshold)
    peptide_graph = make_nxgraph_from_sparse_matrix(edit_distance_matrix, 
                                                    edit_distance_threshold, 
                                                    isolated_nodes)
    dominating_set = min_weighted_dominating_set(peptide_graph, post_process=False)
    return [True if index_table.get(peptide) in dominating_set else False for peptide in peptides]