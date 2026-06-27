"""Pure-Python maximum-parsimony tree inference.

Parsimony scoring uses Fitch's algorithm; the tree search is stepwise addition
followed by nearest-neighbour-interchange (NNI) hill-climbing. Missing/ambiguous
bases are treated as the full base set, so they never add cost — the standard
handling for gappy SNP alignments.

The search runs on a compact internal unrooted-tree representation (an adjacency
map) for speed, then exports the result as a DendroPy tree so it shares output
plumbing with neighbor-joining. :func:`score_tree` scores any DendroPy tree and
is the canonical public scorer.
"""

from __future__ import annotations

import dendropy

from .snps import VALID_BASES, SnpMatrix

FULL = frozenset(VALID_BASES)


def _leaf_state(base: str) -> frozenset[str]:
    return frozenset([base]) if base in FULL else FULL


# -- Fitch scoring on a DendroPy tree (public API) ---------------------------


def score_tree(tree: dendropy.Tree, m: SnpMatrix) -> int:
    """Total Fitch parsimony score of ``tree`` over all sites in ``m``."""
    index = {name: i for i, name in enumerate(m.sample_names)}
    seed = tree.seed_node
    total = 0
    for j in range(m.n_loci):
        total += _fitch_dendropy(seed, index, m, j)
    return total


def _fitch_dendropy(node, index, m: SnpMatrix, j: int) -> int:
    """Post-order Fitch over a DendroPy subtree; returns the score for site j."""
    score = 0

    def rec(n) -> frozenset[str]:
        nonlocal score
        children = n.child_nodes()
        if not children:
            return _leaf_state(m.matrix[index[n.taxon.label], j])
        sets = [rec(c) for c in children]
        inter = frozenset.intersection(*sets)
        if inter:
            return inter
        score += 1
        return frozenset.union(*sets)

    rec(node)
    return score


# -- Internal unrooted-tree search -------------------------------------------


class _Tree:
    """Minimal unrooted binary tree as ``node -> set(neighbours)``.

    Leaves are integers ``0..n-1`` (sample indices); internal nodes are integers
    ``>= n``.
    """

    def __init__(self, n: int) -> None:
        self.n = n
        self.adj: dict[int, set[int]] = {}
        self._next_internal = n

    def _new_internal(self) -> int:
        node = self._next_internal
        self._next_internal += 1
        self.adj[node] = set()
        return node

    def _link(self, a: int, b: int) -> None:
        self.adj.setdefault(a, set()).add(b)
        self.adj.setdefault(b, set()).add(a)

    def _unlink(self, a: int, b: int) -> None:
        self.adj[a].discard(b)
        self.adj[b].discard(a)

    def edges(self) -> list[tuple[int, int]]:
        seen = set()
        out = []
        for a, nbrs in self.adj.items():
            for b in nbrs:
                key = (a, b) if a < b else (b, a)
                if key not in seen:
                    seen.add(key)
                    out.append(key)
        return out

    # Fitch scoring on the internal tree -----------------------------------

    def score(self, columns: list[list[str]]) -> int:
        """Score all sites. ``columns[j][leaf]`` is the base char for that leaf."""
        root = next(iter(self.adj))
        total = 0
        for col in columns:
            total += self._fitch(root, col)
        return total

    def _fitch(self, root: int, col: list[str]) -> int:
        score = 0
        # iterative post-order to avoid recursion overhead/limits
        stack = [(root, -1, False)]
        state: dict[int, frozenset[str]] = {}
        while stack:
            node, parent, processed = stack.pop()
            children = [c for c in self.adj[node] if c != parent]
            if not children:
                state[node] = _leaf_state(col[node])
                continue
            if not processed:
                stack.append((node, parent, True))
                for c in children:
                    stack.append((c, node, False))
            else:
                sets = [state[c] for c in children]
                inter = frozenset.intersection(*sets)
                if inter:
                    state[node] = inter
                else:
                    score += 1
                    state[node] = frozenset.union(*sets)
        return score

    def insert_leaf_on_edge(self, leaf: int, edge: tuple[int, int]) -> int:
        """Insert ``leaf`` by splitting ``edge`` with a new internal node."""
        u, v = edge
        w = self._new_internal()
        self._unlink(u, v)
        self._link(u, w)
        self._link(v, w)
        self._link(w, leaf)
        return w

    def to_newick(self, names: list[str]) -> str:
        # root at an internal node for a rooted-looking newick (topology only)
        root = next(node for node in self.adj if node >= self.n)

        def rec(node: int, parent: int) -> str:
            children = [c for c in self.adj[node] if c != parent]
            if not children:
                return names[node]
            inner = ",".join(rec(c, node) for c in children)
            return f"({inner})"

        return rec(root, -1) + ";"


def _build_columns(m: SnpMatrix) -> list[list[str]]:
    return [[m.matrix[i, j] for i in range(m.n_samples)] for j in range(m.n_loci)]


def _stepwise_addition(m: SnpMatrix, columns: list[list[str]]) -> _Tree:
    n = m.n_samples
    tree = _Tree(n)
    # seed star with first three taxa (or two for n==2)
    hub = tree._new_internal()
    for leaf in range(min(3, n)):
        tree._link(hub, leaf)
    for leaf in range(3, n):
        best_edge = None
        best_score = None
        for edge in tree.edges():
            w = tree.insert_leaf_on_edge(leaf, edge)
            sc = tree.score(columns)
            # undo
            u, v = edge
            tree._unlink(w, leaf)
            tree._unlink(u, w)
            tree._unlink(v, w)
            del tree.adj[w]
            tree._next_internal -= 1
            tree._link(u, v)
            if best_score is None or sc < best_score:
                best_score = sc
                best_edge = edge
        tree.insert_leaf_on_edge(leaf, best_edge)
    return tree


def _nni_search(tree: _Tree, columns: list[list[str]], max_rounds: int = 20) -> _Tree:
    """Greedy NNI hill-climbing minimising parsimony score."""
    best = tree.score(columns)
    for _ in range(max_rounds):
        improved = False
        for u, v in tree.edges():
            if u < tree.n or v < tree.n:
                continue  # internal edges only
            u_nbrs = [x for x in tree.adj[u] if x != v]
            v_nbrs = [x for x in tree.adj[v] if x != u]
            if len(u_nbrs) != 2 or len(v_nbrs) != 2:
                continue
            # two NNI moves: swap one u-side subtree with one v-side subtree
            for a in (u_nbrs[0],):
                for b in (v_nbrs[0], v_nbrs[1]):
                    tree._unlink(u, a)
                    tree._unlink(v, b)
                    tree._link(u, b)
                    tree._link(v, a)
                    sc = tree.score(columns)
                    if sc < best:
                        best = sc
                        improved = True
                    else:
                        # revert
                        tree._unlink(u, b)
                        tree._unlink(v, a)
                        tree._link(u, a)
                        tree._link(v, b)
        if not improved:
            break
    return tree


def parsimony_tree(m: SnpMatrix) -> tuple[dendropy.Tree, int]:
    """Infer a maximum-parsimony tree; returns (tree, parsimony_score)."""
    columns = _build_columns(m)
    tree = _stepwise_addition(m, columns)
    if m.n_samples >= 4:
        tree = _nni_search(tree, columns)
    score = tree.score(columns)
    newick = tree.to_newick(list(m.sample_names))
    dp = dendropy.Tree.get(data=newick, schema="newick", taxon_namespace=dendropy.TaxonNamespace())
    return dp, score
