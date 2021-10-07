"""A graph implementation where data is held for each node and edge."""
from collections import defaultdict, deque
from copy import copy
from typing import Any, Callable, Set, Sequence


class Graph:
    """A graph implementation.

    Nodes are stored in a dictionary, `self.nodes` with the key being the node
    id and the value being an arbitrary payload of data. Edges are likewise
    stored in a dictionary, `self.edges` with the key being the source node id
    and the value being another dictionary that has the target node id as the
    key and the value as an arbitrary payload of data.
    """

    def __init__(self):
        """Create a Graph object."""
        self.nodes = {}
        self.edges = defaultdict(dict)
        self.current_node_id = 0

    def add_node(self, **payload) -> int:
        """Add a node with a payload of keyword arguments.

        Returns
        -------
        node_id: int
            The node id that was just added to the graph.
        """
        payload.update(
            {'id': self.current_node_id},
        )
        self.nodes[self.current_node_id] = payload
        node_id = self.current_node_id
        self.current_node_id += 1
        return node_id

    def _is_valid_node_id(self, node_id: int):
        """Raise a `ValueError` if node id is not in the Graph."""
        if node_id not in self.nodes:
            raise ValueError(
                'The node id: {0} is not in the graph.'.format(node_id),
            )

    def add_edge(self, source_id: int, target_id: int, **payload) -> None:
        """Add an edge with a payload of keyword arguments.

        Parameters
        ----------
        source_id: int
            The id of the source node for this edge.
        target_id: int
            The id of the target node for this edge.

        """
        payload.update(
            {'source': source_id, 'target': target_id}
        )
        self.edges[source_id][target_id] = payload

    def is_edge(self, source_id: int, target_id: int) -> bool:
        """Return True if an edge exists between source id and target id."""
        self._is_valid_node_id(source_id)
        self._is_valid_node_id(target_id)
        if source_id in self.edges and target_id in self.edges[source_id]:
            return True
        return False

    def get_edge(self, source_id: int, target_id: int) -> dict:
        """Return the edge payload between the source id and target id."""
        if source_id in self.edges and target_id in self.edges[source_id]:
            return self.edges[source_id][target_id]
        else:
            raise ValueError(
                'The edge from {0} to {1} does not exist'.format(
                    source_id, target_id,
                )
            )

    def update_edge(
        self,
        source_id: int,
        target_id: int,
        payload_key: Any,
        update_function=None,
    ) -> None:
        """Update the value at `payload_key` for the edge from `source_id` to `target_id`.

        If `payload_key` is already in the payload and `update_function` is
        provided, then the value of `payload_key` will be passed to
        `update_function` as an argument and the return value of
        `update_function` will be assigned to the payload. If not
        `update_function` is provided, then the identity function is used (the
        payload value is not changed).

        If `payload_key` is not in the payload, then `update_function` will be
        passed `None` as its argument. If you believe this is the case, try
        setting `update_function` to `lambda _: <intended value>`.
        """
        def identity(x):
            return x
        if update_function is None:
            update_function = identity
        edge = self.get_edge(source_id, target_id)
        if edge and payload_key in edge:
            self.edges[source_id][target_id].update(
                {payload_key: update_function(edge[payload_key])},
            )
        elif payload_key not in edge:
            self.edges[source_id][target_id].update(
                {payload_key: update_function(None)},
            )

    def update_node(
        self, node_id: int, payload_key: Any, update_function: Callable = None,
    ) -> None:
        """Update the value at `payload_key` for `node_id`.

        If `payload_key` is already in the payload and `update_function` is
        provided, then the value of `payload_key` will be passed to
        `update_function` as an argument and the return value of
        `update_function` will be assigned to the payload. If not
        `update_function` is provided, then the identity function is used (the
        payload value is not changed).

        If `payload_key` is not in the payload, then `update_function` will be
        passed `None` as its argument. If you believe this is the case, try
        setting `update_function` to `lambda _: <intended value>`.

        Parameters
        ----------
        node_id: int
        """
        def identity(x):
            return x
        if update_function is None:
            update_function = identity
        node = self.get_node(node_id)
        if node and payload_key in node:
            self.nodes[node_id].update(
                {payload_key: update_function(node[payload_key])},
            )
        elif payload_key not in node:
            self.nodes[node_id].update(
                {payload_key: update_function(None)}
            )

    def get_node(self, node_id: int) -> dict:
        """Return the payload associated with node id."""
        self._is_valid_node_id(node_id)
        return self.nodes[node_id]

    def find_nodes(self, *payload_keys_and_values):
        """Find node with keys and values in its payload.

        The arguments should be tuples with the first element as the payload
        key and the second element as the corresponding payload value for that
        key. For example,
            ```
            graph.find_nodes(('A', 1), ('B', 2))
            ```
        looks for all nodes that have a value of 1 for the key 'A' and a value
        of 2 for the key 'B'.
        """
        node_payloads = []
        for node_id, payload in self.nodes.items():
            payload_match = all(
                payload_key in payload and
                payload[payload_key] == payload_value
                for payload_key, payload_value in payload_keys_and_values
            )
            if payload_match:
                node_payloads += [payload]
        return node_payloads

    def get_outgoing_neighbors(self, node_id: int) -> Sequence[dict]:
        """Return a list of the node payloads for all neighors that have an edge from node id."""
        return [
            self.nodes[e['target']]
            for e in self.edges[node_id].values()
        ]

    def get_incoming_neighbors(self, node_id: int) -> Sequence[dict]:
        """Return a list of the node payloads for all neighbors that have an edge to node id."""
        return [
            self.nodes[e[node_id]['source']]
            for e in self.edges.values()
            if node_id in e
        ]

    def is_path(self, start_id: int, stop_id: int) -> bool:
        """Return True if there is a path between start id and stop id and False otherwise."""
        if self.get_path(start_id, stop_id):
            return True
        return False

    def get_path(
        self, start_id: int, stop_id: int, include_node_func: Callable = None,
    ) -> Sequence[int]:
        """Return the path between start id and stop id.

        Parameters
        ----------
        start_id: int
            The node id to start the path.
        stop_id: int
            The node id to end the path.
        include_node_func: function
            A function that takes one parameter which is the node payload and
            returns a bool that if True, the node is included in the search and
            if False, the node is not included in the search. The default value
            is None, and this indicates that the function will ignore the
            payload and always return True, thus including all nodes.

        Returns
        -------
        path: Sequence[int]
            A list of node ids that represents a path from start id to stop id.
            Note that both the start and stop ids are included as respectively
            the first and last elements of the path.
        """
        def depth_first_search(current_id, nodes_in_path, visited_nodes):
            visited_nodes.add(current_id)
            nodes_in_path += [current_id]
            if current_id == stop_id or nodes_in_path[-1] == stop_id:
                return nodes_in_path
            if current_id in self.edges:
                for neighbor_id in self.edges[current_id].keys():
                    if (
                        neighbor_id not in visited_nodes and
                        include_node_func(self.get_node(neighbor_id))
                    ):
                        path = depth_first_search(
                            neighbor_id, copy(nodes_in_path), visited_nodes,
                        )
                        if path:
                            return path

        self._is_valid_node_id(start_id)
        self._is_valid_node_id(stop_id)
        if include_node_func is None:
            include_node_func = lambda x: True
        path = depth_first_search(start_id, [], set())
        if (
            path is not None and
            len(path) >= 2 and
            path[0] == start_id and
            path[-1] == stop_id
        ):
            return path
        return []

    def get_all_paths(
        self, start_id: int, stop_id: int,
    ) -> Sequence[Sequence[int]]:
        """Return all of the paths between start id and stop id.

        Parameters
        ----------
        start_id: int
            The node id to start the path.
        stop_id: int
            The node id to end the path.

        Returns
        -------
        paths: Sequence[Sequence[int]]
            A list of lists of node ids that each represent a path from start id
            to stop id. Note that both the start and stop ids are included as
            respectively the first and last elements of the path.
        """
        self._is_valid_node_id(start_id)
        self._is_valid_node_id(stop_id)
        queue = deque()
        queue.append([start_id])
        paths = []
        while queue:
            path = queue.popleft()
            if path[-1] == stop_id:
                paths += [path]
            for neighbor_id in [
                n['id'] for n in self.get_outgoing_neighbors(path[-1])
            ]:
                if neighbor_id not in path:
                    new_path = path.copy()
                    new_path += [neighbor_id]
                    queue.append(new_path)
        return paths

    def filter_paths(
        self, paths: Sequence[Sequence[int]], is_node_valid: Callable,
    ) -> Sequence[Sequence[int]]:
        """Return paths where all nodes of the nodes are valid.

        Parameters
        ----------
        paths: Sequence[Sequence[int]]
            A list of lists of ints that represent the paths represented as
            node ids.
        is_node_valid: Callable
            A function that takes the payload of a node as a parameter and
            returns a bool representing whether the node should be considered
            valid or not.

        Returns
        -------
        All of the paths that are valid.
        """
        return [
            path
            for path in paths
            if all(is_node_valid(self.get_node(node_id)) for node_id in path)
        ]

    def collapse_components(self, split_func: Callable = None) -> None:
        """Calculate which components of the graph can be collapsed.

        A set of nodes can be collapsed, if there is exactly one edge incoming
        and outgoing from each node. The `split_func` parameter can be used to
        add conditions as to when the group should be split.

        Parameters
        ----------
        split_func: function
            A function that has one parameter, which is the node payload and
            returns a bool whether or not the group needs to be split; True if
            the group should be split at this node and False otherwise.
        """
        unvisited_nodes = set(self.nodes.keys())
        if split_func is None:
            split_func = lambda _: True
        self._collapse_components(
            min(unvisited_nodes), set(), unvisited_nodes, split_func,
        )

    def _collapse_components(
        self,
        current_node_id: int,
        current_node_group: Set,
        unvisited_nodes: Set,
        split_func: Callable,
    ) -> None:
        unvisited_nodes.remove(current_node_id)
        out_neighbors = self.get_outgoing_neighbors(current_node_id)
        in_neighbors = self.get_incoming_neighbors(current_node_id)
        current_group_complete = False
        # Expand the current group, or finish the group
        if len(out_neighbors) <= 1 and len(in_neighbors) <= 1 and split_func(
            self.get_node(current_node_id),
        ):
            current_node_group.add(current_node_id)
        elif len(current_node_group) > 1:
            current_group_complete = True
        elif len(current_node_group) == 1:
            current_node_group = set()

        # Find the next node to search
        if out_neighbors and out_neighbors[0]['id'] in unvisited_nodes:
            next_node_id = out_neighbors[0]['id']
        elif in_neighbors and in_neighbors[0]['id'] in unvisited_nodes:
            next_node_id = out_neighbors[0]['id']
        elif unvisited_nodes:
            next_node_id = min(unvisited_nodes)
            if len(current_node_group) > 1:
                current_group_complete = True
            else:
                current_node_group = set()
        else:
            next_node_id = None

        # If the current group is complete, write the ids to all of the nodes
        if current_group_complete:
            current_node_group = list(current_node_group)
            for collapsed_node_id in current_node_group:
                self.update_node(
                    collapsed_node_id,
                    'collapsedNodes',
                    lambda _: current_node_group,
                )
                self.update_node(
                    collapsed_node_id, 'collapsed', lambda _: False,
                )
            current_node_group = set()

        # If there is another node, continue the collapsing
        if next_node_id:
            self._collapse_components(
                next_node_id, current_node_group, unvisited_nodes, split_func,
            )
