from CRISPResso2.utils.graph import Graph


def test_nodes():
    graph = Graph()
    assert graph.add_node(name='A') == 0
    assert graph.add_node(name='B') == 1
    assert graph.add_node(name='C') == 2
    assert graph.find_nodes(
        ('name', 'A'), ('id', 0),
    ) == [{'name': 'A', 'id': 0}]
    graph.update_node(0, 'name', lambda _: 'D')
    assert graph.find_nodes(
        ('name', 'D'), ('id', 0),
    ) == [{'name': 'D', 'id': 0}]
    graph.update_node(1, 'newField', lambda _: 'new value')
    assert graph.get_node(1) == {
        'name': 'B', 'id': 1, 'newField': 'new value',
    }


def test_edges():
    graph = Graph()
    node_a_id = graph.add_node()
    node_b_id = graph.add_node()
    graph.add_edge(node_a_id, node_b_id, testNum=0)
    assert graph.is_edge(node_a_id, node_b_id)
    assert graph.get_edge(node_a_id, node_b_id) == {
        'source': node_a_id, 'target': node_b_id, 'testNum': 0,
    }
    assert graph.get_outgoing_neighbors(node_a_id) == [{'id': node_b_id}]
    assert graph.get_incoming_neighbors(node_b_id) == [{'id': node_a_id}]
    graph.update_edge(node_a_id, node_b_id, 'testNum', lambda x: x + 2)
    assert graph.get_edge(node_a_id, node_b_id) == {
        'source': node_a_id, 'target': node_b_id, 'testNum': 2,
    }
    graph.update_edge(node_a_id, node_b_id, 'newField', lambda _: 'new value')
    assert graph.get_edge(node_a_id, node_b_id) == {
        'source': node_a_id,
        'target': node_b_id,
        'testNum': 2,
        'newField': 'new value',
    }


def test_paths():
    graph = Graph()
    node_0 = graph.add_node()
    node_1 = graph.add_node()
    node_2 = graph.add_node()
    node_3 = graph.add_node()
    node_4 = graph.add_node()
    graph.add_edge(node_0, node_1)
    graph.add_edge(node_1, node_2)
    graph.add_edge(node_2, node_3)
    graph.add_edge(node_3, node_4)
    assert graph.get_path(node_0, node_1) == [node_0, node_1]
    assert graph.is_path(node_0, node_2)
    assert graph.is_path(node_0, node_3)
    assert graph.is_path(node_0, node_4)
    assert not graph.is_path(node_4, node_0)

    node_5 = graph.add_node()
    assert not graph.is_path(node_4, node_5)
    graph.add_edge(node_2, node_5)
    assert graph.get_path(node_0, node_5) == [node_0, node_1, node_2, node_5]

    graph.add_edge(node_0, node_5)
    assert graph.get_all_paths(node_0, node_5) == [
        [node_0, node_5], [node_0, node_1, node_2, node_5],
    ]
