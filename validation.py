def valid_graph(vertices, edges, vertices_matrix):
    for v in vertices:
        if v not in vertices_matrix:
            raise Exception('vertices list and vertices_matrix is not equal')

    for u, v in edges:
        if u not in vertices_matrix or v not in vertices_matrix[u]:
            raise Exception('vertices list and vertices_matrix is not equal')

    for u in vertices_matrix:
        for v in vertices_matrix[u]:
            if (u, v) not in edges and (v, u) not in edges:
                raise Exception('edge list and vertices_matrix is not equal')
