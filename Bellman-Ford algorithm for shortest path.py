def bellman_ford(graph, source):
    distance = [float("Inf")] * (len(graph) + 1)
    distance[source] = 0
    pred = [None] * (len(graph) + 1)

    for _ in range(len(graph)):
        for node in range(1, len(graph) + 1):
            for neighbour, cost in graph[node]:
                if distance[node] != float("Inf") and distance[node] + cost < distance[neighbour]:
                    distance[neighbour] = distance[node] + cost
                    pred[neighbour] = node

    for node in range(1, len(graph) + 1):
        for neighbour, cost in graph[node]:
            if distance[node] != float("Inf") and distance[node] + cost < distance[neighbour]:
                return "Графът съдържа цикъл с отрицателна тегловна сума"

    return distance, pred

def print_path(pred, source, target):
    path = []
    while target is not None:
        path.append(target)
        target = pred[target]
    path.reverse()
    print(" -> ".join(str(node) for node in path))

graph = {1: [(2, 2), (3, 2), (4, 5)],
         2: [(5, 8), (6, 12)],
         3: [(4, 3), (5, 6)],
         4: [(5, 10), (7, 15)],
         5: [(6, 6), (7, 8), (8, 7)],
         6: [(8, 4)],
         7: [(8, 5)],
         8: []}

source = 1
target = 8

distance, pred = bellman_ford(graph, source)
print("Дължината на най-късият път е: ", distance[len(distance) - 1])
print("Най-късият път е: ")
print_path(pred, source, target)
