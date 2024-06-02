#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
//#include "MyQueque.h"

using namespace std;

#define UNIMPLMENTED 7
const char UNCOLORED_BIPARTITIE = 10;
const unsigned int MAX_CAPACITY_UINT = 4294967295 - 1;


class NodeQueue {
public:
    int data;
    NodeQueue* next;

    NodeQueue(int val) {
        data = val;
        next = nullptr;
    }
};
class MyQueue {
private:
    NodeQueue* front;
    NodeQueue* rear;
    int size;

public:
    MyQueue() {
        front = nullptr;
        rear = nullptr;
        size = 0;
    }

    // Enqueue function to add an element 
    void enqueue(int val) {
        NodeQueue* newNode = new NodeQueue(val);
        if (isEmpty()) {
            front = newNode;
        }
        else {
            rear->next = newNode;
        }
        rear = newNode;
        size++;
    }

    // Dequeue function to remove and return 
    int dequeue() {
        if (isEmpty()) {
            return -1;
        }
        NodeQueue* tmp = front;
        int dequeuedValue = tmp->data;
        front = front->next;
        delete tmp;
        size--;
        if (front == nullptr) {
            rear = nullptr;
        }
        return dequeuedValue;
    }

    bool isEmpty() const {
        return front == nullptr;
    }
};


void PrintGraph(int** graph, int* nr_of_edges, int graph_size) {
    cout << endl;
    for (int i = 0; i < graph_size; i++) {
        cout << "Vertex " << i << ": ";
        for (int j = 0; j < nr_of_edges[i]; j++) {
            cout << graph[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void swap(int*& a, int*& b) {
    int* temp = a;
    a = b;
    b = temp;
}

void heapify(int** arr, int n, int i) {
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;
    //change it other if deacrasing
    if (left < n && (arr[left][1] < arr[largest][1] ||
        (arr[left][1] == arr[largest][1] && arr[left][0] > arr[largest][0]))) {
        largest = left;
    }

    ///change it other if deacrasing in the paretneces
    if (right < n && (arr[right][1] < arr[largest][1] ||
        (arr[right][1] == arr[largest][1] && arr[right][0] > arr[largest][0]))) {
        largest = right;
    }

    if (largest != i) {
        swap(arr[i], arr[largest]);
        heapify(arr, n, largest);
    }
}

void heapSort(int** arr, int n) {
    for (int i = n / 2 - 1; i >= 0; i--) {
        heapify(arr, n, i);
    }

    for (int i = n - 1; i >= 0; i--) {
        swap(arr[0], arr[i]);
        heapify(arr, i, 0);
    }
}

void PrintDegreeSequence(int** nr_of_edges_sorted, int graph_size) {
    heapSort(nr_of_edges_sorted, graph_size);
    /*
    for (int i = 0; i < graph_size; i++) {
        cout << "Edge: " << nr_of_edges_sorted[i][0] << " degre: " << nr_of_edges_sorted[i][1] << endl;
    }
    cout << endl;*/

    for (int i = 0; i < graph_size; i++) {
        printf("%d ", nr_of_edges_sorted[i][1]);
    }
    printf("\n");
}

bool DFS_Components_Bipartite(int** graph, int* nr_of_edges, char* colored_graph, bool* visited, int start) {
    bool bipartite = true;
    MyQueue que;
    que.enqueue(start);

    while (!que.isEmpty()) {
        int current = que.dequeue();
        visited[current] = true;

        for (int i = 0; i < nr_of_edges[current]; i++) {
            if (colored_graph[graph[current][i]] == UNCOLORED_BIPARTITIE) {
                colored_graph[graph[current][i]] = 1 - colored_graph[current];
                que.enqueue(graph[current][i]);
            }
            else if (colored_graph[graph[current][i]] == colored_graph[current]) {
                bipartite = false;
            }
            else if (graph[current][i] == current) {
                bipartite = false;
            }
        }
    }

    return bipartite;
}

void ComponentsAndBipartite(int** graph, int* nr_of_edges, int graph_size, int* nr_of_components) {
    bool* visited = new bool[graph_size];
    char* colored_graph = new char[graph_size];
    for (int i = 0; i < graph_size; i++) {
        visited[i] = false;
        colored_graph[i] = UNCOLORED_BIPARTITIE;
    }
    //int nr_of_components = 0;
    bool bipartite = true;

    for (int i = 0; i < graph_size; i++) {
        if (!visited[i]) {
            colored_graph[i] = 1;
            if (!DFS_Components_Bipartite(graph, nr_of_edges, colored_graph, visited, i)) {
                bipartite = false;
            }
            (*nr_of_components)++;
        }
    }
    printf("%d\n", *nr_of_components);
    if (bipartite) {
        printf("T\n");
    }
    else {
        printf("F\n");
    }
    delete[] visited;
    delete[] colored_graph;
}

void BFS(int** graph, int* nr_of_edges, bool* visited, int start, int* size_of_the_component, MyQueue& vertecies_in_component) {
    MyQueue que;
    que.enqueue(start);
    visited[start] = true;
    vertecies_in_component.enqueue(start);
    //vertecies_in_component[0] = start;
    (*size_of_the_component)++;

    while (!que.isEmpty()) {
        int current = que.dequeue();
        for (int i = 0; i < nr_of_edges[current]; i++) {
            if (!visited[graph[current][i]]) {
                que.enqueue(graph[current][i]);
                vertecies_in_component.enqueue(graph[current][i]);
                //vertecies_in_component[*size_of_the_component] = graph[current][i];
                (*size_of_the_component)++;
                visited[graph[current][i]] = true;

            }
        }
    }
}

unsigned int BFS_AllDistances(int** graph, int* nr_of_edges, int start, unsigned int* dist, unsigned int* max_distances, int size_of_the_component) {
    int* queue = new int[size_of_the_component];
    int front = 0, rear = 0;
    
    dist[start] = 0;
    queue[rear++] = start;
    unsigned int max_distance = 0;
    int visited_count = 1;

    while (front < rear) {
        int node = queue[front++];
        for (int i = 0; i < nr_of_edges[node]; i++) {
            int neighbor = graph[node][i];
            if (dist[neighbor] == MAX_CAPACITY_UINT) {
                dist[neighbor] = dist[node] + 1;
                queue[rear++] = neighbor;
                visited_count++;
                if (dist[neighbor] > max_distance) {
                    max_distance = dist[neighbor];
                }
                // Update the maximum distance for the ending vertex
                if (dist[neighbor] > max_distances[neighbor]) {
                    max_distances[neighbor] = dist[neighbor];
                }
                // Early termination if we visited all nodes in the component
                if (visited_count == size_of_the_component) {
                    delete[] queue;
                    return max_distance;
                }
            }
        }
    }

    delete[] queue;
    return max_distance;
}

void Eccentricity(int** graph, int* nr_of_edges, int graph_size, int nr_of_components) {
    unsigned int* max_distances = new unsigned int[graph_size];
    for (int i = 0; i < graph_size; i++) {
        max_distances[i] = MAX_CAPACITY_UINT;
    }

    bool* visited = new bool[graph_size]();
    int visited_vertices = 0;
    int last_visited = 0;
    for (int current_component = 0; current_component < nr_of_components; current_component++) {
        int size_of_the_component = 0;
        MyQueue vertices_in_component;// = new int[graph_size];

        for (int i = last_visited; i < graph_size; i++) {
            if (max_distances[i] == MAX_CAPACITY_UINT) {
                BFS(graph, nr_of_edges, visited, i, &size_of_the_component, vertices_in_component);
                visited_vertices += size_of_the_component;
                last_visited = i;
                break;
            }
        }
        bool is_not_alone = true;
        if (size_of_the_component == 1) {
            is_not_alone = false;
            max_distances[vertices_in_component.dequeue()] = 0;
        }
        for (int i = 0; i < size_of_the_component && is_not_alone; i++) {
            int current_vertex = vertices_in_component.dequeue();
            if (max_distances[current_vertex] == MAX_CAPACITY_UINT) { // Skip already calculated vertices
                unsigned int* distances = new unsigned int[graph_size];
                for (int j = 0; j < graph_size; j++) {
                    distances[j] = MAX_CAPACITY_UINT;
                }
                unsigned int max_distance = BFS_AllDistances(graph, nr_of_edges, current_vertex, distances, max_distances, size_of_the_component);
                max_distances[current_vertex] = max_distance;
                delete[] distances;
            }
        }

        //delete[] vertices_in_component;

        // Early termination if all vertices in the graph have been visited
        if (visited_vertices == graph_size) {
            break;
        }
    }

    for (int i = 0; i < graph_size; ++i) {
        printf("%d ", max_distances[i]);
    }
    printf("\n");

    delete[] visited;
    delete[] max_distances;
}

void Greedy(int** graph, const int* nr_of_edges, const int graph_size) {
    int* result = new int[graph_size];
    bool* available = new bool[graph_size];
    result[0] = 1;

    for (int i = 1; i < graph_size; i++) {
        result[i] = -1;
        available[i] = false;
    }

    for (int i = 0; i < graph_size; i++) {

        for (int j = 0; j < nr_of_edges[i]; j++) {
            if (result[graph[i][j]] != -1) {
                available[result[graph[i][j]]] = true;
            }
        }

        int color;
        for (color = 1; color < graph_size + 1; color++)
            if (available[color] == false)
                break;

        result[i] = color;

        for (int j = 0; j < nr_of_edges[i]; j++) {
            if (result[graph[i][j]] != -1) {
                available[result[graph[i][j]]] = false;
            }
        }
    }
    for (int i = 0; i < graph_size; ++i) {
        printf("%d ", result[i]);
    }
    printf("\n");
    delete[] result;
    delete[] available;
}

void LF(int** graph, int** nr_of_edges_sorted, int* nr_of_edges, const int graph_size) {
    int* result = new int[graph_size];
    bool* available = new bool[graph_size];


    for (int i = 0; i < graph_size; i++) {
        result[i] = -1;
        available[i] = false;
    }

    int f = nr_of_edges_sorted[0][0];
    result[f] = 1;

    for (int i = 0; i < graph_size; i++) {
        int vertex = nr_of_edges_sorted[i][0];

        for (int j = 0; j < nr_of_edges[vertex]; j++) {
            int adjacent_vertex = graph[vertex][j];
            if (result[adjacent_vertex] != -1) {
                available[result[adjacent_vertex]] = true;
            }
        }

        int color;
        for (color = 1; color <= graph_size; color++) {
            if (!available[color]) {
                break;
            }
        }

        result[vertex] = color;

        for (int j = 0; j < nr_of_edges[vertex]; j++) {
            int adjacent_vertex = graph[vertex][j];
            if (result[adjacent_vertex] != -1) {
                available[result[adjacent_vertex]] = false;
            }
        }
    }

    for (int i = 0; i < graph_size; ++i) {
        printf("%d ", result[i]);
    }
    printf("\n");

    delete[] result;
    delete[] available;
}

void SLF(int** graph, int** nr_of_edges_sorted, int* nr_of_edges, const int graph_size, const int nr_of_components) {
    int* result = new int[graph_size];
    bool* available = new bool[nr_of_edges_sorted[0][1]];
    int* saturation = new int[graph_size];
    bool* visited = new bool[graph_size];


    // Initialize result, available, saturation, and visited arrays
    for (int i = 0; i < graph_size; i++) {
        result[i] = -1;
        saturation[i] = 0;
        visited[i] = false;
    }

    // Color remaining vertices
    for (int i = 0; i < nr_of_components; i++) {
        //starting points for each component
        int first_vertex = -1;
        //tutaj mozna isc po size of the component
        for (int j = 0; j < graph_size; j++) {
            int p = nr_of_edges_sorted[j][0];
            if (!visited[nr_of_edges_sorted[j][0]]) {
                first_vertex = nr_of_edges_sorted[j][0];
                visited[first_vertex] = true;
                result[first_vertex] = 0;
                break;
            }
        }
        //zmienic przydzielanie saturacji
        for (int j = 0; j < nr_of_edges[first_vertex]; j++) {
            int neighbor = graph[first_vertex][j];
            if (result[neighbor] == -1) {
                saturation[neighbor]++;
            }
        }
        cout << endl;
        cout << "Commponent: " << i << endl;
        cout << "Firest vertex of the component Vertex index: " << first_vertex << " Saturation: " << saturation[first_vertex] << " Degree: " << nr_of_edges[first_vertex] << endl;
        int max_color = 1;
        while (true && nr_of_edges[first_vertex] != 0) {
            int vertex_to_color = -1;
            int max_saturation = -1;
            int max_degree = -1;

            // Find the vertex with the highest degree among those with saturation greater than zero

            for (int j = 0; j < graph_size; j++) {
                // szukaj posrod nie odwiedzonych i ktore otrzymaly juz saturacje tzn byly juz widziane
                if (!visited[j] && saturation[j] > 0 && result[j] == -1) {
                    cout << "Vertex index: " << j << " Saturation: " << saturation[j] << " Degree: " << nr_of_edges[j] << endl;
                    if (saturation[j] > max_saturation) {
                        max_degree = nr_of_edges[j];
                        max_saturation = saturation[j];
                        vertex_to_color = j;
                    }
                    else if (saturation[j] == max_saturation) {
                        if (nr_of_edges[j] > max_degree) {
                            max_degree = nr_of_edges[j];
                            vertex_to_color = j;
                        }
                        else if (nr_of_edges[j] == max_degree) {
                            if (j < vertex_to_color) {
                                vertex_to_color = j;
                            }
                        }
                    }
                }
            }
            if (vertex_to_color == -1) {
                break;
            }

            // Find the smallest available color
            for (int j = 0; j < nr_of_edges_sorted[0][1]; j++) {
                available[j] = true;
            }
            for (int j = 0; j < nr_of_edges[vertex_to_color]; j++) {
                int neighbor = graph[vertex_to_color][j];
                if (result[neighbor] != -1) {
                    available[result[neighbor]] = false;
                }
            }

            int color;
            for (color = 0; color < graph_size; color++) {
                if (available[color]) {
                    break;
                }
            }
            cout << endl;
            cout << "Vertex chosen: " << vertex_to_color << " Color: " << color << " Saturation: " << saturation[vertex_to_color] << " Degree: " << nr_of_edges[vertex_to_color] << endl;
            cout << endl;
            // Assign the color and update the saturation of the neighbors
            //cout << vertex_to_color << " ";
            result[vertex_to_color] = color;
            visited[vertex_to_color] = true;


            /*
            if (color == max_color) {
                for (int k = 0; k < nr_of_edges[vertex_to_color]; k++) {
                    int neighbour = graph[vertex_to_color][k];
                    saturation[neighbour]++;
                }
                max_color++;
            }
            else {
                for (int k = 0; k < nr_of_edges[vertex_to_color]; k++) {
                    int neighbour = graph[vertex_to_color][k];
                    int color_repeat = 0;
                    for (int l = 0; l < nr_of_edges[neighbour]; l++) {
                        int neighbour_of_neighbour = graph[neighbour][l];
                        if (result[neighbour_of_neighbour] == color) {
                            color_repeat++;
                            break;
                        }
                    }
                    if (color_repeat > 1) {
                        saturation[neighbour]++;
                    }
                }
            }*/
            //jesli 

            //ide do sasiadow kolorowanego verexu

        }
        //cout << endl;

    }

    // Output the result
    for (int i = 0; i < graph_size; i++) {
        cout << result[i] + 1 << " ";
    }
    cout << endl;

    // Clean up
    delete[] result;
    delete[] available;
    delete[] saturation;
    delete[] visited;
}

void C4(int** graph, int* nr_of_edges, const int graph_size) {
    long long int C4_sum = 0;
    bool* not_possible = new bool[graph_size]();

    for (int main = 0; main < graph_size; main++) {
        if (nr_of_edges[main] < 2) {
            not_possible[main] = true;
            continue;
        }

        for (int i = 0; i < nr_of_edges[main]; i++) {

            int neighbour_one = graph[main][i];
            if (nr_of_edges[neighbour_one] < 2 || not_possible[neighbour_one]) {
                continue;
            }

            for (int j = i + 1; j < nr_of_edges[main]; j++) {
                int neighbour_two = graph[main][j];
                if (nr_of_edges[neighbour_two] < 2 || not_possible[neighbour_two]) {
                    continue;
                }
                bool* neighbours = new bool[graph_size]();
                for (int k = 0; k < nr_of_edges[neighbour_one]; k++) {
                    neighbours[graph[neighbour_one][k]] = true;
                }
                for (int k = 0; k < nr_of_edges[neighbour_two]; k++) {
                    if (neighbours[graph[neighbour_two][k]] && graph[neighbour_two][k] > main) {
                        C4_sum++;
                    }
                }
                delete[] neighbours;
            }


        }
        not_possible[main] = true;
    }
    delete[] not_possible;
    // Print the number of 4-cycles
    cout << C4_sum << endl;
}

void ComplementEdges(const int* nr_of_edges, const int graph_size) {
    long long int sum = 0;
    for (int i = 0; i < graph_size; i++) {
        sum += nr_of_edges[i];
    }
    long long int total = (long long int)graph_size * ((long long int)graph_size - 1) / 2 - sum / 2;
    printf("%lld\n", total);
}

int main() {
    int nr_of_graphs;
    scanf("%d", &nr_of_graphs);

    for (int current_graph = 0; current_graph < nr_of_graphs; current_graph++) {
        int graph_size;
        scanf("%d", &graph_size);
        int** graph = (int**)malloc(graph_size * sizeof(int*));

        int* nr_of_edges = (int*)malloc(graph_size * sizeof(int));
        for (int current_line = 0; current_line < graph_size; current_line++) {
            scanf("%d", &nr_of_edges[current_line]);
            graph[current_line] = (int*)malloc(nr_of_edges[current_line] * sizeof(int));
            for (int current_edge = 0; current_edge < nr_of_edges[current_line]; current_edge++) {
                int connection;
                scanf("%d", &connection);
                graph[current_line][current_edge] = connection - 1;
            }
        }
        int nr_of_components = 0;
        int** nr_of_edges_sorted = new int* [graph_size];
        for (int i = 0; i < graph_size; i++) {
            nr_of_edges_sorted[i] = new int[2];
            nr_of_edges_sorted[i][0] = i;//vertex index
            nr_of_edges_sorted[i][1] = nr_of_edges[i];
        }
        //cout << endl;
        //PrintGraph(graph,nr_of_edges, graph_size);
        PrintDegreeSequence(nr_of_edges_sorted, graph_size);
        ComponentsAndBipartite(graph, nr_of_edges, graph_size, &nr_of_components);
        Eccentricity(graph, nr_of_edges, graph_size, nr_of_components);//cout << "?" << endl;//the eccentricity of vertices (within the components)
        cout << "?" << endl;//planarity;
        Greedy(graph, nr_of_edges, graph_size); //cout << "?" << endl;//greedy (the vertex order according to its number)
        LF(graph, nr_of_edges_sorted, nr_of_edges, graph_size);//cout << "?" << endl;//LF method (ties are solved by the vertex number)
        cout << "?" << endl;//SLF(graph, nr_of_edges_sorted, nr_of_edges, graph_size, nr_of_components);//cout << "?" << endl;//SLF method 
        cout << "?" << endl;//C4(graph, nr_of_edges, graph_size);//the number of different C4 subgraphs
        ComplementEdges(nr_of_edges, graph_size);//cout << "?" << endl;//the number of the graph complement's edges
        //cout << endl;
        for (int i = 0; i < graph_size; i++) {
            delete[] graph[i];
            delete[] nr_of_edges_sorted[i];
        }
        delete[] graph;
        delete[] nr_of_edges;
        delete[] nr_of_edges_sorted;
    }

    return 0;
}

