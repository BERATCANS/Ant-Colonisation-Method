import java.util.ArrayList;
import java.util.List;
public class Ant {
    private List<Integer> tour;
    private boolean[] visited;

    public Ant(int numNodes) {
        tour = new ArrayList<>();
        visited = new boolean[numNodes];
        for (int i = 0; i < numNodes; i++) {
            visited[i] = false;
        }
    }
    public List<Integer> getTour() {
        return tour;
    }

    public void visitCity(int node) {
        tour.add(node);
        visited[node] = true;
    }
    public boolean visited(int node) {
        return visited[node];
    }
}