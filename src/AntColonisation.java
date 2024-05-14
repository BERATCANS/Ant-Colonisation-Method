/**
 * @author Beratcan Dogan, Student ID: 2021400132
 * @date 13.05.2024
 *this code find the shortest path with two methods one is ant colony other is brute-force.
 */
import java.awt.*;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;


public class AntColonisation {
    /**
     * Main method come together all my methods and my project starts from here.
     * @param args It is necessary to run code.
     * @throws FileNotFoundException If my project cannot find files it gives this error.
     */
    public static void main(String[] args) throws FileNotFoundException {
        double startTime = System.currentTimeMillis(); //to calc time
        ArrayList<Node> nodeList = readNodesFromFile("input05.txt"); //calls file reading function and keep nodes in nodeList
        int chosenMethod = 2;
        if(chosenMethod==2){//ant colony optimisation method
            int antColonyOutput = 1; //to draw two types of stdDraw I defined variable
            int N = 300; //number of iterations
            int M = 150; //number of ants
            double pheromoneIntensity = 0.1; //initial intensity
            double alpha = 0.9;
            double beta = 1.5;
            double degradationFactor = 0.9; //to decrease intensity in every iteration
            double Q = 0.0001; //to increase intensity in every ant route
            double[][] pheromones = new double[nodeList.size()][nodeList.size()]; // to use stdAntColony function I need to define pheromones in main

            List<Integer> shortestPath = antColony(nodeList,pheromones, N, M,pheromoneIntensity,alpha,beta,degradationFactor,Q); // it calls antColony method I have to give all variables because of writing variables beginning of main.
            Integer[] shortestRoute = objectToArray(shortestPath); //to faster the brute-force method I used array but in ant colony I used array object so ı need to convert to use calc distance
            double shortestDistance = calculateRouteDistance(shortestRoute, nodeList);//calls function which helps to calculate route distance

            System.out.println("Method: Ant Colony Optimisation Method");
            consoleWriting(shortestDistance,shortestRoute,startTime); // to avoid 2 times I call consoleWriting method
            if (antColonyOutput == 1){
                stdPart(nodeList, shortestRoute); //first type of drawing
            }
            else if(antColonyOutput == 2){
                stdAntColony(nodeList, pheromones); //second type of drawing
            }
        }
        else if (chosenMethod == 1){ //brute-force Method
            double[] shortest = {Double.MAX_VALUE}; //to avoid not changing in function. I define double[]
            Integer[] shortestRoute = new Integer[nodeList.size()-1];
            Integer[] indexes = new Integer[nodeList.size()-1];
            for (int l = 0; l < indexes.length; l++) {
                indexes[l] = l+1;
            }
            permute(shortestRoute, nodeList, indexes,shortest, 0); //calls permute function and changing shortestRoute
            double shortestDistance = shortest[0];
            System.out.println("Method: Brute-Force Method");
            consoleWriting(shortestDistance,shortestRoute,startTime);
            stdPart(nodeList, shortestRoute);
        }


    }
    /**
     * Executes the Ant Colony Optimization algorithm to find the shortest path between nodes.
     *
     * @param nodes The list of nodes representing locations.
     * @param pheromones The pheromone intensity matrix.
     * @param numAnts The number of ants to be used in the algorithm.
     * @param numIterations The number of iterations the algorithm will run for.
     * @param pheromoneIntensity The initial intensity of pheromones.
     * @param alpha The pheromone influence parameter.
     * @param beta The heuristic information influence parameter.
     * @param degradationFactor The rate at which pheromones degrade over time.
     * @param Q The amount of pheromone deposited by each ant.
     * @return The shortest path found by the Ant Colony Optimization algorithm.
     */
    public static List<Integer> antColony(List<Node> nodes,double[][] pheromones, int numAnts, int numIterations, double pheromoneIntensity, double alpha, double beta, double degradationFactor, double Q) {
        int numNodes = nodes.size();
        for (int i = 0; i < numNodes; i++) {
            for (int j = 0; j < numNodes; j++) {
                pheromones[i][j] = pheromoneIntensity; //gives initial intensity value.
            }
        }

        List<Integer> shortestPath = null;
        double shortestDistance = Double.MAX_VALUE;
        for (int iter = 0; iter < numIterations; iter++) {
            List<Ant> ants = createAnts(numNodes, numAnts);
            for (Ant ant : ants) {
                for (int i = 0; i < numNodes; i++) { //executes each ant according to the probabilities arising from the edge value
                    int nextCity = selectNextNode(ant, nodes, pheromones,alpha,beta);
                    ant.visitCity(nextCity);
                }
                updatePheromones(pheromones, ant, nodes, Q);
            }

            for (Ant ant : ants) {
                List<Integer> tour = ant.getTour();
                Integer[] route = objectToArray(tour);
                double tourDistance = calculateRouteDistance(route, nodes);

                if (tourDistance < shortestDistance) { // gets the shortest path and distance every iteration.
                    shortestDistance = tourDistance;
                    shortestPath = new ArrayList<>(tour);
                }
            }
            for (int i = 0; i < pheromones.length; i++) { //multiply with degradation factor after every iteration.
                for (int j = 0; j < pheromones[i].length; j++) {
                    pheromones[i][j] *= (degradationFactor);
                }
            }
        }
        return shortestPath;
    }
    /**
     * Creates a list of ants for the Ant Colony Optimization algorithm.
     *
     * @param numNodes The number of nodes or cities in the problem.
     * @param numAnts The number of ants to be created.
     * @return A list of ants initialized with the specified number of nodes.
     */
    public static List<Ant> createAnts(int numNodes, int numAnts) {
        List<Ant> ants = new ArrayList<>();
        for (int i = 0; i < numAnts; i++) {
            ants.add(new Ant(numNodes));
        }
        return ants;
    }
    /**
     * Selects the next node to be visited by an ant based on pheromone levels and heuristic information.
     *
     * @param ant The ant currently traversing the graph.
     * @param nodes The list of nodes representing cities or locations.
     * @param pheromones The pheromone intensity matrix.
     * @param alpha The pheromone influence parameter.
     * @param beta The heuristic information influence parameter.
     * @return The index of the next node to be visited.
     */
    public static int selectNextNode(Ant ant, List<Node> nodes, double[][] pheromones,double alpha,double beta) {
        List<Integer> tour = ant.getTour();
        if (tour.isEmpty()) {
            return 0;
        }
        else {
            int currentNode = tour.get(tour.size() - 1);
            double[] probabilities = calculateProbabilities(currentNode, ant, nodes, pheromones,alpha,beta); //gets probabality of each unvisited node
            double rand = Math.random();
            double tempSum = 0;
            for (int i = 0; i < probabilities.length; i++) {  // choose randomly node
                tempSum += probabilities[i];
                if (rand < tempSum) {
                    return i;
                }
            }
            return -1; //it gives missing return statement.If I do not write this.
        }
    }
    /**
     * Calculates the probabilities of selecting each node as the next node to be visited by an ant.
     *
     * @param currentNode The index of the current node.
     * @param ant The ant currently traversing the graph.
     * @param nodes The list of nodes representing cities or locations.
     * @param pheromones The pheromone intensity matrix.
     * @param alpha The pheromone influence parameter.
     * @param beta The heuristic information influence parameter.
     * @return An array of probabilities for selecting each node.
     */
    public static double[] calculateProbabilities(int currentNode, Ant ant, List<Node> nodes, double[][] pheromones,double alpha,double beta) {
        int numNodes = nodes.size();
        double[] probabilities = new double[numNodes];
        double total = 0;

        for (int i = 0; i < numNodes; i++) {
            if (!ant.visited(i)) {
                double pheromone = Math.pow(pheromones[currentNode][i], alpha);
                double distance = Math.pow(1.0 / nodes.get(currentNode).distanceTo(nodes.get(i)), beta);
                probabilities[i] = pheromone * distance;  //gives probability to each node by pheromone and distance values
                total += probabilities[i];
            }
        }
        for (int i = 0; i < numNodes; i++) {
            probabilities[i] /= total; //gives a value 0 to 1
        }
        return probabilities;
    }
    /**
     * Updates the pheromone levels on the edges of the graph based on the paths taken by ants.
     *
     * @param pheromones The pheromone intensity matrix.
     * @param ant ant that have traversed the graph.
     * @param nodes The list of nodes representing cities or locations.
     * @param Q The amount of pheromone deposited by each ant.
     */
    public static void updatePheromones(double[][] pheromones, Ant ant, List<Node> nodes,double Q) {
        List<Integer> tour = ant.getTour();
        Integer[] route = objectToArray(tour);
        double pheromoneDelta = Q / calculateRouteDistance(route, nodes);  //increase pheromone in each ant walk
        for (int i = 0; i < tour.size() - 1; i++) {
            int node1 = tour.get(i);
            int node2 = tour.get(i + 1);
            pheromones[node1][node2] += pheromoneDelta;
            pheromones[node2][node1] += pheromoneDelta;
            }
        }
    /**
     * Permutes through all possible routes to find the shortest path using the brute-force method.
     *
     * @param shortestRoute The array representing the shortest route found so far.
     * @param nodes The list of nodes representing cities or locations.
     * @param arr The array to be permuted.
     * @param currentDistance The current shortest distance found.
     * @param k The current index of the array being permuted.
     */
    public static void permute(Integer[] shortestRoute, ArrayList<Node> nodes, Integer[] arr, double[] currentDistance, int k) {
        if (k == arr.length) {
            double distance = calculateRouteDistance(arr, nodes);
            if (distance < currentDistance[0]) {
                currentDistance[0] = distance;
                System.arraycopy(arr, 0, shortestRoute, 0, arr.length);
            }
        }
        else {
            for (int i = k; i < arr.length; i++) {
                Integer temp = arr[i];
                arr[i] = arr[k];
                arr[k] = temp;
                permute(shortestRoute, nodes, arr, currentDistance, k + 1);
                temp = arr[k];
                arr[k] = arr[i];
                arr[i] = temp;
            }
        }
    }
    /**
     * Calculates the total distance of a route based on the Euclidean distance between nodes.
     *
     * @param route The array representing the route.
     * @param nodes The list of nodes representing cities or locations.
     * @return The total distance of the route.
     */
    public static double calculateRouteDistance(Integer[] route, List<Node> nodes) {
        double distance = 0;
        int prevIdx = 0;
        for (int nodeIdx : route) {
            distance += Math.sqrt(Math.pow(nodes.get(nodeIdx).getX() - nodes.get(prevIdx).getX(), 2) +
                    Math.pow(nodes.get(nodeIdx).getY() - nodes.get(prevIdx).getY(), 2));
            prevIdx = nodeIdx;
        }
        distance += Math.sqrt(Math.pow(nodes.getFirst().getX() - nodes.get(prevIdx).getX(), 2) +
                Math.pow(nodes.getFirst().getY() - nodes.get(prevIdx).getY(), 2));
        return distance;
    }
    /**
     * Reads the node coordinates from a file and creates a list of nodes.
     *
     * @param fileName The name of the file containing the node coordinates.
     * @return The list of nodes representing locations.
     * @throws FileNotFoundException If the specified input file cannot be found.
     */
    public static ArrayList<Node> readNodesFromFile(String fileName) throws FileNotFoundException {
        File file = new File(fileName);
        Scanner scanner = new Scanner(file);
        int i = 0;
        ArrayList<Node> nodeList = new ArrayList<>();
        while (scanner.hasNextLine()) {
            i++;
            String line = scanner.nextLine();
            String[] parts = line.split(",");
            double x = Double.parseDouble(parts[0]);
            double y = Double.parseDouble(parts[1]);
            Node node = new Node(x, y, i);
            nodeList.add(node);
        }
        scanner.close();
        return nodeList;
    }
    /**
     * Prints the shortest path and the time taken to find it to the console.
     *
     * @param shortestDistance The shortest distance found.
     * @param shortestRoute The shortest route found.
     * @param startTime The start time of the algorithm execution to calculate time.
     */
    public static void consoleWriting(double shortestDistance,Integer[] shortestRoute,double startTime) {
        System.out.printf("Shortest Distance: %.5f \n", shortestDistance);
        System.out.print("Shortest path: [1, ");
        for (int routeOrder : shortestRoute) {
            System.out.print(routeOrder + 1 + ", ");
        }
        System.out.println("1]");
        double endTime = System.currentTimeMillis();
        double time = (endTime - startTime) / 1000;
        System.out.printf("Time it takes to find the shortest path: %.2f seconds.", time);
    }
    /**
     * Draws the nodes and the shortest path on a canvas for visualization.
     *
     * @param nodeList The list of nodes representing locations.
     * @param shortestRoute The shortest route found.
     */
    public static void stdPart(ArrayList<Node> nodeList,Integer[] shortestRoute){
        StdDraw.enableDoubleBuffering();
        int canvasWidth = 800;
        int canvasHeight = 800;
        StdDraw.setCanvasSize(canvasWidth, canvasHeight);
        StdDraw.setXscale(0, 1);
        StdDraw.setYscale(0, 1);
        double circleRadius = 0.02;
        StdDraw.clear(StdDraw.WHITE);
        StdDraw.setPenColor(StdDraw.BLACK);
        StdDraw.setFont(new Font("Serif", Font.BOLD, 14));
        StdDraw.setPenRadius(0.003);
        int prevIdx = 0;
        for(int n: shortestRoute){
            StdDraw.line(nodeList.get(prevIdx).getX(), nodeList.get(prevIdx).getY(), nodeList.get(n).getX(), nodeList.get(n).getY());
            prevIdx = n;
        }
        StdDraw.line(nodeList.get(prevIdx).getX(), nodeList.get(prevIdx).getY(), nodeList.getFirst().getX(), nodeList.getFirst().getY());

        for(Node node:nodeList){
            if (node.getOrder()==1) {
                StdDraw.setPenColor(StdDraw.PRINCETON_ORANGE);
            } else {
                StdDraw.setPenColor(StdDraw.LIGHT_GRAY);
            }
            StdDraw.filledCircle(node.getX(), node.getY(), circleRadius);
            StdDraw.setPenColor(StdDraw.BLACK);
            StdDraw.text(node.getX(), node.getY(), String.valueOf(node.getOrder()));
        }
        StdDraw.show();
    }
    /**
     * Draws the nodes and the pheromone trails on a canvas for visualization in the Ant Colony Optimization method.
     *
     * @param nodeList The list of nodes representing locations.
     * @param pheromones The pheromone intensity matrix.
     */
    public static void stdAntColony(ArrayList<Node> nodeList, double[][] pheromones) {
        StdDraw.enableDoubleBuffering();
        int canvasWidth = 800;
        int canvasHeight = 800;
        StdDraw.setCanvasSize(canvasWidth, canvasHeight);
        StdDraw.setXscale(0, 1);
        StdDraw.setYscale(0, 1);
        double circleRadius = 0.02;
        StdDraw.clear(StdDraw.WHITE);
        StdDraw.setPenColor(StdDraw.BLACK);
        StdDraw.setFont(new Font("Serif", Font.BOLD, 14));
        double maxValue = 0;
        for (int i = 0; i < pheromones.length; i++) {
            double sum = 0;
            for (int j = 0; j < pheromones.length; j++) {
                sum += pheromones[i][j];
                if (maxValue < pheromones[i][j]) {
                    maxValue = pheromones[i][j]; // get max value to standardize line thickness
                }
            }
            double average = sum / pheromones.length/2; // to do not draw every line
            double penRadius = 0.01 / maxValue;
            for (int j = 0; j < pheromones.length; j++) {
                if (pheromones[i][j] > average) {
                    StdDraw.setPenRadius(penRadius * pheromones[i][j]);
                    StdDraw.line(nodeList.get(i).getX(), nodeList.get(i).getY(), nodeList.get(j).getX(), nodeList.get(j).getY());
                }
            }

            for (Node node : nodeList) {
                StdDraw.setPenColor(StdDraw.LIGHT_GRAY);
                StdDraw.filledCircle(node.getX(), node.getY(), circleRadius);
                StdDraw.setPenColor(StdDraw.BLACK);
                StdDraw.text(node.getX(), node.getY(), String.valueOf(node.getOrder()));
            }
            StdDraw.show();
        }
    }

    /**
     * converts List<Integer> to Integer[] to use same distance calculation method with brute-force.
     * @param tour route of ant in array list.
     * @return route of ant in array.
     */
    public static Integer[] objectToArray(List<Integer> tour){
        Integer[] route = new Integer[tour.size()-1];
        for (int i=0;i<tour.size()-1;i++){
            route[i] = tour.get(i+1);
        }
        return route;
    }

}


