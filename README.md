# Ant-Colonisation-Method
This Java program, named "AntColonisation", implements two algorithms for finding the shortest path between nodes in a graph representing locations: the Ant Colony Optimization method and the Brute-Force method.

Ant Colony Optimization Method:  

-This method utilizes a heuristic based on the pheromone trails left by artificial ants traversing the graph.  
-The ants probabilistically choose the next node to visit based on both the intensity of the pheromone trail and a heuristic function representing the desirability of each potential node.  
-The pheromone levels are updated based on the paths taken by the ants, with the aim of reinforcing the paths of shorter distances over time.   
-Parameters such as the number of iterations, the number of ants, pheromone intensity, and degradation factor are customizable to fine-tune the optimization process.  
-The program visualizes the nodes and the pheromone trails on a canvas for visualization.  

Brute-Force Method:
-This method exhaustively searches through all possible routes to find the shortest path by calculating the distance for each permutation.  
-It systematically generates permutations of node indices and calculates the total distance for each permutation.  
-The shortest path and its distance are updated whenever a shorter route is found during the search process.  
-The program also reads node coordinates from a file, performs necessary computations, and visualizes the results using the StdDraw library.  

Here is the table comparing the execution times and the results of the shortest paths found by the two codes:
<img width="926" alt="Screenshot 2024-05-13 at 05 08 48" src="https://github.com/BERATCANS/Ant-Colonisation-Method/assets/101462943/f2733975-96e3-4f06-839a-22a89f563700">  


Here is the samples of outputs of code: <img width="912" alt="ss of 4 ant-colony 05 36 23" src="https://github.com/BERATCANS/Ant-Colonisation-Method/assets/101462943/ce5d8f4c-042a-43a5-8cf5-b2d4f7ebb111">
<img width="912" alt="ss of 5 ant-colony2" src="https://github.com/BERATCANS/Ant-Colonisation-Method/assets/101462943/8a424346-2799-4a81-b9cd-4b8b8158ffb5">
