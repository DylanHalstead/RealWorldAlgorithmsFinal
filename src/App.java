import java.lang.String;
import java.util.*;
import java.io.IOException;
import bridges.base.GraphAdjList;
import bridges.connect.Bridges;
import bridges.connect.DataSource;
import bridges.base.Edge;
import bridges.base.Element;
import bridges.validation.RateLimitException;
import bridges.data_src_dependent.City;

public class App {

    /*
     * Lab 5 main function, builds test graph, the MST of said graph, and MST of US
     * cities with populations based on popThresholds
     * 
     * @param args command line arguments (args[0] assignemnt id, args[1] username,
     * args[2] api key)
     * 
     * @throws IOException when failing reading the graph to visualize
     * 
     * @throws RateLimitException when exceeding the rate limit of the bridges api
     */
    public static void main(String[] args) throws IOException, RateLimitException {

        Bridges bridges = new Bridges(Integer.parseInt(args[0]), args[1], args[2]);
        GraphAdjList<String, String, Double> graph = new GraphAdjList<>();

        double edgeThickness = 3.0;
        double vertexSize = 10.0;
        bridges.setCoordSystemType("albersusa");
        bridges.setMapOverlay(true);
        int[] popThresholds = { 640_000 };
        for (int popThreshold : popThresholds) {
            DataSource ds = bridges.getDataSource();
            HashMap<String, String> params = new HashMap<String, String>();
            params.put("min_pop", String.valueOf(popThreshold));

            Vector<City> cities;
            try {
                cities = ds.getUSCitiesData(params);
            } catch (Exception e) {
                System.out.println("Error: " + e.getMessage());
                return;
            }

            graph = createCityGraph(cities);

            Map<Element<String>, Element<String>> childtoParent = fringePrimsAlgorithm(graph,
                    graph.getVertex("Charlotte_NC"));
            if (childtoParent.equals(null)) {
                System.out.println("Graph is not connected");
                return;
            }

            double totalCost = 0.;
            for (Element<String> vertex : childtoParent.keySet()) {
                Element<String> parent = childtoParent.get(vertex);
                if (parent != null) {
                    totalCost += graph.getEdgeData(vertex.getLabel(), parent.getLabel());
                }
            }

            if (childtoParent.size() > 200) {
                vertexSize = 3.0;
                edgeThickness = 1.0;
            } else {
                vertexSize = 10.0;
                edgeThickness = 3.0;
            }
            GraphAdjList<String, String, Double> MSTGraph = createPathGraph(graph, childtoParent, vertexSize,
                    edgeThickness);
            // Below needed to display larger graphs. Smaller than 40k has too many vertices
            // (LARGE_GRAPH_VERT_SIZE = 1000). When loading on website it crashes when type
            // is "graph-webgl" which "largegraph" is converted to; think this is a bridges
            // bug?
            MSTGraph.forceSmallVisualization(true);
            bridges.setDataStructure(MSTGraph);
            bridges.setTitle(
                    String.format("MST (US Cities) | Population Threshold: %d  | Vertices: %d | Edges: %d | Cost: %f",
                            popThreshold,
                            MSTGraph.getVertices().size(), childtoParent.size() - 1, totalCost));
            try {
                bridges.visualize();
            } catch (Exception e) {
                System.out.println(e);
            }
        }
    }

    /*
     * Prims algorithm (fringe vertices version) to find a minimum spanning tree of
     * a graph
     * 
     * @param graph the graph to find the MST starting from
     * 
     * @param startVtx the starting vertex
     * 
     * @return a map of child-to-parent relationships in the MST
     * 
     */
    static Map<Element<String>, Element<String>> fringePrimsAlgorithm(GraphAdjList<String, String, Double> graph,
            Element<String> startVtx) {
        Set<Element<String>> visited = new HashSet<>();
        Set<Element<String>> fringe = new HashSet<>();
        Map<Element<String>, Double> minDist = new HashMap<>();
        Map<Element<String>, Element<String>> childtoParent = new HashMap<>();

        for (Element<String> vertex : graph.getVertices().values()) {
            minDist.put(vertex, Double.MAX_VALUE);
        }

        minDist.put(startVtx, 0.);
        fringe.add(startVtx);
        childtoParent.put(startVtx, null);

        while (!fringe.isEmpty()) {
            Element<String> currentVertex = fringe.iterator().next();
            for (Element<String> vertex : fringe) {
                if (minDist.get(vertex) < minDist.get(currentVertex)) {
                    currentVertex = vertex;
                }
            }
            fringe.remove(currentVertex);

            if (minDist.get(currentVertex) == Double.MAX_VALUE) {
                System.out.println("Graph is not connected");
                return null;
            }

            visited.add(currentVertex);

            for (Edge<String, Double> edge : graph.outgoingEdgeSetOf(currentVertex.getLabel())) {
                Element<String> neighbor = graph.getVertex(edge.getTo());
                if (!visited.contains(neighbor) && minDist.get(neighbor) > edge.getEdgeData()) {
                    minDist.put(neighbor, edge.getEdgeData());
                    childtoParent.put(neighbor, currentVertex);
                    fringe.add(neighbor);
                }
            }
        }

        return childtoParent;
    }

    /*
     * Visualize the path of a graph using child-to-parent relationships; create a
     * new graph with only the edges in the path
     * 
     * @param originalGraph the graph with all vertices and edges
     * 
     * @param childtoParent a map of child-to-parent relationships in the search
     * 
     * @param vertexSize the size of the vertices in the path
     * 
     * @param edgeThickness the thickness of the edges in the path
     * 
     * @return the graph with the path visualized; removing all other edges
     */
    static GraphAdjList<String, String, Double> createPathGraph(GraphAdjList<String, String, Double> originalGraph,
            Map<Element<String>, Element<String>> childtoParent, double vertexSize, double edgeThickness) {
        GraphAdjList<String, String, Double> pathGraph = new GraphAdjList<>();
        for (Element<String> vertex : childtoParent.keySet()) {
            pathGraph.addVertex(vertex.getLabel(), vertex.getValue());
            Element<String> pathVtx = pathGraph.getVertex(vertex.getLabel());
            pathVtx.setLocation(vertex.getLocationX(), vertex.getLocationY());
            pathVtx.setColor("red");
            pathVtx.setSize(vertexSize);
            Element<String> parent = childtoParent.get(vertex);
            if (parent != null) {
                Element<String> pathParentVtx;
                if (!pathGraph.getVertices().containsKey(parent.getLabel())) {
                    pathGraph.addVertex(parent.getLabel(), parent.getValue());
                    pathParentVtx = pathGraph.getVertex(parent.getLabel());
                    pathParentVtx.setLocation(parent.getLocationX(), parent.getLocationY());
                    pathParentVtx.setColor("red");
                    pathParentVtx.setSize(vertexSize);
                }
                pathParentVtx = pathGraph.getVertex(parent.getLabel());
                pathGraph.addEdge(vertex.getLabel(), parent.getLabel(),
                        originalGraph.getEdgeData(vertex.getLabel(), parent.getLabel()));
                pathGraph.getLinkVisualizer(pathVtx.getLabel(), pathParentVtx.getLabel()).setColor("red");
                pathGraph.getLinkVisualizer(pathVtx.getLabel(), pathParentVtx.getLabel()).setThickness(edgeThickness);
                pathGraph.getLinkVisualizer(pathVtx.getLabel(), pathParentVtx.getLabel()).setLabel(
                        String.valueOf(originalGraph.getEdgeData(vertex.getLabel(), parent.getLabel())));
            }
        }
        return pathGraph;
    }

    /*
     * Create the city graphs of the US cities with population > 100000. One full
     * graph and one of the graph's MST
     * 
     * Note that when dealing with a large set of cities, creating the graph eats up
     * a large amount of ram and may crash with out of memory (ran into with min_pop
     * less than 40k). Fixed by manualling giving java more ram
     * "java -Xmx8192m <filename>" gives the jvm 8gb to work with
     * 
     * @param cities the cities to build the graph from
     * 
     * @return the graph with the cities as vertices and edges to all other
     * cities (complete graph); have weight as distance between them
     */
    static GraphAdjList<String, String, Double> createCityGraph(Vector<City> cities) {
        GraphAdjList<String, String, Double> graph = new GraphAdjList<>();
        for (City city : cities) {
            String cityLabel = String.format("%s_%s", city.getCity(), city.getState());
            graph.addVertex(cityLabel, cityLabel);
            graph.getVertex(cityLabel).setLocation(city.getLongitude(), city.getLatitude());
        }

        for (City city1 : cities) {
            String city1Label = String.format("%s_%s", city1.getCity(), city1.getState());
            for (City city2 : cities) {
                String city2Label = String.format("%s_%s", city2.getCity(), city2.getState());
                if (!city1.equals(city2)) {
                    double dist = getDist(city1.getLatitude(), city1.getLongitude(),
                            city2.getLatitude(), city2.getLongitude());
                    graph.addEdge(city1Label, city2Label, dist);
                }
            }
        }

        // set edge labels
        for (String v : graph.getVertices().keySet()) {
            for (Edge<String, Double> edge : graph.outgoingEdgeSetOf(v)) {
                String src = edge.getFrom(), dest = edge.getTo();
                String l = String.valueOf(graph.getEdgeData(src, dest));
                graph.getLinkVisualizer(src, dest).setLabel(l);
            }
        }

        return graph;
    }

    /*
     * Calculate the distance, in meters, between two points based on their latitude
     * and longitude
     * 
     * @param lat1 the latitude of the first point
     * 
     * @param long1 the longitude of the first point
     * 
     * @param lat2 the latitude of the second point
     * 
     * @param long2 the longitude of the second point
     * 
     * @return the distance between the two points in meters
     */
    static double getDist(double lat1, double long1, double lat2, double long2) {
        // uses the haversine formula
        final int R = 6371000; // meters
        final double phi1 = Math.toRadians(lat1);
        final double phi2 = Math.toRadians(lat2);
        final double delPhi = Math.toRadians((lat2 - lat1));
        final double delLambda = Math.toRadians((long2 - long1));

        final double a = Math.sin(delPhi / 2) * Math.sin(delPhi / 2)
                + Math.cos(phi1) * Math.cos(phi2)
                        * Math.sin(delLambda / 2) * Math.sin(delLambda / 2);
        final double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
        return R * c; // meters
    }
}
