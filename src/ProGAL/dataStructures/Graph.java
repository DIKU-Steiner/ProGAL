package ProGAL.dataStructures;

import java.util.ArrayList;
import java.util.Stack;

//TODO: This class is copyrighted by the authors. We need to check if its ok to use this! 

/**
 *  The <tt>Graph</tt> class represents a graph of vertices named 0 through <em>V</em> - 1.
 *  It supports the following two primary operations: add an edge to the graph,
 *  iterate over all of the vertices adjacent to a vertex. It also provides
 *  methods for returning the number of vertices <em>V</em> and the number
 *  of edges <em>E</em>. Parallel edges and self-loops are permitted.
 *  <p>
 *  This implementation uses an adjacency-lists representation, which 
 *  is a vertex-indexed array of {@link Bag} objects.
 *  All operations take constant time (in the worst case) except
 *  iterating over the vertices adjacent to a given vertex, which takes
 *  time proportional to the number of such vertices.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/41undirected">Section 4.1</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class Graph {
    private int V;             // number of vertices
    private int E;             // number of edges 
    private boolean directed;  // TRUE if the graph is directed
    private ArrayList<ArrayList<Integer>> adj = new ArrayList<ArrayList<Integer>>();   // list of adjacency lists
    
    /**
     * Initializes an empty graph with <tt>V</tt> vertices and 0 edges.
     * @param V the number of vertices
     * @throws java.lang.IllegalArgumentException if <tt>V</tt> < 0
     */
    public Graph(int V, boolean directed) {
        if (V < 0) throw new IllegalArgumentException("Number of vertices must be nonnegative");
        this.V = V;
        this.E = 0;
        this.directed = directed;
        for (int v = 0; v < V; v++) adj.add(new ArrayList<Integer>());
    }

    /**
     * Initializes a new graph that is a deep copy of <tt>G</tt>.
     * @param G the graph to copy
     */
    public Graph(Graph G) {
        this(G.V(), G.directed);
        this.E = G.E();
        for (int v = 0; v < G.V(); v++) {
            // reverse so that adjacency list is in same order as original
            Stack<Integer> reverse = new Stack<Integer>();
            for (int w : G.adj.get(v)) {
                reverse.push(w);
            }
            for (int w : reverse) {
                adj.get(v).add(w);
            }
        }
    }

    /**
     * Returns the number of vertices in the graph.
     * @return the number of vertices in the graph
     */
    public int V() { return V; }

    /**
     * Returns the number of edges in the graph.
     * @return the number of edges in the graph
     */
    public int E() { return E; }

    /**
     * Adds the edge v-w to the graph.
     * @param v start vertex of the edge
     * @param w end vertex of the edge
     * @throws java.lang.IndexOutOfBoundsException unless both 0 <= v < V and 0 <= w < V
     */
    public void addEdge(int v, int w) {
        if (v < 0 || v >= V) throw new IndexOutOfBoundsException();
        if (w < 0 || w >= V) throw new IndexOutOfBoundsException();
        E++;
        adj.get(v).add(w);
        if (!directed) adj.get(w).add(v);
    }

    public void removeEdge(int v, int w){
    	if (v < 0 || v >= V) throw new IndexOutOfBoundsException();
        if (w < 0 || w >= V) throw new IndexOutOfBoundsException();
        boolean removed = adj.get(v).remove((Integer)w);
        if (!directed) adj.get(w).remove((Integer)v);
        if(removed) E--;
    }
    
    /**
     * Returns TRUE if there is an edge from vertex <tt>v</tt> to vertex <tt>w</tt>
     * @param v 
     * @param w
     * @return TRUE if there is an edge from vertex <tt>v</tt> to vertex <tt>w</tt>
     */
    public boolean adjacent(int v, int w) {
        if (v < 0 || v >= V) throw new IndexOutOfBoundsException();
        if (w < 0 || w >= V) throw new IndexOutOfBoundsException();
    	return adj.get(v).contains(w);
    }

    /**
     * Returns the vertices adjacent to vertex <tt>v</tt>.
     * @return the vertices adjacent to vertex <tt>v</tt> as an Iterable
     * @param v the vertex
     * @throws java.lang.IndexOutOfBoundsException unless 0 <= v < V
     */
    public Iterable<Integer> adj(int v) {
        if (v < 0 || v >= V) throw new IndexOutOfBoundsException();
        return adj.get(v);
    }

    /** Returns list of vertices at most depth edges away from vertex q. Vertex q is not in the returned list */
    public ArrayList<Integer> breadthFirst(int q, int depth) {
    	ArrayList<Integer> vertices = new ArrayList<Integer>();
    	vertices.add(q);
    	int level = 0;
    	int u;
    	int lastInLevel = vertices.size();
    	int i = 0;
    	while (level < depth) {
    		level++;
    		lastInLevel = vertices.size();
    		while (i < lastInLevel) {
    			u = vertices.get(i);
    			for (int v : adj(u)) 
    				if (!vertices.contains(v)) vertices.add(v);
    			i++;
    		}
    	}
    	vertices.remove(0);
    	return vertices;
    }
    
    /** returns list vertices at most depth edges away from a subset of vertices.*/
    public ArrayList<Integer> breadthFirst(ArrayList<Integer> qVertices, int depth) {
    	ArrayList<Integer> vertices = new ArrayList<Integer>();
    	ArrayList<Integer> bVertices = new ArrayList<Integer>();
    	for (int q : qVertices) {
    		bVertices = breadthFirst(q, depth);
    		for (int b : bVertices) 
    			if (!vertices.contains(b)) vertices.add(b);
    	}
    	return vertices;
    }
    
    /**
     * Returns a string representation of the graph.
     * This method takes time proportional to <em>E</em> + <em>V</em>.
     * @return the number of vertices <em>V</em>, followed by the number of edges <em>E</em>,
     *    followed by the <em>V</em> adjacency lists
     */
    public String toString() {
        StringBuilder s = new StringBuilder();
        String NEWLINE = System.getProperty("line.separator");
        s.append(V + " vertices, " + E + " edges " + NEWLINE);
        for (int v = 0; v < V; v++) {
            s.append(v + ": ");
            for (int w : adj.get(v)) {
                s.append(w + " ");
            }
            s.append(NEWLINE);
        }
        return s.toString();
    }


    
    /**
     * Unit tests the <tt>Graph</tt> data type.
     */
    public static void main(String[] args) {
        Graph G = new Graph(4, false);
        G.addEdge(0, 3);
        G.addEdge(2, 1);
    }
}

