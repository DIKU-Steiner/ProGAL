package ProGAL.dataStructures;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ProGAL.dataStructures.Pair;

/**
 * A weighted graph implemented using adjacency lists. At construction it can be specified if the graph should be 
 * directed or undirected (standard is undirected). Each vertex is represented by an index (does not need to 
 * be contigous). Addition and removal of vertices and edges is supported.
 * 
 * If a weight-function is not specified the weight of each edge is assumed to be 1.
 *   
 * @author R.Fonseca
 */
public class WeightedGraph {
	private final Set<Integer> vertices = new HashSet<Integer>();
	private final Map<Integer, List<Edge>> adjacencies = new HashMap<Integer, List<Edge>>();
	public final boolean directed;
	public WeightFunction weightFunction;
	
	public WeightedGraph(int vertices){
		this();
		
		for(int i=0;i<vertices;i++)
			addVertex(i);
	}
	
	public WeightedGraph(){
		this(false);
	}
	
	public WeightedGraph(boolean directed){
		this(directed, null);
	}
	
	public WeightedGraph(boolean directed, WeightFunction wFunc){
		this.directed = directed;
		this.weightFunction = wFunc;
	}
	
	/**
	 * Add the vertex with specified index. If index already exist in the graph nothing happens. 
	 * @param index 
	 */
	public void addVertex(int index){
		if(!vertices.contains(index)){
			vertices.add(index);
			adjacencies.put(index, new LinkedList<Edge>());
		}
	}
	
	/**
	 * Add an edge between specified vertices. If either of the vertices is not in the graph nothing will 
	 * be done and <code>null</code> will be returned. Note that multiple edges can exist between a pair of 
	 * vertices.  
	 * @param u Index of first vertex
	 * @param v Index of second vertex
	 * @return the constructed edge
	 */
	public Edge addEdge(int u, int v){
		if( !(vertices.contains(u) && vertices.contains(v)) ) return null;
		
		Edge e = new Edge(u,v);
		adjacencies.get(u).add(e);
		if(!directed)
			adjacencies.get(v).add(e);
		
		return e;
	}
	
	public Edge getEdge(int u, int v){
		if( !(vertices.contains(u) && vertices.contains(v)) ) return null;

		for(Edge e: adjacencies.get(u))
			if( (e.fst==v && e.snd==u) || (e.snd==v && e.fst==u) ) return e;
		
		if(!directed)
			for(Edge e: adjacencies.get(v))
				if( (e.fst==v && e.snd==u) || (e.snd==v && e.fst==u) ) return e;
		
		return null;
	}
	
	/**
	 * Remove the specified edge
	 * @param edge
	 * @return true iff <code>edge</code> was successfully removed
	 */
	public boolean removeEdge(Edge edge){
		int u = edge.fst, v = edge.snd;
		if(!adjacencies.containsKey(u)) return false;
		boolean found = adjacencies.get(u).remove(edge);
		
		if(!directed && found)
			adjacencies.get(v).remove(edge);
		
		return found;
	}
	
	/**
	 * Remove the specified vertex. If there are edges connected they will be removed as well.
	 * @param index
	 * @return true iff the vertex previously existed in the graph
	 */
	public boolean removeVertex(int index){
		if(!vertices.contains(index)) return false;
		
		List<Edge> removed = adjacencies.remove(index);
		if(!directed){
			for(Edge e: removed){
				adjacencies.get(e.opposite(index)).remove(e);
			}
		}

		vertices.remove((Integer)index);
		return true;
	}
	
	public Set<Integer> getVertices(){
		return new HashSet<Integer>(vertices);
	}
	
	/**
	 * Get all edges adjacent to the specified vertex.
	 * @param v
	 * @return a list containing all edges adjacent to <code>v</code>
	 */
	public List<Edge> getAdjacentEdges(int v){
		return new ArrayList<Edge>(adjacencies.get(v));
	}
	
	/**
	 * Get all edges in this graph.
	 * @return a set of edges
	 */
	public Set<Edge> getEdges(){
		HashSet<Edge> edges = new HashSet<Edge>();
		for(Integer v: vertices){
			edges.addAll(adjacencies.get(v));
		}
		return edges;
	}
	
	/**
	 * Calculate the summed weight of all edges.
	 * @return a double indicating the sum of all edge-weights. Might be infinity.
	 */
	public double getWeight(){
		double ret = 0;
		for(Edge e: getEdges())
			ret+=e.weight();
		return ret;
	}
	
	/** Return the set of all vertices reachable from v */
	public Set<Integer> reachable(int v){
		Set<Integer> ret = new HashSet<Integer>();
		LinkedList<Integer> open = new LinkedList<Integer>();
		open.add(v);
		while(!open.isEmpty()){
			int visit = open.poll();
			ret.add(visit);
			
			for(Edge e: adjacencies.get(visit)){
				int o = e.opposite(visit);
				if(!ret.contains(o)) open.push(o);
			}
			
		}
		
		return ret;
	}
	
	public String toString(){
		StringBuilder sb = new StringBuilder();
		sb.append("Graph[\n");
		sb.append(" Vertices: ");
		sb.append(vertices.toString());
		sb.append('\n');
		sb.append(" Edges: ");
		sb.append(getEdges().toString());
		sb.append('\n');
		sb.append("]");
		return sb.toString();
	}
	
	public WeightedGraph clone(){
		WeightedGraph ret = new WeightedGraph(directed);
		ret.weightFunction = weightFunction;
		for(Integer v: vertices) ret.addVertex(v);
		for(Edge e: getEdges()) ret.addEdge(e.fst, e.snd);
		return ret;
	}
	
	public static WeightedGraph createFullyConnectedGraph(int vertices){
		WeightedGraph g = new WeightedGraph(vertices);
		for(int i=0;i<vertices;i++){
			for(int j=i+1;j<vertices;j++){
				g.addEdge(i, j);
			}
		}
		return g;
	}
	
	
	public class Edge extends Pair<Integer, Integer>{
		protected Edge(int v1, int v2){
			super(v1,v2);
		}
		
		public int opposite(int v){
			if(fst==v) return snd;
			if(snd==v) return fst;
			throw new RuntimeException("Edge "+this+" not connected to "+v);
		}
		
		public double weight(){
			if(weightFunction==null) return 1;
			return weightFunction.w(this);
		}
		
		public String toString(){
			return String.format("(%d,%d)", fst,snd);
		}
	}
	
	public static interface WeightFunction{
		public double w(Pair<Integer,Integer> p);
	}

}
