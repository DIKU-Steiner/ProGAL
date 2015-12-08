package ProGAL.dataStructures.shortestPath;

import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.PriorityQueue;

import ProGAL.dataStructures.WeightedGraph;
import ProGAL.dataStructures.WeightedGraph.Edge;

public class Dijkstra {
	private final Map<Integer, Double> dMap = new HashMap<Integer, Double>();
	private final Map<Integer, Integer> pMap = new HashMap<Integer, Integer>();
	
	public Dijkstra(WeightedGraph g, int source){
		initialize(g, source);
		PriorityQueue<Integer> Q = new PriorityQueue<Integer>(g.getVertices().size(), new Comparator<Integer>(){
			@Override
			public int compare(Integer o1, Integer o2) {
				return Double.compare(dMap.get(o1), dMap.get(o2));
			}});
		Q.addAll(g.getVertices());
		
		while(!Q.isEmpty()){
			int u = Q.poll();
			for( Edge e: g.getAdjacentEdges(u) ){
				if(relax(u,e.opposite(u), e.weight()))
					Q.add(e.opposite(u));//Instead of decreasing keys we just add the element to queue again.
			}
		}
	}
	
	private void initialize(WeightedGraph g, int source){
		for(Integer v: g.getVertices()){
			dMap.put(v, Double.POSITIVE_INFINITY);
			pMap.put(v, null);
		}
		dMap.put(source, 0.0);
	}
	
	private boolean relax(int u, int v, double w){
//		System.out.printf("relax(%d,%d,%.3f)\n",u,v,w);
		double duw = dMap.get(u)+w;
		double dv = dMap.get(v);
		if(dv>duw){
			dMap.put(v, duw);
			pMap.put(v, u);
			return true;
		}
		return false;
	}
	
	/**
	 * Return the distance from the source (specified at construction) to a certain vertex.
	 * @return the shortest path from source to v
	 */
	public double getDistance(int v){
		return dMap.get(v);
	}
	
	/**
	 * Return the shortest path from the source to the specified vertex. If no path exists <code>null</code> is 
	 * returned.
	 * @param v a vertex in the graph
	 * @return a linked list containing all vertices from (including) the source to (including) <code>v</code>
	 */
	public LinkedList<Integer> getShortestPath(int v){
		LinkedList<Integer> ret = new LinkedList<Integer>();
		Integer p = v;
		while(p!=null){
			ret.addFirst(p);
			p = pMap.get(p);
		}
		if(dMap.get(v)==0.0) return ret; //If v is the source just return
		if(ret.size()<2) return null;
		return ret;
	}
	
	public static void main(String[] args) {
		WeightedGraph g = new WeightedGraph(false);
		for(int i=0;i<5;i++) g.addVertex(i);
		g.addEdge(0,1);
		g.addEdge(1,3);
		g.addEdge(0,4);
		g.addEdge(2,4);
		Dijkstra sp = new Dijkstra(g, 3);
		System.out.println(sp.getDistance(4));
		
//		WeightedGraph g = new WeightedGraph(true);
//		g.addVertex(0);
//		g.addVertex(1);
//		g.addVertex(2);
//		g.addVertex(3);
//		g.addVertex(4);
//		g.addVertex(5);
//		g.addEdge(0, 1);
//		g.addEdge(0, 3);
//		g.addEdge(1, 2);
//		g.addEdge(1, 3);
//		g.addEdge(2, 4);
//		g.addEdge(3, 1);
//		g.addEdge(3, 2);
//		g.addEdge(3, 4);
//		g.addEdge(4, 2);
//		g.addEdge(4, 0);
//		g.addEdge(4, 5);
//		g.weightFunction = new WeightedGraph.WeightFunction(){
//			@Override
//			public double w(Pair<Integer, Integer> p) {
//				if(p.fst==0 && p.snd==1) return 10;
//				if(p.fst==0 && p.snd==3) return 5;
//				if(p.fst==1 && p.snd==2) return 1;
//				if(p.fst==1 && p.snd==3) return 2;
//				if(p.fst==2 && p.snd==4) return 4;
//				if(p.fst==3 && p.snd==1) return 3;
//				if(p.fst==3 && p.snd==2) return 1;
//				if(p.fst==3 && p.snd==4) return 2;
//				if(p.fst==4 && p.snd==2) return 6;
//				if(p.fst==4 && p.snd==0) return 7;
//				if(p.fst==4 && p.snd==5) return 1;
//				return Double.POSITIVE_INFINITY;
//			}
//		};
//		Dijkstra sp = new Dijkstra(g, 5);
//		System.out.println(sp.getDistance(0));
	}

}
