package ProGAL.dataStructures;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class Kruskal {
	private final List<Edge> mstEdges = new ArrayList<Edge>();
	private double length = 0.0; 
	
	public Kruskal(Collection<? extends Vertex> V, Collection<? extends Edge> E){
		List<Edge> sortedEdges = new ArrayList<Edge>(E);
		Collections.sort(sortedEdges, new Comparator<Edge>(){
			public int compare(Edge e0, Edge e1) {
				return Double.compare(e0.getDist(),e1.getDist());
			}});
		
		UnionFind<Vertex> ds = new UnionFind<Vertex>();
		for(Edge e: sortedEdges){
			if(  ds.find(e.getV1())!=ds.find(e.getV2())  ){
				length+=e.getDist();
				mstEdges.add(e);
				ds.union(e.getV1(), e.getV2());
			}
		}
	}
	
	public List<Edge> getEdges(){ return new ArrayList<Edge>(mstEdges); }
	public double getLength(){ return length; }

	public static interface Vertex{}
	
	public static interface Edge{
		public Vertex getV1();
		public Vertex getV2();
		public double getDist();
	}
}
