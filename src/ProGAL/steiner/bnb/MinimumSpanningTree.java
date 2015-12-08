package ProGAL.steiner.bnb;

import java.util.ArrayList;
import java.util.List;

import ProGAL.dataStructures.Kruskal;
import ProGAL.geom2d.delaunay.DTWithBigPoints;
import ProGAL.geom2d.delaunay.DelaunayTriangulationPawel;
import ProGAL.geomNd.Point;

public class MinimumSpanningTree{
	private final List<Vertex> vertices;
	private final List<Edge> edges;
	private final double length;
	
	public MinimumSpanningTree(Point[] sites) {
		vertices = createVertices(sites);
		edges = createEdges(vertices);
		Kruskal mst = new Kruskal(vertices,edges);
		List<Kruskal.Edge> edgeList = mst.getEdges();

		edges.clear();
		for(Kruskal.Edge ke: edgeList){
			edges.add((Edge)ke);
		}
		
		length = mst.getLength();
	}
	
	public double getLength(){ return length; }
	
	public int[][] getEdges(){
		int[][] ret = new int[edges.size()][2];
		for(int e=0;e<edges.size();e++){
			Edge edge = edges.get(e);
			ret[e][0] = edge.v1.idx;
			ret[e][1] = edge.v2.idx;
		}
		return ret;
	}

	private static List<Vertex> createVertices(Point[] sites){
		List<Vertex> ret = new ArrayList<Vertex>();
		for(int i=0;i<sites.length; i++) 
			ret.add(new Vertex(sites[i], i));
		return ret;
	}
	private static List<Edge> createEdges(List<Vertex> vertices){
		List<Edge> ret = new ArrayList<Edge>(); 
		for(Vertex v1: vertices){
			for(Vertex v2: vertices){
				if(v1==v2) continue;
				ret.add(   new Edge( v1, v2, v1.p.distance(v2.p) )   );
			}
		}
		return ret;
	}
	private static class Vertex implements Kruskal.Vertex{
		final Point p;
		final int idx;
		Vertex(Point p, int idx){ this.p = p; this.idx = idx; }
	}
	
	private static class Edge implements Kruskal.Edge{
		private Vertex v1,v2;
		private double dist;
		public Edge(Vertex v1, Vertex v2, double dist){ 
			this.v1 = v1;
			this.v2 = v2;
			this.dist = dist;
		}
		public Vertex getV1() {	return v1;	}
		public Vertex getV2() {	return v2;	}
		public double getDist() {  return dist;	}
	}
}
