package ProGAL.dataStructures.tests;

import static org.junit.Assert.*;
import org.junit.Test;
import java.util.ArrayList;

import ProGAL.dataStructures.Kruskal;

public class KruskalTest {
	private static class Vertex implements Kruskal.Vertex{}
	
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
	
	
	@Test
	public void testKruskal() {
		ArrayList<Vertex> vertices = new ArrayList<Vertex>();
		ArrayList<Edge> edges = new ArrayList<Edge>(); 
		
		for(int i=0;i<5;i++)
			vertices.add(new Vertex());
		edges.add(new Edge(vertices.get(0), vertices.get(1), 3));
		edges.add(new Edge(vertices.get(0), vertices.get(4), 1));
		edges.add(new Edge(vertices.get(1), vertices.get(2), 5));
		edges.add(new Edge(vertices.get(1), vertices.get(4), 4));
		edges.add(new Edge(vertices.get(2), vertices.get(3), 2));
		edges.add(new Edge(vertices.get(2), vertices.get(4), 6));
		edges.add(new Edge(vertices.get(3), vertices.get(4), 7));
		
		Kruskal mst = new Kruskal(vertices,edges);
		assertEquals(11, mst.getLength(), 0.0000001);
		assertEquals(4, mst.getEdges().size());
		assertTrue( mst.getEdges().contains( edges.get(0) ) );
		assertTrue( mst.getEdges().contains( edges.get(1) ) );
		assertTrue( mst.getEdges().contains( edges.get(2) ) );
		assertTrue( mst.getEdges().contains( edges.get(4) ) );
		
	}

}
