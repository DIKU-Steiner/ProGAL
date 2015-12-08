package ProGAL.geom3d.complex.delaunayComplex.tests;

import static org.junit.Assert.*;

import org.junit.Test;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.PointList;
import ProGAL.geom3d.complex.CEdge;
import ProGAL.geom3d.complex.delaunayComplex.DelaunayComplex;

public class DelaunayComplexTest {

	@Test
	public void testDelaunayComplex() {
		PointList pl = new PointList();
		pl.add(new Point(1,0,0));
		pl.add(new Point(2,0,0));
		pl.add(new Point(1,2,0));
		pl.add(new Point(1,0,3));
		pl.add(new Point(3,2,3));
		
		DelaunayComplex dc = new DelaunayComplex(pl);
		assertTrue(dc.checkTetrahedra());
		
		//The triangulation should contain two tetrahedra (and corresponding 0, 1 and 2-simplices)
		assertEquals(2,dc.getTetrahedra().size());
		assertEquals(7,dc.getTriangles().size());
		assertEquals(9,dc.getEdges().size());
		assertEquals(5,dc.getVertices().size());
		
		//All of the vertices should be (almost) on top of the points and come in the same order 
		assertTrue(pl.get(0).distance(dc.getVertices().get(0))<0.0001);
		assertTrue(pl.get(1).distance(dc.getVertices().get(1))<0.0001);
		assertTrue(pl.get(2).distance(dc.getVertices().get(2))<0.0001);
		assertTrue(pl.get(3).distance(dc.getVertices().get(3))<0.0001);
		assertTrue(pl.get(4).distance(dc.getVertices().get(4))<0.0001);

		assertEquals(3, dc.getVertices().get(0).getAdjacentEdges().size());
		assertEquals(4, dc.getVertices().get(1).getAdjacentEdges().size());
		assertEquals(4, dc.getVertices().get(2).getAdjacentEdges().size());
		assertEquals(4, dc.getVertices().get(3).getAdjacentEdges().size());
		assertEquals(3, dc.getVertices().get(4).getAdjacentEdges().size());

		//Test triangles adjacent to edges
		for(CEdge e: dc.getVertices().get(0).getAdjacentEdges())
			assertEquals(2, e.getAdjacentTriangles().size());
		for(CEdge e: dc.getVertices().get(1).getAdjacentEdges())
			assertTrue(e.getAdjacentTriangles().size()==2||e.getAdjacentTriangles().size()==3);
		for(CEdge e: dc.getVertices().get(2).getAdjacentEdges())
			assertTrue(e.getAdjacentTriangles().size()==2||e.getAdjacentTriangles().size()==3);
		for(CEdge e: dc.getVertices().get(4).getAdjacentEdges())
			assertEquals(2, e.getAdjacentTriangles().size());
		
	}
	
}
