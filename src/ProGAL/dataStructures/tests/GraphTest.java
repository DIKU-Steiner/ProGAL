package ProGAL.dataStructures.tests;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.junit.Test;

import ProGAL.dataStructures.Graph;

public class GraphTest {

	@Test
	public void testGraphIntBoolean() {
		Graph g = new Graph(5, true);
		g.addEdge(0, 1);
		g.addEdge(1, 2);
		g.addEdge(1, 0);
		assertTrue(g.adjacent(0, 1));
		assertTrue(g.adjacent(1, 2));
		assertTrue(g.adjacent(1, 0));
		assertFalse(g.adjacent(2, 1));
		assertFalse(g.adjacent(4, 1));
		assertFalse(g.adjacent(1, 4));
		assertFalse(g.adjacent(4, 3));
		assertEquals(5, g.V());
		assertEquals(3, g.E());
		
		g = new Graph(5, false);
		g.addEdge(0, 1);
		g.addEdge(1, 2);
		g.addEdge(1, 0);
		assertTrue(g.adjacent(0, 1));
		assertTrue(g.adjacent(1, 2));
		assertTrue(g.adjacent(1, 0));
		assertTrue(g.adjacent(2, 1));
		assertFalse(g.adjacent(4, 1));
		assertFalse(g.adjacent(1, 4));
		assertFalse(g.adjacent(4, 3));
		assertEquals(5, g.V());
		assertEquals(3, g.E());
	}

	@Test
	public void testGraphGraph() {
		Graph g1 = new Graph(5, true);
		g1.addEdge(0, 1);
		g1.addEdge(1, 2);
		g1.addEdge(1, 0);
		
		Graph g2 = new Graph(g1);
		assertTrue(g2.adjacent(0, 1));
		assertTrue(g2.adjacent(1, 2));
		assertTrue(g2.adjacent(1, 0));
		assertFalse(g2.adjacent(2, 1));
		assertFalse(g2.adjacent(4, 1));
		assertFalse(g2.adjacent(1, 4));
		assertFalse(g2.adjacent(4, 3));
		assertEquals(5, g2.V());
		assertEquals(3, g2.E());
	}

	@Test
	public void testV() {
		Graph g = new Graph(5, false);
		assertEquals(5, g.V());
		g.addEdge(0, 1);
		g.addEdge(1, 2);
		g.addEdge(1, 0);
		assertEquals(5, g.V());
	}

	@Test
	public void testE() {
		Graph g = new Graph(5, false);
		assertEquals(0, g.E());
		g.addEdge(0, 1);
		g.addEdge(1, 2);
		g.addEdge(1, 0);
		assertEquals(3, g.E());
		g.removeEdge(1, 0);
		assertEquals(2, g.E());
		g.removeEdge(1, 0);
		assertEquals(1, g.E());
	}

	@Test
	public void testAddEdge() {
		//Identical to Graph(int,int) test
		Graph g = new Graph(5, true);
		g.addEdge(0, 1);
		g.addEdge(1, 2);
		g.addEdge(1, 0);
		assertTrue(g.adjacent(0, 1));
		assertTrue(g.adjacent(1, 2));
		assertTrue(g.adjacent(1, 0));
		assertFalse(g.adjacent(2, 1));
		assertFalse(g.adjacent(4, 1));
		assertFalse(g.adjacent(1, 4));
		assertFalse(g.adjacent(4, 3));
		assertEquals(5, g.V());
		assertEquals(3, g.E());
		
		g = new Graph(5, false);
		g.addEdge(0, 1);
		g.addEdge(1, 2);
		g.addEdge(1, 0);
		assertTrue(g.adjacent(0, 1));
		assertTrue(g.adjacent(1, 2));
		assertTrue(g.adjacent(1, 0));
		assertTrue(g.adjacent(2, 1));
		assertFalse(g.adjacent(4, 1));
		assertFalse(g.adjacent(1, 4));
		assertFalse(g.adjacent(4, 3));
		assertEquals(5, g.V());
		assertEquals(3, g.E());
	}

	@Test
	public void testRemoveEdge() {
		Graph g = new Graph(5, false);
		assertEquals(0, g.E());
		g.addEdge(0, 1);
		g.addEdge(1, 2);
		g.addEdge(1, 0);
		assertEquals(3, g.E());
		g.removeEdge(1, 0);
		assertEquals(2, g.E());
		g.removeEdge(1, 0);
		assertEquals(1, g.E());
		
		g = new Graph(5, true);
		assertEquals(0, g.E());
		g.addEdge(0, 1);
		g.addEdge(1, 2);
		g.addEdge(1, 0);
		assertEquals(3, g.E());
		g.removeEdge(1, 0);
		assertEquals(2, g.E());
		g.removeEdge(1, 0);
		assertEquals(2, g.E());
		g.removeEdge(0, 1);
		assertEquals(1, g.E());
	}

	@Test
	public void testAdjacent() {
		Graph g = new Graph(5, false);
		assertFalse(g.adjacent(0, 1));assertFalse(g.adjacent(0, 2));assertFalse(g.adjacent(0, 3));
		g.addEdge(0, 1);
		assertTrue(g.adjacent(0, 1));assertFalse(g.adjacent(0, 2));assertFalse(g.adjacent(0, 3));
		g.addEdge(1, 2);
		assertTrue(g.adjacent(0, 1));assertFalse(g.adjacent(0, 2));assertFalse(g.adjacent(0, 3));
		assertTrue(g.adjacent(1, 2));assertTrue(g.adjacent(2, 1));
		g.addEdge(1, 0);
		assertTrue(g.adjacent(0, 1));assertFalse(g.adjacent(0, 2));assertFalse(g.adjacent(0, 3));
		

		g = new Graph(5, true);
		assertFalse(g.adjacent(0, 1));assertFalse(g.adjacent(0, 2));assertFalse(g.adjacent(0, 3));
		g.addEdge(0, 1);
		assertTrue(g.adjacent(0, 1));assertFalse(g.adjacent(0, 2));assertFalse(g.adjacent(0, 3));
		g.addEdge(1, 2);
		assertTrue(g.adjacent(0, 1));assertFalse(g.adjacent(0, 2));assertFalse(g.adjacent(0, 3));
		assertTrue(g.adjacent(1, 2));assertFalse(g.adjacent(2, 1));
		g.addEdge(1, 0);
		assertTrue(g.adjacent(0, 1));assertFalse(g.adjacent(0, 2));assertFalse(g.adjacent(0, 3));
	}

	@Test
	public void testAdj() {
		Graph g = new Graph(5, false);
		g.addEdge(0, 1);
		g.addEdge(1, 2);
		g.addEdge(1, 0);

		List<Integer> adj = new ArrayList<Integer>();
		for(Integer v: g.adj(0)) adj.add(v);
		Collections.sort(adj);
		assertEquals(Arrays.asList(1,1), adj);
		
		adj.clear();
		for(Integer v: g.adj(1)) adj.add(v);
		Collections.sort(adj);
		assertEquals(Arrays.asList(0,0,2), adj);

		adj.clear();
		for(Integer v: g.adj(2)) adj.add(v);
		Collections.sort(adj);
		assertEquals(Arrays.asList(1), adj);
		
		adj.clear();
		for(Integer v: g.adj(3)) adj.add(v);
		assertTrue(adj.isEmpty());
	
		
		g = new Graph(5, true);
		g.addEdge(0, 1);
		g.addEdge(1, 2);
		g.addEdge(1, 0);

		adj = new ArrayList<Integer>();
		for(Integer v: g.adj(0)) adj.add(v);
		Collections.sort(adj);
		assertEquals(Arrays.asList(1), adj);
		
		adj.clear();
		for(Integer v: g.adj(1)) adj.add(v);
		Collections.sort(adj);
		assertEquals(Arrays.asList(0,2), adj);

		adj.clear();
		for(Integer v: g.adj(2)) adj.add(v);
		assertTrue(adj.isEmpty());
		
		adj.clear();
		for(Integer v: g.adj(3)) adj.add(v);
		assertTrue(adj.isEmpty());
	}

	@Test
	public void testBreadthFirstIntInt() {
		Graph g = new Graph(5, false);
		g.addEdge(0, 1);
		g.addEdge(1, 2);
		g.addEdge(1, 0);
		
		ArrayList<Integer> gComp = g.breadthFirst(0, 1000);
		Collections.sort(gComp);
		assertEquals(Arrays.asList(1,2), gComp);
		
		gComp = g.breadthFirst(0, 1);
		Collections.sort(gComp);
		assertEquals(Arrays.asList(1), gComp);

		gComp = g.breadthFirst(2, 1000);
		Collections.sort(gComp);
		assertEquals(Arrays.asList(0,1), gComp);
		
		gComp = g.breadthFirst(3, 1000);
		assertTrue(gComp.isEmpty());
		
		
		g = new Graph(5, true);
		g.addEdge(0, 1);
		g.addEdge(1, 2);
		g.addEdge(1, 0);
		
		gComp = g.breadthFirst(0, 1000);
		Collections.sort(gComp);
		assertEquals(Arrays.asList(1,2), gComp);
		
		gComp = g.breadthFirst(0, 1);
		Collections.sort(gComp);
		assertEquals(Arrays.asList(1), gComp);

		gComp = g.breadthFirst(2, 1000);
		assertTrue(gComp.isEmpty());
		
		gComp = g.breadthFirst(3, 1000);
		assertTrue(gComp.isEmpty());
	}

}
