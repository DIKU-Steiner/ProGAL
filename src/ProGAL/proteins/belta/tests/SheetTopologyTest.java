package ProGAL.proteins.belta.tests;

import static org.junit.Assert.*;

import java.util.Arrays;

import org.junit.Test;

import ProGAL.proteins.belta.BetaTopology;
import ProGAL.proteins.belta.SecondaryStructure;
import ProGAL.proteins.belta.SheetTopology;

public class SheetTopologyTest {

	@Test
	public void testSheetTopology() {
		SecondaryStructure ss = new SecondaryStructure(" EE EE EE EE EE ");
		BetaTopology bTop = new BetaTopology(ss);
		bTop.setPaired(0, 1);
		bTop.setPaired(2, 3);
		bTop.setPaired(3, 4);
		SheetTopology sTop = new SheetTopology(bTop, 0);
		assertTrue(sTop.secondaryStructure==ss);
		assertEquals(2, sTop.strands.size());
		assertEquals(0, sTop.strandPairs.get(0).strand1);
		assertEquals(1, sTop.strandPairs.get(0).strand2);
		assertFalse(sTop.strandPairs.get(0).parallel);
		
		sTop = new SheetTopology(bTop, 1);
		assertEquals(2, sTop.strands.size());
		assertEquals(0, sTop.strandPairs.get(0).strand1);
		assertEquals(1, sTop.strandPairs.get(0).strand2);
		assertFalse(sTop.strandPairs.get(0).parallel);
		
		sTop = new SheetTopology(bTop, 3);
		assertEquals(3, sTop.strands.size());
		assertEquals(2, sTop.strandPairs.get(0).strand1);
		assertEquals(3, sTop.strandPairs.get(0).strand2);
		assertEquals(3, sTop.strandPairs.get(1).strand1);
		assertEquals(4, sTop.strandPairs.get(1).strand2);
		assertFalse(sTop.strandPairs.get(0).parallel);
		assertFalse(sTop.strandPairs.get(1).parallel);
		
		bTop.setNotPaired(3, 4);
		bTop.setPaired(4, 3);
		sTop = new SheetTopology(bTop, 3);
		assertEquals(3, sTop.strands.size());
		assertEquals(2, sTop.strandPairs.get(0).strand1);
		assertEquals(3, sTop.strandPairs.get(0).strand2);
		assertEquals(3, sTop.strandPairs.get(1).strand1);
		assertEquals(4, sTop.strandPairs.get(1).strand2);
		assertFalse(sTop.strandPairs.get(0).parallel);
		assertTrue(sTop.strandPairs.get(1).parallel);
		
	}

	@Test
	public void testContainsStrand() {
		SecondaryStructure ss = new SecondaryStructure(" EE EE EE EE EE ");
		BetaTopology bTop = new BetaTopology(ss);
		bTop.setPaired(0, 1);
		bTop.setPaired(2, 3);
		bTop.setPaired(4, 3);
		SheetTopology sTop = new SheetTopology(bTop,3);
		assertFalse(sTop.containsStrand(0));
		assertFalse(sTop.containsStrand(1));
		assertTrue(sTop.containsStrand(2));
		assertTrue(sTop.containsStrand(3));
		assertTrue(sTop.containsStrand(4));
	}

	@Test
	public void testGetStrandOrder() {
		SecondaryStructure ss = new SecondaryStructure(" EE EE EE EE EE ");
		BetaTopology bTop = new BetaTopology(ss);
		bTop.setPaired(0, 2);
		bTop.setPaired(1, 3);
		bTop.setPaired(4, 3);
		SheetTopology sTop = new SheetTopology(bTop,3);
		assertTrue(Arrays.equals(new int[]{1,3,4}, sTop.getStrandOrder()));
		sTop = new SheetTopology(bTop,0);
		assertTrue(Arrays.equals(new int[]{0,2}, sTop.getStrandOrder()));
		
		bTop = new BetaTopology(ss);
		bTop.setPaired(0, 2);
		bTop.setPaired(1, 3);
		bTop.setPaired(4, 1);
		sTop = new SheetTopology(bTop,3);
		assertTrue(Arrays.equals(new int[]{3,1,4}, sTop.getStrandOrder()));
		sTop = new SheetTopology(bTop,0);
		assertTrue(Arrays.equals(new int[]{0,2}, sTop.getStrandOrder()));
	}

	@Test
	public void testGetNormalizedStrandOrder() {
		SecondaryStructure ss = new SecondaryStructure(" EE EE EE EE EE ");
		BetaTopology bTop = new BetaTopology(ss);
		bTop.setPaired(0, 2);
		bTop.setPaired(1, 3);
		bTop.setPaired(4, 1);
		SheetTopology sTop = new SheetTopology(bTop,3);
		assertTrue(Arrays.equals(new int[]{1,0,2}, sTop.getNormalizedStrandOrder()));
		sTop = new SheetTopology(bTop,0);
		assertTrue(Arrays.equals(new int[]{0,1}, sTop.getNormalizedStrandOrder()));	
	}

	@Test
	public void testGetParTopology() {
		SecondaryStructure ss = new SecondaryStructure(" EE EE EE EE EE ");
		BetaTopology bTop = new BetaTopology(ss);
		bTop.setPaired(0, 2);
		bTop.setPaired(1, 3);
		bTop.setPaired(4, 1);
		SheetTopology sTop = new SheetTopology(bTop,3);
		assertTrue(Arrays.equals(new int[]{0,1,1}, sTop.getStrandOrientation()));
		sTop = new SheetTopology(bTop,0);
		assertTrue(Arrays.equals(new int[]{1,0}, sTop.getStrandOrientation()));
	}

}
