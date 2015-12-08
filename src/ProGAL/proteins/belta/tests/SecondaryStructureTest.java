package ProGAL.proteins.belta.tests;

import static org.junit.Assert.*;

import org.junit.Test;

import ProGAL.proteins.belta.SSType;
import ProGAL.proteins.belta.SecondaryStructure;
import ProGAL.proteins.belta.SecondaryStructure.SSSegment;

public class SecondaryStructureTest {

	@Test
	public void testSecondaryStructure() {
		//Constructing from 3-class
		SecondaryStructure ss = new SecondaryStructure("  HHHH EEEEE EEEE   HHHHHH HHHHHHH");
		assertEquals(10, ss.segments.length);
		assertEquals(2, ss.segments[0].length);
		assertEquals(0, ss.segments[0].start);
		assertEquals(2, ss.segments[0].end);
		assertEquals(SSType.COIL, ss.segments[0].type);
		assertEquals(6,  ss.segments[7].length);
		assertEquals(20, ss.segments[7].start);
		assertEquals(26, ss.segments[7].end);
		assertEquals(SSType.HELIX, ss.segments[7].type);
		for(int i=0;i<ss.segments.length;i++) 
			assertEquals(i, ss.segments[i].segmentIndex);

		//Constructing from DSSP
		ss = new SecondaryStructure("HHHHHHHHT  TTSSSEEEGGGHHHHHHHTT");
		assertEquals(5, ss.segments.length);
		assertEquals(8, ss.segments[0].length);
		assertEquals(0, ss.segments[0].start);
		assertEquals(8, ss.segments[0].end);
		assertEquals(SSType.HELIX, ss.segments[0].type);
		assertEquals(8,  ss.segments[1].length);
		assertEquals(8,  ss.segments[1].start);
		assertEquals(16, ss.segments[1].end);
		assertEquals(SSType.COIL, ss.segments[1].type);
		for(int i=0;i<ss.segments.length;i++) 
			assertEquals(i, ss.segments[i].segmentIndex);
	}

	@Test
	public void testGetStrands() {
		SecondaryStructure ss = new SecondaryStructure("  HHHH EEEEE EEEE   HHHHHH HHHHHHH");
		SSSegment[] strands = ss.getStrands();
		assertEquals(2, strands.length);
		assertEquals(3, strands[0].segmentIndex);
		assertEquals(5, strands[1].segmentIndex);
	}

	@Test
	public void testGetHelices() {
		SecondaryStructure ss = new SecondaryStructure("  HHHH EEEEE EEEE   HHHHHH HHHHHHH");
		SSSegment[] helices = ss.getHelices();
		assertEquals(3, helices.length);
		assertEquals(1, helices[0].segmentIndex);
		assertEquals(7, helices[1].segmentIndex);
		assertEquals(9, helices[2].segmentIndex);
	}

	@Test
	public void testGetCoils() {
		SecondaryStructure ss = new SecondaryStructure("  HHHH EEEEE EEEE   HHHHHH HHHHHHH");
		SSSegment[] coils = ss.getCoils();
		assertEquals(5, coils.length);
		assertEquals(0, coils[0].segmentIndex);
		assertEquals(2, coils[1].segmentIndex);
		assertEquals(4, coils[2].segmentIndex);
		assertEquals(6, coils[3].segmentIndex);
		assertEquals(8, coils[4].segmentIndex);
	}

	@Test
	public void testMatches(){
		SecondaryStructure ss1 = new SecondaryStructure("  HHHH EEEEE EEEE   HHHHHH HHHHHHH");
		SecondaryStructure ss2 = new SecondaryStructure("       EEEEE EEEE                 ");
		assertTrue(ss1.matches(ss2));
		assertTrue(ss2.matches(ss1));
		
		ss1 = new SecondaryStructure("  HHHH EEEEE EEEE   HHHHHH HHHHHHH");
		ss2 = new SecondaryStructure("    EEEEE      EEEE               ");
		assertTrue(ss1.matches(ss2));
		assertTrue(ss2.matches(ss1));
		
		ss1 = new SecondaryStructure("  HHHH EEEEE EEEE   HHHHHH HHHHHHH");
		ss2 = new SecondaryStructure("    EEEEE EEEE                    ");
		assertTrue(ss1.matches(ss2));
		assertTrue(ss2.matches(ss1));
		
		ss1 = new SecondaryStructure("  HHHH EEEEE EEEE   HHHHHH HHHHHHH");
		ss2 = new SecondaryStructure("     EEEEE EEEE     EEE           ");
		assertFalse(ss1.matches(ss2));
		assertFalse(ss2.matches(ss1));
		
		ss1 = new SecondaryStructure("  HHHH EEEEE EEEE   HHHHHH HHHHHHH");
		ss2 = new SecondaryStructure("     EEEEE HHHHH     EEE          ");
		assertFalse(ss1.matches(ss2));
		assertFalse(ss2.matches(ss1));
	}

	@Test
	public void testRespects(){
		SecondaryStructure ss1 = new SecondaryStructure("  HHHH EEEEE EEEE   HHHHHH HHHHHHH");
		SecondaryStructure ss2 = new SecondaryStructure("       EEEEE EEEE                 ");
		assertTrue(ss1.respects(ss2));
		assertTrue(ss2.respects(ss1));
		
		ss1 = new SecondaryStructure("  HHHH EEEEE EEEE   HHHHHH HHHHHHH");
		ss2 = new SecondaryStructure("    EEEEE      EEEE               ");
		assertTrue(ss1.respects(ss2));
		assertTrue(ss2.respects(ss1));
		
		ss1 = new SecondaryStructure("  HHHH EEEEE EEEE   HHHHHH HHHHHHH");
		ss2 = new SecondaryStructure("    EEEEE EEEE                    ");
		assertTrue(ss1.respects(ss2));
		assertTrue(ss2.respects(ss1));
		
		ss1 = new SecondaryStructure("  HHHH EEEEE EEEE   HHHHHH HHHHHHH");
		ss2 = new SecondaryStructure("     EEEEE EEEE     EEE           ");
		assertTrue(ss1.respects(ss2));
		assertFalse(ss2.respects(ss1));
		
		ss1 = new SecondaryStructure("  HHHH EEEEE EEEE   HHHHHH HHHHHHH");
		ss2 = new SecondaryStructure("     EEEEE HHHHH     EEE          ");
		assertFalse(ss1.respects(ss2));
		assertFalse(ss2.respects(ss1));
	}

}
