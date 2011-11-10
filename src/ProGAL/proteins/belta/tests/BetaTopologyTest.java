package ProGAL.proteins.belta.tests;

import static org.junit.Assert.*;

import java.util.List;

import org.junit.Test;

import ProGAL.proteins.belta.BetaTopology;
import ProGAL.proteins.belta.SecondaryStructure;
import ProGAL.proteins.belta.SheetTopology;

public class BetaTopologyTest {

	@Test
	public void testBetaTopology(){
		SecondaryStructure ss = new SecondaryStructure("  EEEE EEEE  EEE  EEE  ");
		BetaTopology bTop = new BetaTopology(ss);
		assertEquals(4, bTop.N);
		assertEquals(ss, bTop.secondaryStructure);
	}

	@Test
	public void testIsValid() {
		SecondaryStructure ss = new SecondaryStructure("  EEEE EEEE  EEE  EEE  ");
		BetaTopology bTop = new BetaTopology(ss);
		
		//There are four rules that makes a Beta topology valid: 
		// 1: No pair of strands are paired both parallel and antiparallel
		// 2: No strand is paired with itself
		// 3: Each strand has at least 1 partner
		// 4: Each strand has at most 2 partners
		//This is a holdout test where a matrix is set up to violate exactly one of these
		//rules in turn.
		
		//Valid matrices
		bTop.setPaired(0,1);
		bTop.setPaired(1,2);
		bTop.setPaired(2,3);
		assertTrue(bTop.isValid());
		bTop.setPaired(3,0);
		assertTrue(bTop.isValid());
		bTop.setNotPaired(0,1);
		bTop.setNotPaired(2,3);
		assertTrue(bTop.isValid());
		
		//Violating rule 1
		bTop = new BetaTopology(ss);
		bTop.setPaired(0,3);
		bTop.setPaired(3,0);
		bTop.setPaired(1,2);
		assertFalse(bTop.isValid());
		bTop.setPaired(2,1);
		assertFalse(bTop.isValid());
		
		//Violating rule 2
		bTop = new BetaTopology(ss);
		bTop.setPaired(2,0);
		bTop.setPaired(2,3);
		bTop.setPaired(1,1);
		assertFalse(bTop.isValid());
		
		//Violating rule 3
		bTop = new BetaTopology(ss);
		assertFalse(bTop.isValid());
		bTop.setPaired(0,1);
		assertFalse(bTop.isValid());
		bTop.setPaired(1,2);
		assertFalse(bTop.isValid());

		//Violating rule 4
		bTop = new BetaTopology(ss);
		bTop.setPaired(0,1);
		bTop.setPaired(0,2);
		bTop.setPaired(3,0);
		assertFalse(bTop.isValid());
		bTop.setNotPaired(3,0);
		bTop.setPaired(0,3);
		assertFalse(bTop.isValid());
		
	}

	@Test
	public void testGetSheets() {
		SecondaryStructure ss = new SecondaryStructure("  EEEE EEEE  EEE  EEE  ");
		BetaTopology bTop = new BetaTopology(ss);
		//An exception is expected if topology is not valid
		try{ bTop.getSheets(); assertTrue(false); }catch(RuntimeException exc){}
		
		//Just check the number of sheets. The actual sheet construction is tested in 
		//the SheetTopology-constructor.
		bTop.setPaired(0, 1);
		bTop.setPaired(2, 3);
		List<SheetTopology> sheets = bTop.getSheets();
		assertEquals(2, sheets.size());
		bTop.setPaired(1, 2);
		sheets = bTop.getSheets();
		assertEquals(1, sheets.size());
		
	}

	@Test
	public void testMatches() {
		SecondaryStructure ss1 = new SecondaryStructure("  EEEE EEEE  EEE  EEE  ");
		BetaTopology bTop1 = new BetaTopology(ss1);
		SecondaryStructure ss2 = new SecondaryStructure(" EEEE EEEE  EEE  EEEEE ");
		BetaTopology bTop2 = new BetaTopology(ss2);

		bTop1.setPaired(0, 1);
		bTop1.setPaired(2, 3);
		bTop2.setPaired(0, 1);
		bTop2.setPaired(2, 3);
		assertTrue(bTop1.matches(bTop2)); 
		assertTrue(bTop2.matches(bTop1));
		bTop1.setPaired(1, 2);
		assertFalse(bTop1.matches(bTop2)); 
		assertFalse(bTop2.matches(bTop1));
		bTop2.setPaired(2, 1);
		assertFalse(bTop1.matches(bTop2)); 
		assertFalse(bTop2.matches(bTop1));
		
		ss1 = new SecondaryStructure("  EEEE EEEE  EEE  EEE    ");
		bTop1 = new BetaTopology(ss1);
		ss2 = new SecondaryStructure("  EEE EEEE  EEE      EEE ");
		bTop2 = new BetaTopology(ss2);
		bTop1.setPaired(0, 1);
		bTop1.setPaired(2, 3);
		bTop2.setPaired(0, 1);
		bTop2.setPaired(2, 3);
		assertFalse(bTop1.matches(bTop2)); 
		assertFalse(bTop2.matches(bTop1));
	}

	@Test
	public void testRespects() {
		SecondaryStructure ss1 = new SecondaryStructure("  EEEE EEEE  EEE  EEE  ");
		BetaTopology bTop1 = new BetaTopology(ss1);
		SecondaryStructure ss2 = new SecondaryStructure(" EEEE EEEE  EEE  EEEEE ");
		BetaTopology bTop2 = new BetaTopology(ss2);

		bTop1.setPaired(0, 1);
		bTop1.setPaired(2, 3);
		bTop2.setPaired(0, 1);
		bTop2.setPaired(2, 3);
		assertTrue(bTop1.respects(bTop2));
		assertTrue(bTop2.respects(bTop1));
		bTop1.setPaired(1, 2);
		assertFalse(bTop1.respects(bTop2));
		assertTrue(bTop2.respects(bTop1));
		bTop2.setPaired(2, 1);
		assertFalse(bTop1.respects(bTop2));
		assertFalse(bTop2.respects(bTop1));
		
		ss1 = new SecondaryStructure("  EEEE EEEE  EEE  EEE    ");
		bTop1 = new BetaTopology(ss1);
		ss2 = new SecondaryStructure("  EEE EEEE  EEE      EEE ");
		bTop2 = new BetaTopology(ss2);
		bTop1.setPaired(0, 1);
		bTop1.setPaired(2, 3);
		bTop2.setPaired(0, 1);
		bTop2.setPaired(2, 3);
		assertFalse(bTop1.respects(bTop2));
		assertFalse(bTop2.respects(bTop1));
	
	}

}
