package ProGAL.geom2d.tests;

import static org.junit.Assert.*;

import org.junit.Test;

import ProGAL.geom2d.Vector;

public class PointTest {

	@Test
	public void testAddVector() {
		fail("Not yet implemented");
	}

	@Test
	public void testAddDoubleDouble() {
		fail("Not yet implemented");
	}

	@Test
	public void testAddThisPoint() {
		fail("Not yet implemented");
	}

	@Test
	public void testAddThisDoubleDouble() {
		fail("Not yet implemented");
	}

	@Test
	public void testSubtract() {
		fail("Not yet implemented");
	}

	@Test
	public void testSubtractThisVector() {
		fail("Not yet implemented");
	}

	@Test
	public void testSubtractThisPoint() {
		fail("Not yet implemented");
	}

	@Test
	public void testArea() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetSquaredDistance() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetSquaredDistancePoint() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetSquaredDistanceLineSegment() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetDistance() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetSignedAngle() {
		fail("Not yet implemented");
	}

	@Test
	public void testInCircle() {
		fail("Not yet implemented");
	}

	@Test
	public void testPolarAngle() {
		fail("Not yet implemented");
	}

	@Test
	public void testPolarAngleSin() {
		fail("Not yet implemented");
	}

	@Test
	public void testPolarAngleCos() {
		fail("Not yet implemented");
	}

	@Test
	public void testLeftTurn() {
		fail("Not yet implemented");
	}

	@Test
	public void testRightTurn() {
		fail("Not yet implemented");
	}

	@Test
	public void testMidPoint() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetBisector() {
		fail("Not yet implemented");
	}

	@Test
	public void testEqualsPoint() {
		fail("Not yet implemented");
	}

	@Test
	public void testRotationDouble() {
		Vector v = new Vector(2,1);
		v.rotateThis(Math.PI/2.0);
		assertEquals(-1, v.x(), 0.000001);
		assertEquals( 2, v.y(), 0.000001);
		
		v = new Vector(1,0);
		v.rotateThis(-Math.PI/4.0);
		assertEquals( 0.70710678118, v.x(), 0.00001);
		assertEquals(-0.70710678118, v.y(), 0.00001);
	}

	@Test
	public void testRotationClone() {
		fail("Not yet implemented");
	}

	@Test
	public void testRotationDoubleDouble() {
		fail("Not yet implemented");
	}

	@Test
	public void testRotationPointDouble() {
		fail("Not yet implemented");
	}

}
