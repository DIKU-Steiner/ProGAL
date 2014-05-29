package ProGAL.geom2d.tests;

import static org.junit.Assert.*;

import org.junit.Test;

import ProGAL.geom2d.Triangle;
import ProGAL.math.Constants;

public class TriangleTest {

	@Test
	public void testTriangle() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetAltitude() {
		fail("Not yet implemented");
	}

	@Test
	public void testCalculateArea() {
		assertEquals( 3.0, Triangle.calculateArea(2, 3, Math.sqrt(13)), Constants.EPSILON );
		assertEquals( Math.sqrt(7.5*2.5*2.5*2.5), Triangle.calculateArea(5.0,5.0,5.0), Constants.EPSILON );
	}

	@Test
	public void testCalculateHeight() {
		assertEquals( 2.0, Triangle.calculateHeight(Math.sqrt(13), Math.sqrt(8), 5.0), Constants.EPSILON );
	}

	@Test
	public void testGetCos() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetCircumCircle() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetSteinerPoint() {
		fail("Not yet implemented");
	}

	@Test
	public void testInCircumCircle() {
		fail("Not yet implemented");
	}

}
