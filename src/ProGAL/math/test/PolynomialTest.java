package ProGAL.math.test;

import static org.junit.Assert.*;

import org.junit.Test;

import ProGAL.math.Polynomial;

public class PolynomialTest {

	@Test
	public void testCalcRootsDoubleDoubleDouble() {
		double[] roots = Polynomial.calcRoots(0,0, 0);
		assertEquals(1, roots.length);
		assertEquals(Double.NaN,roots[0],0.00001);

		roots = Polynomial.calcRoots(3,0,0);
		assertEquals(1,roots.length);
		assertEquals(0.0,roots[0],0.00001);

		roots = Polynomial.calcRoots(-5,0,0);
		assertEquals(1,roots.length);
		assertEquals(0.0,roots[0],0.00001);
		
		roots = Polynomial.calcRoots(2,0,-3);
		assertEquals(2,roots.length);
		assertEquals(-Math.sqrt(3.0/2.0),roots[0],0.00001);//2x^2 = 3 -> x = +- sqrt(3/2)
		assertEquals(Math.sqrt(3.0/2.0),roots[1],0.00001);//2x^2 = 3 -> x = +- sqrt(3/2)
	}

	@Test
	public void testCalcRootsDoubleDoubleDoubleDouble() {
		fail("Not yet implemented");
	}

}
