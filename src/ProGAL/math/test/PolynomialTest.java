package ProGAL.math.test;

import static org.junit.Assert.*;

import java.util.Arrays;

import org.junit.Test;

import ProGAL.math.Constants;
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
		//Test that cubic generalizes to square
		double[] roots = Polynomial.calcRoots(0,0,0, 0);
//		assertEquals(1, roots.length);
		assertEquals(Double.NaN,roots[0],0.00001);
//		assertTrue(roots[3]>0);

		roots = Polynomial.calcRoots(0,3,0,0);
		assertEquals(1,roots.length);
		assertEquals(0.0,roots[0],0.00001);

		roots = Polynomial.calcRoots(0,-5,0,0);
		assertEquals(1,roots.length);
		assertEquals(0.0,roots[0],0.00001);
		
		roots = Polynomial.calcRoots(0,2,0,-3);
		assertEquals(2,roots.length);
		assertEquals(-Math.sqrt(3.0/2.0),roots[0],0.00001);//2x^2 = 3 -> x = +- sqrt(3/2)
		assertEquals(Math.sqrt(3.0/2.0),roots[1],0.00001);//2x^2 = 3 -> x = +- sqrt(3/2)
		
		
		//Now test some simple cubic roots
		roots = Polynomial.calcRoots(1, 10, 0, 0);
		assertEquals(2,roots.length);
		assertEquals(-10, roots[0],Constants.EPSILON);
		assertEquals(0, roots[1],Constants.EPSILON);

		roots = Polynomial.calcRoots(1, 0, 1, 0);
		assertEquals(1,roots.length);
		assertEquals(0, roots[0],Constants.EPSILON);
		

		roots = Polynomial.calcRoots(1, 0, -3, 0);
		assertEquals(3,roots.length);
		assertEquals(-Math.sqrt(3), roots[0],Constants.EPSILON);
		assertEquals(0, roots[1],Constants.EPSILON);
		assertEquals(Math.sqrt(3), roots[2],Constants.EPSILON);
		
		roots = Polynomial.calcRoots(1, 0, 0, 10);
		System.out.println(Arrays.toString(roots));
		assertEquals(1,roots.length);
		assertEquals(Math.cbrt(-10d), roots[0],Constants.EPSILON);
		
	}

}
