package ProGAL.geomNd.tests;

import static org.junit.Assert.*;

import org.junit.Test;

import ProGAL.geomNd.*;
import ProGAL.math.Constants;

public class PointTest {

	@Test
	public void testPointDoubleArray() {
		double[] arr = new double[]{0.1,0.2,0.3,0.5};
		Point p = new Point(arr);
		assertEquals(0.1, p.get(0), Constants.EPSILON);
		assertEquals(0.2, p.get(1), Constants.EPSILON);
		assertEquals(0.3, p.get(2), Constants.EPSILON);
		assertEquals(0.5, p.get(3), Constants.EPSILON);
		assertEquals(4, p.getDimensions());
		
		arr[1] = 2.2;
		assertEquals(2.2, p.get(1), Constants.EPSILON);
	}

	@Test
	public void testPointInt() {
		Point p = new Point(5);
		assertEquals(0.0, p.get(0), Constants.EPSILON);
		assertEquals(0.0, p.get(1), Constants.EPSILON);
		assertEquals(0.0, p.get(2), Constants.EPSILON);
		assertEquals(0.0, p.get(3), Constants.EPSILON);
		assertEquals(0.0, p.get(4), Constants.EPSILON);
		assertEquals(4, p.getDimensions());
	}

	@Test
	public void testGet() {
		//Basically tested by the constructors. If specification changes this should be updated
	}

	@Test
	public void testGetCoord() {
		//Basically tested by the constructors. If specification changes this should be updated
	}

	@Test
	public void testSet() {
		double[] arr = new double[]{0.1,0.2,0.3,0.5};
		Point p = new Point(arr);
		assertEquals(0.1, p.get(0), Constants.EPSILON);
		assertEquals(0.2, p.get(1), Constants.EPSILON);
		assertEquals(0.3, p.get(2), Constants.EPSILON);
		assertEquals(0.5, p.get(3), Constants.EPSILON);
		assertEquals(4, p.getDimensions());
		
		p.set(1, 2.2);
		assertEquals(2.2, p.get(1), Constants.EPSILON);
		assertEquals(2.2, arr[1], Constants.EPSILON);
		
	}

	@Test
	public void testGetDistanceSquared() {
		Point p = new Point(new double[]{3,2,1,4});
		Point q = new Point(new double[]{-3,1,1,4});
		assertEquals(37, p.distanceSquared(q), Constants.EPSILON);
		assertEquals(0, p.distanceSquared(p), Constants.EPSILON);
		

		q = new Point(new double[]{10,20,30,40});
		assertEquals(7*7+18*18+29*29+36*36, p.distanceSquared(q), Constants.EPSILON);
	}

	@Test
	public void testGetDistance() {
		Point p = new Point(new double[]{3,2,1,4});
		Point q = new Point(new double[]{-3,1,1,4});
		assertEquals(Math.sqrt(37), p.distance(q), Constants.EPSILON);
		assertEquals(0, p.distance(p), Constants.EPSILON);
		

		q = new Point(new double[]{10,20,30,40});
		assertEquals(Math.sqrt(7*7+18*18+29*29+36*36), p.distanceSquared(q), Constants.EPSILON);
	}

	@Test
	public void testGetMidpoint() {
		Point p = new Point(new double[]{3,2,1,4});
		Point q = new Point(new double[]{0,2,1,4});
		Point midpoint = Point.getMidpoint(p,q);
		assertTrue(midpoint.distanceSquared(p)==midpoint.distanceSquared(q));
		assertTrue(	Point.collinear(p, q, midpoint) );
		
		q = p.clone();
		assertTrue(Point.getMidpoint(p, q).equals(p));
	}


	@Test
	public void testGetAngle() {
		//First a simple right rectangle in the xy-plane
		Point p = new Point(new double[]{1,0,0,3});
		Point q = new Point(new double[]{0,0,0,3});
		Point r = new Point(new double[]{0,1,0,3});
		assertEquals(Math.PI/2,	Point.getAngle(p, q, r), Constants.EPSILON);
		assertEquals(Math.PI/4,	Point.getAngle(q, r, p), Constants.EPSILON);
		assertEquals(Math.PI/4,	Point.getAngle(r, p, q), Constants.EPSILON);
		
		//Collinear points
		p = new Point(new double[]{1,1,1,3});
		q = new Point(new double[]{2,2,2,3});
		r = new Point(new double[]{4,4,4,3});
		assertEquals(Math.PI, 	Point.getAngle(p, q, r), Constants.EPSILON);
		assertEquals(0,			Point.getAngle(q, r, p), Constants.EPSILON);
		assertEquals(0,			Point.getAngle(r, p, q), Constants.EPSILON);
		
		//TODO: These are just 3D-examples translated into 4D. Make some better multi-dimensional tests
	}
	
	@Test
	public void testToVector() {
		Point p = new Point(new double[]{10,20,3.1415,1});
		Vector v = p.toVector();
		assertTrue(v.equals(new Vector(new double[]{10,20,3.1415,1})));
		p.set(0, 0);
		assertTrue(v.equals(new Vector(new double[]{10,20,3.1415,1})));//Ensure that the underlying array has been cloned
	}

}
