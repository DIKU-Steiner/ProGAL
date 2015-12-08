package ProGAL.geom3d.volumes.tests;

import static org.junit.Assert.*;

import org.junit.Test;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.volumes.Lens;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.math.Constants;

public class LensTest {

	@Test
	public void testLens() {
		Sphere s0 = new Sphere(new Point(0,1,0),2);
		Sphere s1 = new Sphere(new Point(2,1,0),1);
		Lens l = new Lens(s0, s1);
		double focalDist1 = l.getFocalDistance(0);
		double r1 = l.getRadius();
		s0.getCenter().set(0, -1);
		double focalDist2 = l.getFocalDistance(0);
		double r2 = l.getRadius();
		//Checks that the s0 and s1 pointers are not used in Lens
		assertEquals(0.0, focalDist1-focalDist2, Constants.EPSILON);
		assertEquals(0.0, r1-r2, Constants.EPSILON);
		
		try{
			s0 = new Sphere(new Point(0,1,0),Math.sqrt(2));
			s1 = new Sphere(new Point(3,1,0),Math.sqrt(2));
			l = new Lens(s0, s1);//Should fail with IllegalArgumentException
			assertFalse(true);
		}catch(IllegalArgumentException exc){
			assertTrue(true);
		}
		try{
			s0 = new Sphere(new Point(0,1,0),Math.sqrt(2));
			s1 = new Sphere(new Point(1,1,0),0.3);//Contained in s0
			l = new Lens(s0, s1);//Should fail with IllegalArgumentException
			assertFalse(true);
		}catch(IllegalArgumentException exc){
			assertTrue(true);
		}
		
		
	}

	@Test
	public void testGetRadius() {
		Sphere s0 = new Sphere(new Point(0,1,0),Math.sqrt(2));
		Sphere s1 = new Sphere(new Point(2,1,0),Math.sqrt(2));
		Lens l = new Lens(s0, s1);
		assertEquals(1, l.getRadius(), Constants.EPSILON);
		
		s1.getCenter().setX(3);
		s1.setRadius(Math.sqrt(5));
		l = new Lens(s0, s1);
		assertEquals(1, l.getRadius(), Constants.EPSILON);
	}

	@Test
	public void testGetFocalDistance() {
		Sphere s0 = new Sphere(new Point(0,1,0),Math.sqrt(2));
		Sphere s1 = new Sphere(new Point(2,1,0),Math.sqrt(2));
		Lens l = new Lens(s0, s1);
		assertEquals(1, l.getFocalDistance(0), Constants.EPSILON);
		assertEquals(1, l.getFocalDistance(1), Constants.EPSILON);
		
		s1.getCenter().setX(3);
		s1.setRadius(Math.sqrt(5));
		l = new Lens(s0, s1);
		assertEquals(1, l.getFocalDistance(0), Constants.EPSILON);
		assertEquals(2, l.getFocalDistance(1), Constants.EPSILON);
	}

	@Test
	public void testGetCenter() {
		Sphere s0 = new Sphere(new Point(0,1,0),Math.sqrt(2));
		Sphere s1 = new Sphere(new Point(2,1,0),Math.sqrt(2));
		Lens l = new Lens(s0, s1);
		assertEquals(new Point(1,1,0), l.getCenter());
		
		s1.getCenter().setX(3);
		s1.setRadius(Math.sqrt(5));
		l = new Lens(s0, s1);
		assertEquals(new Point(1,1,0), l.getCenter());
		
	}

	@Test
	public void testGetVolume() {
		Sphere a0 = new Sphere(new Point(0,1,0), 1);
		Sphere a1 = new Sphere(new Point(1,1,0), 1);
		Lens A = new Lens(a0,a1);
		assertTrue(A.getVolume()<2 && A.getVolume()>1);//Not the test I'm proudest of
	}

	@Test
	public void testDistance() {
		Sphere a0 = new Sphere(new Point(0,1,0), Math.sqrt(2));
		Sphere a1 = new Sphere(new Point(2,1,0), Math.sqrt(2));
		Lens A = new Lens(a0,a1), a = new Lens(a1, a0);

		Sphere b0,b1;
		Lens B,b;
		
		//Case 1: Distance should be the sphere distance from a0 to b1.
		b0 = new Sphere(new Point(2,1,0), Math.sqrt(2));
		b1 = new Sphere(new Point(4,1,0), Math.sqrt(2));
		B = new Lens(b0,b1); b = new Lens(b1,b0);
		assertEquals(4-Math.sqrt(2)*2, A.distance(B), Constants.EPSILON);
		assertEquals(4-Math.sqrt(2)*2, A.distance(b), Constants.EPSILON);
		assertEquals(4-Math.sqrt(2)*2, a.distance(B), Constants.EPSILON);
		assertEquals(4-Math.sqrt(2)*2, a.distance(b), Constants.EPSILON);
		
		//Case 2: Distance should be from the equator of B to the sphere of A (a0). 
		b0 = new Sphere(new Point(2,4,0), Math.sqrt(2));
		b1 = new Sphere(new Point(6,4,0), Math.sqrt(10));
		B = new Lens(b0,b1); b = new Lens(b1,b0);
		assertEquals(Math.sqrt(13)-Math.sqrt(2), A.distance(B), Constants.EPSILON);
		assertEquals(Math.sqrt(13)-Math.sqrt(2), A.distance(b), Constants.EPSILON);
		assertEquals(Math.sqrt(13)-Math.sqrt(2), a.distance(B), Constants.EPSILON);
		assertEquals(Math.sqrt(13)-Math.sqrt(2), a.distance(b), Constants.EPSILON);
		
		//Case 3: Distance should be from the equator of B to the equator of A. 
		b0 = new Sphere(new Point(0.5,-2,0), Math.sqrt(2));
		b1 = new Sphere(new Point(4.5,-2,0), Math.sqrt(10));
		B = new Lens(b0,b1); b = new Lens(b1,b0);
		assertEquals(Math.sqrt(1.25), A.distance(B), Constants.EPSILON);
		assertEquals(Math.sqrt(1.25), A.distance(b), Constants.EPSILON);
		assertEquals(Math.sqrt(1.25), a.distance(B), Constants.EPSILON);
		assertEquals(Math.sqrt(1.25), a.distance(b), Constants.EPSILON);

		//Case 4: Overlap but not between discs
		b0 = new Sphere(new Point(0.5,1,0), Math.sqrt(2));
		b1 = new Sphere(new Point(2.5,1,0), Math.sqrt(2));
		B = new Lens(b0,b1); b = new Lens(b1,b0);
		assertEquals(0, A.distance(B), Constants.EPSILON);
		assertEquals(0, A.distance(b), Constants.EPSILON);
		assertEquals(0, a.distance(B), Constants.EPSILON);
		assertEquals(0, a.distance(b), Constants.EPSILON);
		
		//Case 5: Overlap between disc and sphere
		b0 = new Sphere(new Point(0.01,2.8,0), Math.sqrt(2));
		b1 = new Sphere(new Point(4.01,2.8,0), Math.sqrt(10));
		B = new Lens(b0,b1); b = new Lens(b1,b0);
		assertEquals(0, A.distance(B), Constants.EPSILON);
		assertEquals(0, A.distance(b), Constants.EPSILON);
		assertEquals(0, a.distance(B), Constants.EPSILON);
		assertEquals(0, a.distance(b), Constants.EPSILON);
		
		//Case 6: Overlap between discs
		b0 = new Sphere(new Point(1.01,-1,-1), Math.sqrt(2));
		b1 = new Sphere(new Point(1.01,-1,3), Math.sqrt(10));
		B = new Lens(b0,b1); b = new Lens(b1,b0);
		assertEquals(0, A.distance(B), Constants.EPSILON);
		assertEquals(0, A.distance(b), Constants.EPSILON);
		assertEquals(0, a.distance(B), Constants.EPSILON);
		assertEquals(0, a.distance(b), Constants.EPSILON);
		
	}

}
