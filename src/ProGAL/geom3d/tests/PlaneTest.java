package ProGAL.geom3d.tests;

import static org.junit.Assert.*;

import org.junit.Test;

import ProGAL.geom3d.Line;
import ProGAL.geom3d.Plane;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.Vector;
import ProGAL.math.Constants;

public class PlaneTest {

	@Test
	public void testPlaneVector() {
		Plane p = new Plane(new Vector(0,3,3));
		assertEquals(new Point(0,0,0), p.getPoint());
		assertEquals(new Vector(0,3,3), p.getNormal());
	}
	@Test
	public void testPlanePointVector() {
		Plane p = new Plane(new Point(1,1,1), new Vector(0,3,3));
		assertEquals(new Point(1,1,1), p.getPoint());
		assertEquals(new Vector(0,3,3), p.getNormal());
	}

	@Test
	public void testPlanePointPointPoint() {
		Plane p = new Plane(new Point(1,1,1), new Point(2,1,1), new Point(3,2,1));
		assertEquals(new Point(1,1,1), p.getPoint());
		assertTrue(	new Vector(0,0,1).equals(p.getNormal().normalize()) ||  
					new Vector(0,0,-1).equals(p.getNormal().normalize())  );
	}

	@Test
	public void testGetPoint() {
		Point point = new Point(1,1,1);
		Plane p = new Plane(point, new Vector(0,3,3));
		assertTrue(point==p.getPoint());
		assertEquals(new Point(1,1,1), p.getPoint());
		assertFalse(new Point(1,0,0)==new Point(1,0,0));
	}

	@Test
	public void testGetNormal() {
		Vector norm = new Vector(0,3,3);
		Plane p = new Plane(new Point(1,1,1),norm);
		assertTrue(norm==p.getNormal());
		assertEquals(new Vector(0,3,3),p.getNormal());
	}

	@Test
	public void testSetNormal() {
		Vector norm = new Vector(0,3,3);
		Plane p = new Plane(new Point(1,1,1),norm);
		assertTrue(norm==p.getNormal());
		assertEquals(new Vector(0,3,3),p.getNormal());
		Vector norm2 = new Vector(1,0,0);
		p.setNormal(norm2);
		assertTrue(norm2==p.getNormal());
		assertFalse(norm==p.getNormal());
		assertEquals(new Vector(1,0,0),p.getNormal());
	}

	@Test
	public void testProjectPoint() {
		Plane p = new Plane(new Point(1,1,1), new Vector(0,0,3));
		Point proj = p.projectPoint(new Point(0,0,20));
		assertEquals(new Point(0,0,1), proj);
		
		p = new Plane(new Point(1,1,1), new Vector(2,2,2));
		assertEquals(new Point(1,1,1), p.projectPoint(new Point(0,0,0)));
		assertEquals(new Point(1,1,1), p.projectPoint(new Point(20,20,20)));
		assertEquals(new Point(1,1,1), p.projectPoint(new Point(1,1,1)));
		
		//Now a non-trivial point-projection. Since these are not trivial for me i use a line
		Line l = new Line(new Point(2,0,0), new Vector(1,1,1));
		assertEquals(  p.getIntersection(l), p.projectPoint(l.getPoint(10))  );
		//Also check that projections from the back-side of the plane gets projected correct
		assertEquals(  p.getIntersection(l), p.projectPoint(l.getPoint(-10))  );
	}

	@Test
	public void testAbove() {
		Plane p = new Plane(new Point(1,1,1), new Vector(2,2,2));
		assertEquals(  1, p.above(new Point(1.1,1,1)) );
		assertEquals(  0, p.above(new Point(1,1,1)) );
		assertEquals( -1, p.above(new Point(0.99,1,1)) );
		assertEquals(  0, p.above(p.getPoint()) );
	}

	@Test
	public void testBelow() {
		Plane p = new Plane(new Point(1,1,1), new Vector(2,2,2));
		assertEquals( -1, p.below(new Point(1.1,1,1)) );
		assertEquals(  0, p.below(new Point(1,1,1)) );
		assertEquals(  1, p.below(new Point(0.99,1,1)) );
		assertEquals(  0, p.below(p.getPoint()) );
	}

	@Test
	public void testGetDistance() {
		Plane p = new Plane(new Point(1,1,1), new Vector(2,2,2));
		assertEquals(Math.sqrt(3)*2, p.getDistance(new Point(3,3,3)), Constants.EPSILON);
		assertEquals(Math.sqrt(3), p.getDistance(new Point(0,0,0)), Constants.EPSILON);
		assertEquals(0, p.getDistance(new Point(3,0,0)), Constants.EPSILON);
	}

	@Test
	public void testgetUnsignedDihedralAngle() {
		Plane p;
		Plane p1 = new Plane(new Point(1,1,1), new Vector(2,2,2));
		p = new Plane(new Point(4,0,0), new Vector(1,1,0));
		assertEquals(Math.PI/4, p1.getUnsignedDihedralAngle(p), Constants.EPSILON);
		p = new Plane(new Point(-2,1,1), new Vector(1,1,0));
		assertEquals(Math.PI/4, p1.getUnsignedDihedralAngle(p), Constants.EPSILON);

		//Ensure that it is unsigned
		p = new Plane(new Point(1,1,1), new Vector(-1,-1,0));
		assertEquals(Math.PI/4, p1.getUnsignedDihedralAngle(p), Constants.EPSILON);
		
	}

	@Test
	public void testGetIntersectionLine() {
		
	}

	@Test
	public void testGetIntersectionParameter() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetIntersectionLineSegment() {
		fail("Not yet implemented");
	}

}
