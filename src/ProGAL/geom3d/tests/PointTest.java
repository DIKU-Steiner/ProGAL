package ProGAL.geom3d.tests;

import static org.junit.Assert.*;

import org.junit.Test;

import ProGAL.geom3d.Plane;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.Vector;
import ProGAL.math.Constants;

public class PointTest {

	@Test
	public void testPoint3dDoubleDoubleDouble() {
		Point p = new Point(0.1, 0.2, 0.3);
		assertEquals(0.1, p.x(), Constants.EPSILON);
		assertEquals(0.1, p.get(0), Constants.EPSILON);
		assertEquals(0.2, p.y(), Constants.EPSILON);
		assertEquals(0.2, p.get(1), Constants.EPSILON);
		assertEquals(0.3, p.z(), Constants.EPSILON);
		assertEquals(0.3, p.get(2), Constants.EPSILON);
	}

	@Test
	public void testPoint3dPoint3d() {
		Point p = new Point(new Point(0.1, 0.2, 0.3));
		assertEquals(0.1, p.x(), Constants.EPSILON);
		assertEquals(0.1, p.get(0), Constants.EPSILON);
		assertEquals(0.2, p.y(), Constants.EPSILON);
		assertEquals(0.2, p.get(1), Constants.EPSILON);
		assertEquals(0.3, p.z(), Constants.EPSILON);
		assertEquals(0.3, p.get(2), Constants.EPSILON);
	}

	@Test
	public void testPoint3dVector3d() {
		Point p = new Point(new Vector(0.1, 0.2, 0.3));
		assertEquals(0.1, p.x(), Constants.EPSILON);
		assertEquals(0.1, p.get(0), Constants.EPSILON);
		assertEquals(0.2, p.y(), Constants.EPSILON);
		assertEquals(0.2, p.get(1), Constants.EPSILON);
		assertEquals(0.3, p.z(), Constants.EPSILON);
		assertEquals(0.3, p.get(2), Constants.EPSILON);
	}

	@Test
	public void testGetCoord() {
		Point p = new Point(new Vector(0.1, 0.2, 0.3));
		assertEquals(0.1, p.x(), Constants.EPSILON);
		assertEquals(0.1, p.get(0), Constants.EPSILON);
		assertEquals(0.2, p.y(), Constants.EPSILON);
		assertEquals(0.2, p.get(1), Constants.EPSILON);
		assertEquals(0.3, p.z(), Constants.EPSILON);
		assertEquals(0.3, p.get(2), Constants.EPSILON);
	}

	@Test
	public void testGet() {
		Point p = new Point(new Vector(0.1, 0.2, 0.3));
		assertEquals(0.1, p.x(), Constants.EPSILON);
		assertEquals(0.1, p.get(0), Constants.EPSILON);
		assertEquals(0.2, p.y(), Constants.EPSILON);
		assertEquals(0.2, p.get(1), Constants.EPSILON);
		assertEquals(0.3, p.z(), Constants.EPSILON);
		assertEquals(0.3, p.get(2), Constants.EPSILON);
	}

	@Test
	public void testGetX() {
		Point p = new Point(new Vector(0.1, 0.2, 0.3));
		assertEquals(0.1, p.x(), Constants.EPSILON);
		assertEquals(0.1, p.get(0), Constants.EPSILON);
		assertEquals(0.2, p.y(), Constants.EPSILON);
		assertEquals(0.2, p.get(1), Constants.EPSILON);
		assertEquals(0.3, p.z(), Constants.EPSILON);
		assertEquals(0.3, p.get(2), Constants.EPSILON);
	}

	@Test
	public void testGetY() {
		Point p = new Point(new Vector(0.1, 0.2, 0.3));
		assertEquals(0.1, p.x(), Constants.EPSILON);
		assertEquals(0.1, p.get(0), Constants.EPSILON);
		assertEquals(0.2, p.y(), Constants.EPSILON);
		assertEquals(0.2, p.get(1), Constants.EPSILON);
		assertEquals(0.3, p.z(), Constants.EPSILON);
		assertEquals(0.3, p.get(2), Constants.EPSILON);
	}

	@Test
	public void testGetZ() {
		Point p = new Point(new Vector(0.1, 0.2, 0.3));
		assertEquals(0.1, p.x(), Constants.EPSILON);
		assertEquals(0.1, p.get(0), Constants.EPSILON);
		assertEquals(0.2, p.y(), Constants.EPSILON);
		assertEquals(0.2, p.get(1), Constants.EPSILON);
		assertEquals(0.3, p.z(), Constants.EPSILON);
		assertEquals(0.3, p.get(2), Constants.EPSILON);
	}

	@Test
	public void testSetCoord() {
		Point p = new Point(new Vector(0.1, 0.2, 0.3));
		p.setCoord(0, 0.4);
		p.setCoord(1, 0.5);
		p.setCoord(2, 0.6);
		assertEquals(0.4, p.x(), Constants.EPSILON);
		assertEquals(0.4, p.get(0), Constants.EPSILON);
		assertEquals(0.5, p.y(), Constants.EPSILON);
		assertEquals(0.5, p.get(1), Constants.EPSILON);
		assertEquals(0.6, p.z(), Constants.EPSILON);
		assertEquals(0.6, p.get(2), Constants.EPSILON);
	}

	@Test
	public void testSet() {
		Point p = new Point(new Vector(0.1, 0.2, 0.3));
		p.set(0, 0.4);
		p.set(1, 0.5);
		p.set(2, 0.6);
		assertEquals(0.4, p.x(), Constants.EPSILON);
		assertEquals(0.4, p.get(0), Constants.EPSILON);
		assertEquals(0.5, p.y(), Constants.EPSILON);
		assertEquals(0.5, p.get(1), Constants.EPSILON);
		assertEquals(0.6, p.z(), Constants.EPSILON);
		assertEquals(0.6, p.get(2), Constants.EPSILON);
	}

	@Test
	public void testSetX() {
		Point p = new Point(new Vector(0.1, 0.2, 0.3));
		p.setX(0.4);
		assertEquals(0.4, p.x(), Constants.EPSILON);
		assertEquals(0.4, p.get(0), Constants.EPSILON);
		assertEquals(0.3, p.get(2), Constants.EPSILON);
	}

	@Test
	public void testSetY() {
		Point p = new Point(new Vector(0.1, 0.2, 0.3));
		p.setY(0.5);
		p.setCoord(0, 0.6);
		assertEquals(0.5, p.y(), Constants.EPSILON);
		assertEquals(0.5, p.get(1), Constants.EPSILON);
	}

	@Test
	public void testSetZ() {
		Point p = new Point(0.1, 0.2, 0.3);
		p.setZ(0.6);
		assertEquals(0.6, p.z(), Constants.EPSILON);
		assertEquals(0.6, p.get(2), Constants.EPSILON);
	}

	@Test
	public void testVectorTo() {
		Point p = new Point(1,1,1);
		Point q = new Point(-1,-1,1);
		assertTrue(new Vector(-2,-2,0).equals(p.vectorTo(q)));
		assertTrue(p.vectorTo(q).multiplyThis(-1).equals(q.vectorTo(p)));
		assertTrue(new Vector(0,0,0).equals(q.vectorTo(q)));
		assertTrue(new Vector(0,0,0).equals(p.vectorTo(p)));
	}

	@Test
	public void testCollinear() {
		//If all points are identical colinear should return true
		Point p = new Point(0,0,0);
		Point q = new Point(0,0,0);
		Point r = new Point(0,0,0);
		assertTrue(Point.collinear(p, q, r));
		assertTrue(Point.collinear(p, p, r));
		assertTrue(Point.collinear(p, p, p));
		
		
		p.setX(1);p.setY(1);p.setZ(2);
		q.setX(4);q.setY(1);q.setZ(2);
		r.setX(-4);r.setY(1);r.setZ(2);
		assertTrue(Point.collinear(p, q, r));
		assertTrue(Point.collinear(r, p, q));
		assertTrue(Point.collinear(q, r, p));
		assertTrue(Point.collinear(p, p, r));
		assertTrue(Point.collinear(p, p, p));
		
		p.setY(1.001);
		assertFalse(Point.collinear(p, q, r));
		assertFalse(Point.collinear(r, p, q));
		assertFalse(Point.collinear(q, r, p));
		assertTrue(Point.collinear(p, p, r));
		assertTrue(Point.collinear(p, p, p));
	}

	@Test
	public void testCoplanar() {
		//If all points are identical coplanar should return true
		Point p = new Point(0,0,1);
		Point q = new Point(0,0,1);
		Point r = new Point(0,0,1);
		Point s = new Point(0,0,1);
		assertTrue(Point.coplanar(p, q, r, s));
		assertTrue(Point.coplanar(p, p, r, s));
		assertTrue(Point.coplanar(p, p, p, s));
		assertTrue(Point.coplanar(p, p, p, p));
		
		//If three points are identical coplanar should return true
		p.setX(0.1);
		assertTrue(Point.coplanar(p, q, r, s));
		assertTrue(Point.coplanar(q, p, r, s));
		assertTrue(Point.coplanar(q, r, p, s));
		assertTrue(Point.coplanar(q, r, s, p));

		//If all points are pairwise identical to another point coplanar should return true
		s.setX(0.1);
		assertTrue(Point.coplanar(p, q, r, s));
		assertTrue(Point.coplanar(q, p, r, s));
		assertTrue(Point.coplanar(q, r, p, s));
		assertTrue(Point.coplanar(q, r, s, p));
		
		//Now for some meaningful true cases
		p = new Point(1,1,1);
		q = new Point(2,2,2);
		r = new Point(3,4,3);
		s = new Point(4,5,4);
		assertTrue(Point.coplanar(p, q, r, s));
		assertTrue(Point.coplanar(s, r, p, q));

		p = new Point(3,	1,		1);
		q = new Point(3.2,	2,		1);
		r = new Point(4,	2.1,	1);
		s = new Point(3.9,	1.01,	1);
		assertTrue(Point.coplanar(p, q, r, s));
		assertTrue(Point.coplanar(s, r, p, q));

		//And finally for some false cases
		p.setZ(1.001);
		assertFalse(Point.coplanar(p, q, r, s));
		assertFalse(Point.coplanar(s, r, p, q));
		q.setZ(1.001);
		assertFalse(Point.coplanar(p, q, r, s));
		assertFalse(Point.coplanar(s, r, p, q));
		
	}

	@Test
	public void testTranslate() {
		Point p = new Point(2,1,1);
		Point q = p.clone();
		p.translateThis(0, 0, 0);
		assertTrue(p.equals(q));
		
		p.translateThis(0.001, 0, 0);
		assertFalse(p.equals(q));
		assertTrue(p.x()==2.001);
		assertTrue(p.y()==1);
		assertTrue(p.z()==1);

		p.translateThis(0.001, 2.001, 3.001);
		assertFalse(p.equals(q));
		assertEquals(2.002, p.x(), Constants.EPSILON);
		assertEquals(3.001, p.y(), Constants.EPSILON);
		assertEquals(4.001, p.z(), Constants.EPSILON);
	}

	@Test
	public void testScale() {
		Point p = new Point(2,1,1);
		Point q = p.clone();
		p.scaleThis(1);
		assertTrue(p.equals(q));

		p.scaleThis(3);
		assertFalse(p.equals(q));
		assertTrue(p.x()==6);
		assertTrue(p.y()==3);
		assertTrue(p.z()==3);

		p.scaleThis(0);
		assertTrue(p.x()==0);
		assertTrue(p.y()==0);
		assertTrue(p.z()==0);

		p.scaleThis(2);
		assertTrue(p.x()==0);
		assertTrue(p.y()==0);
		assertTrue(p.z()==0);
	}

	@Test
	public void testAddThis() {
		Point p = new Point(2,1,1);
		Point q = p.clone();
		p.addThis( new Vector(0, 0, 0) );
		assertTrue(p.equals(q));
		
		Point r = p.addThis( new Vector(0.001, 0, 0) );
		assertFalse(p.equals(q));
		assertTrue(p.equals(r));
		assertTrue(p==r);
		assertTrue(p.x()==2.001);
		assertTrue(p.y()==1);
		assertTrue(p.z()==1);

		p.addThis( new Vector(0.001, 2.001, 3.001) );
		assertFalse(p.equals(q));
		assertEquals(2.002, p.x(), Constants.EPSILON);
		assertEquals(3.001, p.y(), Constants.EPSILON);
		assertEquals(4.001, p.z(), Constants.EPSILON);
	}

	@Test
	public void testAdd() {
		Point p = new Point(2,1,1);
		Point q = p.clone();
		
		Point r = p.add( new Vector(0.001, 0.002, 0.003) );
		assertTrue(p.equals(q));
		assertFalse(p.equals(r));
		assertTrue(p.x()==2);
		assertTrue(p.y()==1);
		assertTrue(p.z()==1);
		assertTrue(r.x()==2.001);
		assertTrue(r.y()==1.002);
		assertTrue(r.z()==1.003);
	}

	@Test
	public void testSubtractThis() {
		Point p = new Point(2,1,1);
		Point q = p.clone();
		p.subtractThis( new Vector(0, 0, 0) );
		assertTrue(p.equals(q));
		
		Point r = p.subtractThis( new Vector(0.001, 0, 0) );
		assertFalse(p.equals(q));
		assertTrue(p.equals(r));
		assertTrue(p==r);
		assertTrue(p.x()==1.999);
		assertTrue(p.y()==1);
		assertTrue(p.z()==1);

		p.subtractThis( new Vector(0.001, 2.001, 3.001) );
		assertFalse(p.equals(q));
		assertEquals( 1.998, p.x(), Constants.EPSILON);
		assertEquals(-1.001, p.y(), Constants.EPSILON);
		assertEquals(-2.001, p.z(), Constants.EPSILON);
	}

	@Test
	public void testSubtract() {
		Point p = new Point(2,1,1);
		Point q = p.clone();
		
		Point r = p.subtract( new Vector(0.001, 0.002, 0.003) );
		assertTrue(p.equals(q));
		assertFalse(p.equals(r));
		assertTrue(p.x()==2);
		assertTrue(p.y()==1);
		assertTrue(p.z()==1);
		assertTrue(r.x()==1.999);
		assertTrue(r.y()==0.998);
		assertTrue(r.z()==0.997);
	}

	@Test
	public void testReflectThroughOrigoThis() {
		Point p = new Point(3,2,1);
		Point r = p.reflectThroughOrigoThis();
		assertTrue(p.equals(r));
		assertTrue(p==r);
		assertTrue(p.equals(new Point(-3,-2,-1)));
	}

	@Test
	public void testGetDistanceSquared() {
		Point p = new Point(3,2,1);
		Point q = new Point(-3,1,1);
		assertEquals(37, p.distanceSquared(q), Constants.EPSILON);
		assertEquals(0, p.distanceSquared(p), Constants.EPSILON);
	}

	@Test
	public void testGetDistance() {
		Point p = new Point(3,2,1);
		Point q = new Point(-3,1,1);
		assertEquals(Math.sqrt(37),	p.distance(q), Constants.EPSILON);
		assertEquals(0,				p.distance(p), Constants.EPSILON);
	}

	@Test
	public void testGetBisector() {
		Point p = new Point(3,2,1);
		Point q = new Point(-3,2,1);
		Plane bisector = Point.getBisector(p, q);
		assertEquals(bisector.getDistance(p), bisector.getDistance(q), Constants.EPSILON);
		assertTrue(	
				p.vectorTo(q).normalizeThis().equals(bisector.getNormal().normalizeThis()) || 
				p.vectorTo(q).normalizeThis().equals(bisector.getNormal().normalizeThis().multiplyThis(-1)) );	
		

		p = new Point(1,1,1);
		bisector = Point.getBisector(p, q);
		assertEquals(bisector.getDistance(p),bisector.getDistance(q), Constants.EPSILON);
		assertTrue(	
				p.vectorTo(q).normalizeThis().equals(bisector.getNormal().normalizeThis()) || 
				p.vectorTo(q).normalizeThis().equals(bisector.getNormal().normalizeThis().multiplyThis(-1)) );
	}

	@Test
	public void testMidpoint() {
		Point p = new Point(3,2,1);
		Point q = new Point(0,2,1);
		Point midpoint = Point.getMidpoint(p,q);
		assertTrue(midpoint.distanceSquared(p)==midpoint.distanceSquared(q));
		assertTrue(	Point.collinear(p, q, midpoint) );
		
		q = p.clone();
		assertTrue(Point.getMidpoint(p, q).equals(p));
	}

	@Test
	public void testGetAngle() {
		//First a simple right rectangle in the xy-plane
		Point p = new Point(1,0,0);
		Point q = new Point(0,0,0);
		Point r = new Point(0,1,0);
		assertEquals(Math.PI/2,	Point.getAngle(p, q, r), Constants.EPSILON);
		assertEquals(Math.PI/4,	Point.getAngle(q, r, p), Constants.EPSILON);
		assertEquals(Math.PI/4,	Point.getAngle(r, p, q), Constants.EPSILON);
		
		//Collinear points
		p = new Point(1,1,1);
		q = new Point(2,2,2);
		r = new Point(4,4,4);
		assertEquals(Math.PI, 	Point.getAngle(p, q, r), Constants.EPSILON);
		assertEquals(0,			Point.getAngle(q, r, p), Constants.EPSILON);
		assertEquals(0,			Point.getAngle(r, p, q), Constants.EPSILON);
		
	}

	@Test
	public void testGetDihedralAngle() {
		//First, axis-aligned directions are chosen
		Point p1 = new Point(0,1,0);
		Point p2 = new Point(0,0,0);
		Point p3 = new Point(1,0,0);
		Point p4 = new Point(1,1,0);
		assertEquals(0,		Point.getDihedralAngle(p1, p2, p3, p4), Constants.EPSILON);
		
		p4 = new Point(1,0,1);
		assertEquals( Math.PI/2,	Point.getDihedralAngle(p1, p2, p3, p4), Constants.EPSILON);
		p4 = new Point(1,-1,0);
		assertEquals(Math.PI,		Point.getDihedralAngle(p1, p2, p3, p4), Constants.EPSILON);
		p4 = new Point(1,0,-1);
		assertEquals(-Math.PI/2,	Point.getDihedralAngle(p1, p2, p3, p4), Constants.EPSILON);
		
		//Non-axis-aligned directions
		p1 = new Point(-1,1,0);
		p4 = new Point(2,1,0);
		assertEquals(0,		Point.getDihedralAngle(p1, p2, p3, p4), Constants.EPSILON);
		
		p4 = new Point(2,0,1);
		assertEquals( Math.PI/2,	Point.getDihedralAngle(p1, p2, p3, p4), Constants.EPSILON);
		p4 = new Point(2,-1,0);
		assertEquals( Math.PI,		Point.getDihedralAngle(p1, p2, p3, p4), Constants.EPSILON);
		p4 = new Point(2,0,-1);
		assertEquals(-Math.PI/2,	Point.getDihedralAngle(p1, p2, p3, p4), Constants.EPSILON);
		
		
	}

	@Test
	public void testDominatesPoint3d() {
		Point p1 = new Point(2,1,0);
		Point p2 = new Point(1,0,0);
		assertTrue(p1.dominates(p2));
		assertFalse(p2.dominates(p1));
		
		p2.setX(1.99999);
		assertTrue(p1.dominates(p2));
		assertFalse(p2.dominates(p1));

		p2.setX(2);
		assertTrue(p1.dominates(p2));
		assertFalse(p2.dominates(p1));

		p2.setY(2);
		assertFalse(p1.dominates(p2));
		assertTrue(p2.dominates(p1));

		p2.setY(1);
		assertFalse(p1.dominates(p2));
		assertFalse(p2.dominates(p1));

		p2.setZ(-0.0001);
		assertTrue(p1.dominates(p2));
		assertFalse(p2.dominates(p1));
	}

	@Test
	public void testDominatesPoint3dIntIntInt() {
		Point p1 = new Point(2,	3,	0.02);
		Point p2 = new Point(1,	3,	0.020001);

		assertTrue(p1.dominates(p2,0,1,2));
		assertFalse(p2.dominates(p1,0,1,2));
		assertTrue(p1.dominates(p2,0,2,1));
		assertFalse(p2.dominates(p1,0,2,1));

		assertTrue(p1.dominates(p2,1,0,2));
		assertFalse(p2.dominates(p1,1,0,2));
		assertFalse(p1.dominates(p2,1,2,0));
		assertTrue(p2.dominates(p1,1,2,0));

		assertFalse(p1.dominates(p2,2,0,1));
		assertTrue(p2.dominates(p1,2,0,1));
		assertFalse(p1.dominates(p2,2,1,0));
		assertTrue(p2.dominates(p1,2,1,0));
	}

	@Test
	public void testEqualsPoint3d() {
		Point p1 = new Point(1.1, 2.2, 3.3);
		Point p2 = p1;
		assertTrue(p1.equals(p2));
		p2 = new Point(1.1,2.2,3.3);
		assertTrue(p1.equals(p2));
		p2 = p1.clone();
		assertTrue(p1.equals(p2));

		p2 = new Point(1.101, 2.2, 3.3);
		assertFalse(p1.equals(p2));
		p2 = new Point(1.1, 2.202, 3.3);
		assertFalse(p1.equals(p2));
		p2 = new Point(1.1, 2.2, 3.303);
		assertFalse(p1.equals(p2));

		p2 = new Point(15,16,17);
		assertFalse(p1.equals(p2));
	}

	@Test
	public void testClone() {
		Point p1 = new Point(10, 20, 3.1415);
		Point p2 = p1.clone();
		assertTrue(p1.equals(p2));
		assertFalse(p1==p2);
	}

	@Test
	public void testToVector() {
		Point p = new Point(10,20,3.1415);
		Vector v = p.toVector();
		assertTrue(v.equals(new Vector(10,20,3.1415)));
	}

}
