package ProGAL.geom3d.tests;

import static org.junit.Assert.*;

import org.junit.Test;

import ProGAL.geom3d.Point3d;
import ProGAL.geom3d.Vector3d;
import ProGAL.math.Constants;

public class Point3dTest {

	@Test
	public void testPoint3dDoubleDoubleDouble() {
		Point3d p = new Point3d(0.1, 0.2, 0.3);
		assertEquals(0.1, p.getX(), Constants.EPSILON);
		assertEquals(0.1, p.get(0), Constants.EPSILON);
		assertEquals(0.2, p.getY(), Constants.EPSILON);
		assertEquals(0.2, p.get(1), Constants.EPSILON);
		assertEquals(0.3, p.getZ(), Constants.EPSILON);
		assertEquals(0.3, p.get(2), Constants.EPSILON);
	}

	@Test
	public void testPoint3dPoint3d() {
		Point3d p = new Point3d(new Point3d(0.1, 0.2, 0.3));
		assertEquals(0.1, p.getX(), Constants.EPSILON);
		assertEquals(0.1, p.get(0), Constants.EPSILON);
		assertEquals(0.2, p.getY(), Constants.EPSILON);
		assertEquals(0.2, p.get(1), Constants.EPSILON);
		assertEquals(0.3, p.getZ(), Constants.EPSILON);
		assertEquals(0.3, p.get(2), Constants.EPSILON);
	}

	@Test
	public void testPoint3dVector3d() {
		Point3d p = new Point3d(new Vector3d(0.1, 0.2, 0.3));
		assertEquals(0.1, p.getX(), Constants.EPSILON);
		assertEquals(0.1, p.get(0), Constants.EPSILON);
		assertEquals(0.2, p.getY(), Constants.EPSILON);
		assertEquals(0.2, p.get(1), Constants.EPSILON);
		assertEquals(0.3, p.getZ(), Constants.EPSILON);
		assertEquals(0.3, p.get(2), Constants.EPSILON);
	}

	@Test
	public void testGetCoord() {
		Point3d p = new Point3d(new Vector3d(0.1, 0.2, 0.3));
		assertEquals(0.1, p.getX(), Constants.EPSILON);
		assertEquals(0.1, p.get(0), Constants.EPSILON);
		assertEquals(0.2, p.getY(), Constants.EPSILON);
		assertEquals(0.2, p.get(1), Constants.EPSILON);
		assertEquals(0.3, p.getZ(), Constants.EPSILON);
		assertEquals(0.3, p.get(2), Constants.EPSILON);
	}

	@Test
	public void testGet() {
		Point3d p = new Point3d(new Vector3d(0.1, 0.2, 0.3));
		assertEquals(0.1, p.getX(), Constants.EPSILON);
		assertEquals(0.1, p.get(0), Constants.EPSILON);
		assertEquals(0.2, p.getY(), Constants.EPSILON);
		assertEquals(0.2, p.get(1), Constants.EPSILON);
		assertEquals(0.3, p.getZ(), Constants.EPSILON);
		assertEquals(0.3, p.get(2), Constants.EPSILON);
	}

	@Test
	public void testGetX() {
		Point3d p = new Point3d(new Vector3d(0.1, 0.2, 0.3));
		assertEquals(0.1, p.getX(), Constants.EPSILON);
		assertEquals(0.1, p.get(0), Constants.EPSILON);
		assertEquals(0.2, p.getY(), Constants.EPSILON);
		assertEquals(0.2, p.get(1), Constants.EPSILON);
		assertEquals(0.3, p.getZ(), Constants.EPSILON);
		assertEquals(0.3, p.get(2), Constants.EPSILON);
	}

	@Test
	public void testGetY() {
		Point3d p = new Point3d(new Vector3d(0.1, 0.2, 0.3));
		assertEquals(0.1, p.getX(), Constants.EPSILON);
		assertEquals(0.1, p.get(0), Constants.EPSILON);
		assertEquals(0.2, p.getY(), Constants.EPSILON);
		assertEquals(0.2, p.get(1), Constants.EPSILON);
		assertEquals(0.3, p.getZ(), Constants.EPSILON);
		assertEquals(0.3, p.get(2), Constants.EPSILON);
	}

	@Test
	public void testGetZ() {
		Point3d p = new Point3d(new Vector3d(0.1, 0.2, 0.3));
		assertEquals(0.1, p.getX(), Constants.EPSILON);
		assertEquals(0.1, p.get(0), Constants.EPSILON);
		assertEquals(0.2, p.getY(), Constants.EPSILON);
		assertEquals(0.2, p.get(1), Constants.EPSILON);
		assertEquals(0.3, p.getZ(), Constants.EPSILON);
		assertEquals(0.3, p.get(2), Constants.EPSILON);
	}

	@Test
	public void testSetCoord() {
		Point3d p = new Point3d(new Vector3d(0.1, 0.2, 0.3));
		p.setCoord(0, 0.4);
		p.setCoord(0, 0.5);
		p.setCoord(0, 0.6);
		assertEquals(0.4, p.getX(), Constants.EPSILON);
		assertEquals(0.4, p.get(0), Constants.EPSILON);
		assertEquals(0.5, p.getY(), Constants.EPSILON);
		assertEquals(0.5, p.get(1), Constants.EPSILON);
		assertEquals(0.6, p.getZ(), Constants.EPSILON);
		assertEquals(0.6, p.get(2), Constants.EPSILON);
	}

	@Test
	public void testSet() {
		Point3d p = new Point3d(new Vector3d(0.1, 0.2, 0.3));
		p.set(0, 0.4);
		p.set(0, 0.5);
		p.set(0, 0.6);
		assertEquals(0.4, p.getX(), Constants.EPSILON);
		assertEquals(0.4, p.get(0), Constants.EPSILON);
		assertEquals(0.5, p.getY(), Constants.EPSILON);
		assertEquals(0.5, p.get(1), Constants.EPSILON);
		assertEquals(0.6, p.getZ(), Constants.EPSILON);
		assertEquals(0.6, p.get(2), Constants.EPSILON);
	}

	@Test
	public void testSetX() {
		Point3d p = new Point3d(new Vector3d(0.1, 0.2, 0.3));
		p.setX(0.4);
		assertEquals(0.4, p.getX(), Constants.EPSILON);
		assertEquals(0.4, p.get(0), Constants.EPSILON);
		assertEquals(0.6, p.get(2), Constants.EPSILON);
	}

	@Test
	public void testSetY() {
		Point3d p = new Point3d(new Vector3d(0.1, 0.2, 0.3));
		p.setY(0.5);
		p.setCoord(0, 0.6);
		assertEquals(0.5, p.getY(), Constants.EPSILON);
		assertEquals(0.5, p.get(1), Constants.EPSILON);
	}

	@Test
	public void testSetZ() {
		Point3d p = new Point3d(new Vector3d(0.1, 0.2, 0.3));
		p.setZ(0.6);
		assertEquals(0.6, p.getZ(), Constants.EPSILON);
		assertEquals(0.6, p.get(2), Constants.EPSILON);
	}

	@Test
	public void testVectorTo() {
		Point3d p = new Point3d(new Vector3d(1,1,1));
		Point3d q = new Point3d(new Vector3d(-1,-1,1));
		assertTrue(new Vector3d(-2,-2,0).equals(p.vectorTo(q)));
		assertTrue(p.vectorTo(q).multiplyThis(-1).equals(q.vectorTo(p)));
		assertTrue(new Vector3d(0,0,0).equals(q.vectorTo(q)));
		assertTrue(new Vector3d(0,0,0).equals(p.vectorTo(p)));
	}

	@Test
	public void testColinear() {
		//If all points are identical colinear should return true
		Point3d p = new Point3d(new Vector3d(0,0,0));
		Point3d q = new Point3d(new Vector3d(0,0,0));
		Point3d r = new Point3d(new Vector3d(0,0,0));
		assertTrue(Point3d.collinear(p, q, r));
		assertTrue(Point3d.collinear(p, p, r));
		assertTrue(Point3d.collinear(p, p, p));
		
		
		p.setX(1);p.setY(1);p.setZ(2);
		q.setX(4);q.setY(1);q.setZ(2);
		r.setX(-4);r.setY(1);r.setZ(2);
		assertTrue(Point3d.collinear(p, q, r));
		assertTrue(Point3d.collinear(r, p, q));
		assertTrue(Point3d.collinear(q, r, p));
		assertTrue(Point3d.collinear(p, p, r));
		assertTrue(Point3d.collinear(p, p, p));
		
		p.setY(1.001);
		assertFalse(Point3d.collinear(p, q, r));
		assertFalse(Point3d.collinear(r, p, q));
		assertFalse(Point3d.collinear(q, r, p));
		assertTrue(Point3d.collinear(p, p, r));
		assertTrue(Point3d.collinear(p, p, p));
	}

	@Test
	public void testCoplanar() {
		fail("Not yet implemented");
	}

	@Test
	public void testIsBehindPoint3dPoint3dVector3d() {
		fail("Not yet implemented");
	}

	@Test
	public void testIsBehindPoint3dPoint3dPoint3dPoint3d() {
		fail("Not yet implemented");
	}

	@Test
	public void testIsInside() {
		fail("Not yet implemented");
	}

	@Test
	public void testIsInCone() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetNormal() {
		fail("Not yet implemented");
	}

	@Test
	public void testTranslate() {
		fail("Not yet implemented");
	}

	@Test
	public void testScale() {
		fail("Not yet implemented");
	}

	@Test
	public void testAddThis() {
		fail("Not yet implemented");
	}

	@Test
	public void testAdd() {
		fail("Not yet implemented");
	}

	@Test
	public void testSubtractThis() {
		fail("Not yet implemented");
	}

	@Test
	public void testSubtract() {
		fail("Not yet implemented");
	}

	@Test
	public void testReflectThroughOrigo() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetDistanceSquared() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetDistance() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetBisector() {
		fail("Not yet implemented");
	}

	@Test
	public void testMidpoint() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetAngle() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetDihedralAngle() {
		fail("Not yet implemented");
	}

	@Test
	public void testDominatesPoint3d() {
		fail("Not yet implemented");
	}

	@Test
	public void testDominatesPoint3dIntIntInt() {
		fail("Not yet implemented");
	}

	@Test
	public void testEqualsPoint3d() {
		fail("Not yet implemented");
	}

	@Test
	public void testClone() {
		fail("Not yet implemented");
	}

	@Test
	public void testToVector() {
		fail("Not yet implemented");
	}

	@Test
	public void testToDoubleArray() {
		fail("Not yet implemented");
	}

}
