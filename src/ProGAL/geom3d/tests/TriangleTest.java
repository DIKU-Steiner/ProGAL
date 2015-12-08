package ProGAL.geom3d.tests;

import static org.junit.Assert.*;

import org.junit.Test;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.Triangle;
import ProGAL.geom3d.Vector;
import ProGAL.math.Constants;

public class TriangleTest {

	@Test
	public void testGetP1() {
		Triangle t = new Triangle(new Point(1,1,1),new Point(2,1,1), new Point(1,2,1));
		assertTrue(t.getP1().equals(new Point(1,1,1)));
		t.getP1().setX(-1);
		assertTrue(t.getP1().equals(new Point(-1,1,1)));
	}

	@Test
	public void testGetP2() {
		Triangle t = new Triangle(new Point(1,1,1),new Point(2,1,1), new Point(1,2,1));
		assertTrue(t.getP2().equals(new Point(2,1,1)));
		t.getP2().setX(-1);
		assertTrue(t.getP2().equals(new Point(-1,1,1)));
	}

	@Test
	public void testGetP3() {
		Triangle t = new Triangle(new Point(1,1,1),new Point(2,1,1), new Point(1,2,1));
		assertTrue(t.getP3().equals(new Point(1,2,1)));
		t.getP3().setX(-1);
		assertTrue(t.getP3().equals(new Point(-1,2,1)));
	}

	@Test
	public void testGetCorner() {
		Triangle t = new Triangle(new Point(1,1,1),new Point(2,1,1), new Point(1,2,1));
		assertTrue(t.getPoint(0).equals(new Point(1,1,1)));
		assertTrue(t.getPoint(1).equals(new Point(2,1,1)));
		assertTrue(t.getPoint(2).equals(new Point(1,2,1)));
		t.getPoint(0).setZ(2);
		t.getPoint(1).setZ(2);
		t.getPoint(2).setZ(2);
		assertTrue(t.getPoint(0).equals(new Point(1,1,2)));
		assertTrue(t.getPoint(1).equals(new Point(2,1,2)));
		assertTrue(t.getPoint(2).equals(new Point(1,2,2)));
	}

	@Test
	public void testGetCenter() {
		Triangle t = new Triangle(new Point(1,1,1),new Point(2,1,1), new Point(1,2,1));
		assertTrue(t.getCenter().distance(t.getPoint(0))<t.getPoint(0).distance(t.getPoint(1)));
		assertTrue(t.getCenter().distance(t.getPoint(0))<t.getPoint(0).distance(t.getPoint(2)));
		assertTrue(t.getCenter().distance(t.getPoint(1))<t.getPoint(1).distance(t.getPoint(2)));
		//Note .. not generally true, but true for this triangle.
	}

	@Test
	public void testGetArea() {
		Triangle t = new Triangle(new Point(1,1,1),new Point(2,1,1), new Point(1,2,1));
		assertEquals(0.5, t.getArea(), Constants.EPSILON);
		t.getP2().setX(4);
		assertEquals(1.5, t.getArea(), Constants.EPSILON);

		//A flat triangle should have area 0
		t = new Triangle(new Point(1,1,1),new Point(2,2,2), new Point(3,3,3));
		assertEquals(0,t.getArea(),Constants.EPSILON);
		t = new Triangle(new Point(1,1,1),new Point(2,2,2), new Point(300,300,300));
		assertEquals(0,t.getArea(),Constants.EPSILON);
	}

	@Test
	public void testGetNormal() {
		Triangle t = new Triangle(new Point(1,1,1),new Point(2,1,1), new Point(1,2,1));
		assertTrue(t.getNormal().equals(new Vector(0,0,1)) || t.getNormal().equals(new Vector(0,0,-1)) );
		
	}

}
