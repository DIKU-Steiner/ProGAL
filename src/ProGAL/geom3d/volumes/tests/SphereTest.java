package ProGAL.geom3d.volumes.tests;

import static org.junit.Assert.*;

import org.junit.Test;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.volumes.Sphere;

public class SphereTest {

	@Test
	public void testSpherePointDouble() {
		Sphere s;
		s = new Sphere(new Point(1,2,3), 4);
		assertEquals(new Point(1,2,3), s.getCenter());
		assertEquals(4.0, s.getRadius(), 0.000001);
	}

	@Test
	public void testSpherePointPointPointPoint() {
		Point p1 = new Point(1,1,2);
		Point p2 = new Point(2,1,2);
		Point p3 = new Point(1,2,2);
		Point p4 = new Point(1,1,3);
		Sphere s = new Sphere(p1,p2,p3,p4);
		assertEquals(new Point(1.5,1.5,2.5), s.getCenter());
		assertEquals(new Point(1.5,1.5,2.5).distance(p1), s.getRadius(), 0.00001);
	}

	@Test
	public void testGetSurfaceArea() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetVolume() {
		fail("Not yet implemented");
	}

	@Test
	public void testIsInsidePoint() {
		fail("Not yet implemented");
	}

	@Test
	public void testIsInsidePointDouble() {
		fail("Not yet implemented");
	}

	@Test
	public void testIsIntersected() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetIntersection() {
		fail("Not yet implemented");
	}

	@Test
	public void testIntersectionParameters() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetIntersectionsCircle() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetIntersectionAngle() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetMinSphereCircle() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetMinSpherePointPoint() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetMinSpherePointPointPoint() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetMinSpherePointPointPointPoint() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetIntersectionsSphereSphereSphere() {
		fail("Not yet implemented");
	}

}
