package ProGAL.geom3d.volumes.tests;

import static org.junit.Assert.*;

import org.junit.Test;

import ProGAL.geom3d.Circle;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.Vector;
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
	public void testGetIntersectionSphereSphere() {
		Sphere s1, s2; 
		Circle i;

		//Intersection circle left of s1
		s1 = new Sphere(new Point(  0,1,0), 3);
		s2 = new Sphere(new Point(2.5,1,0), 4);
		i = Sphere.getIntersection(s1, s2);
		assertEquals(3.0, i.getCenter().add(new Vector(0,i.getRadius(),0)).distance(s1.getCenter()), 0.0000001 ); //Random point on circle to s1.center
		assertEquals(4.0, i.getCenter().add(new Vector(0,i.getRadius(),0)).distance(s2.getCenter()), 0.0000001 ); //Random point on circle to s2.center
		assertEquals( 1.0, Math.abs(new Vector(1,0,0).dot(i.getNormal())), 0.000001); //Correct normal vector
		assertEquals( 1.0, i.getCenter().y(), 0.000001); //Correct center.y
		assertEquals( 0.0, i.getCenter().z(), 0.000001); //Correct center.z

		//Intersection circle right of s2
		s1 = new Sphere(new Point(  0,1,0), 4);
		s2 = new Sphere(new Point(2.5,1,0), 3);
		i = Sphere.getIntersection(s1, s2);
		assertEquals(4.0, i.getCenter().add(new Vector(0,i.getRadius(),0)).distance(s1.getCenter()), 0.0000001 ); //Random point on circle to s1.center
		assertEquals(3.0, i.getCenter().add(new Vector(0,i.getRadius(),0)).distance(s2.getCenter()), 0.0000001 ); //Random point on circle to s2.center
		assertEquals( 1.0, Math.abs(new Vector(1,0,0).dot(i.getNormal())), 0.000001); //Correct normal vector
		assertEquals( 1.0, i.getCenter().y(), 0.000001); //Correct center.y
		assertEquals( 0.0, i.getCenter().z(), 0.000001); //Correct center.z

		//Intersection circle between s1 and s2
		s1 = new Sphere(new Point(  0,1,0), 3);
		s2 = new Sphere(new Point(2.5,1,0), 3);
		i = Sphere.getIntersection(s1, s2);
		assertEquals(3.0, i.getCenter().add(new Vector(0,i.getRadius(),0)).distance(s1.getCenter()), 0.0000001 ); //Random point on circle to s1.center
		assertEquals(3.0, i.getCenter().add(new Vector(0,i.getRadius(),0)).distance(s2.getCenter()), 0.0000001 ); //Random point on circle to s2.center
		assertEquals( 1.0, Math.abs(new Vector(1,0,0).dot(i.getNormal())), 0.000001); //Correct normal vector
		assertEquals( 1.0, i.getCenter().y(), 0.000001); //Correct center.y
		assertEquals( 0.0, i.getCenter().z(), 0.000001); //Correct center.z

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
		Sphere s1,s2,s3;
		Point[] is;
		Point e1, e2;
		
		s1 = new Sphere(new Point(7.32,22.48, 13.53),22.227305);
		s2 = new Sphere(new Point(29.64,53.66,13.14),58.307177);
		s3 = new Sphere(new Point(54.74,18.42,39.94),59.382983);
		e1 = new Point(5.96,0.41,11.26);
		e2 = new Point(-1.47,5.91,25.45);
		is = Sphere.getIntersections(s1, s2, s3);
		System.out.println(is[0]);
		System.out.println(is[1]);
		double diff1 = Math.abs(e1.distance(is[0])) + Math.abs(e2.distance(is[1]));
		double diff2 = Math.abs(e1.distance(is[1])) + Math.abs(e2.distance(is[0]));
		assertEquals(0.0, Math.min(diff1, diff2), 0.01);
	}

}
