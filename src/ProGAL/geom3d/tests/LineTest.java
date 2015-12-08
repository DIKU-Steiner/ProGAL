package ProGAL.geom3d.tests;

import static org.junit.Assert.*;

import org.junit.Test;

import ProGAL.geom3d.Line;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.LineSegment;
import ProGAL.geom3d.PointList;
import ProGAL.geom3d.Vector;
import ProGAL.math.Constants;

public class LineTest {

	@Test
	public void testLineVector() {
		Line l = new Line(new Vector(2,1,1));
		assertTrue(l.getP().equals(new Point(0,0,0)));
		assertTrue(l.getDir().equals(new Vector(2,1,1)));
	}

	@Test
	public void testLinePointVector() {
		Line l = new Line(new Point(0,1,0), new Vector(2,1,1));
		assertTrue(l.getP().equals(new Point(0,1,0)));
		assertTrue(l.getDir().equals(new Vector(2,1,1)));
	}

	@Test
	public void testLineSegment() {
		LineSegment ls = new LineSegment(new Point(0,1,0), new Point(2,1,1));
		Line l = new Line(ls);
		assertTrue(l.getP().equals(new Point(0,1,0)));
		assertTrue(l.getDir().equals(new Vector(2,0,1)));

		//Check that changes to the linesegment does not change the constructed line
		ls.getA().setZ(200);
		ls.getB().setZ(4);
		assertTrue(l.getP().equals(new Point(0,1,0)));
		assertTrue(l.getDir().equals(new Vector(2,0,1)));
	}

	@Test
	public void testGetP() {
		Point p = new Point(0,1,0);
		Line l = new Line(p, new Vector(2,1,1));
		assertTrue(p==l.getP());
	}

	@Test
	public void testGetDir() {
		Vector d = new Vector(2,1,1);
		Line l = new Line(new Point(0,1,0), d);
		assertTrue(d==l.getDir());
	}

	@Test
	public void testGetPoint() {
		Point p = new Point(0,1,0);
		Line l = new Line(p, new Vector(2,1,1));
		assertFalse(p==l.getPoint(0));
		assertTrue(p.equals(l.getPoint(0)));
		assertEquals(new Point(2,2,1), 		l.getPoint(1));
		assertEquals(new Point(-2,0,-1), 	l.getPoint(-1));
	}

	@Test
	public void testOrthogonalProjectionPoint() {
		Line l = new Line(new Point(0,1,0), new Vector(2,0,0));
		Point proj = l.orthogonalProjection(new Point(4,1,20));
		assertEquals(new Point(4,1,0), proj);
		proj = l.orthogonalProjection(new Point(0,0,0));
		assertEquals(new Point(0,1,0), proj);
		proj = l.orthogonalProjection(new Point(-3,1,0));
		assertEquals(new Point(-3,1,0), proj);
	}

	@Test
	public void testOrthogonalProjectionParameter() {
		Line l = new Line(new Point(0,1,0), new Vector(2,0,0));
		Point proj = l.orthogonalProjection(new Point(4,1,20));
		assertEquals(new Point(4,1,0), proj);
		proj = l.orthogonalProjection(new Point(0,0,0));
		assertEquals(new Point(0,1,0), proj);
		proj = l.orthogonalProjection(new Point(-3,1,0));
		assertEquals(new Point(-3,1,0), proj);
	}

	@Test
	public void testOrthogonalProjectionPointList() {
		Line l = new Line(new Point(0,1,0), new Vector(2,0,0));
		PointList pl = new PointList();
		pl.add(new Point(4,1,20));
		pl.add(new Point(0,0,0));
		pl.add(new Point(-3,1,0));
		LineSegment proj = l.orthogonalProjection(pl);
		assertEquals(new LineSegment(new Point(-3,1,0),new Point(4,1,0)), proj);

		//Ensure that reversing the points do not reverse the line-segment
		pl = new PointList();
		pl.add(new Point(-3,1,0));
		pl.add(new Point(0,0,0));
		pl.add(new Point(4,1,20));
		proj = l.orthogonalProjection(pl);
		assertEquals(new LineSegment(new Point(-3,1,0),new Point(4,1,0)), proj);
	}

	@Test
	public void testOrthogonalProjectionParameters() {
		Line l = new Line(new Point(0,1,0), new Vector(2,0,0));
		PointList pl = new PointList();
		pl.add(new Point(4,1,20));
		pl.add(new Point(0,0,0));
		pl.add(new Point(-3,1,0));
		double[] proj = l.orthogonalProjectionParameters(pl);
		assertEquals(2,proj.length);
		assertEquals(-3/2d, proj[0], Constants.EPSILON);
		assertEquals(2d, proj[1], Constants.EPSILON);
		

		//Ensure that reversing the points do not reverse the line-segment
		pl = new PointList();
		pl.add(new Point(-3,1,0));
		pl.add(new Point(0,0,0));
		pl.add(new Point(4,1,20));
		proj = l.orthogonalProjectionParameters(pl);
		assertEquals(2,proj.length);
		assertEquals(-3/2d, proj[0], Constants.EPSILON);
		assertEquals(2d, proj[1], Constants.EPSILON);
	}

	@Test
	public void testGetDistanceSquared() {
		Line l = new Line(new Point(0,1,0), new Vector(2,0,0));
		Point proj = l.orthogonalProjection(new Point(4,1,20));
		assertEquals(proj.distanceSquared(new Point(4,1,20)), l.getDistanceSquared(new Point(4,1,20)), Constants.EPSILON);
		assertEquals(0, l.getDistanceSquared(new Point(0,1,0)), Constants.EPSILON);
		assertEquals(0, l.getDistanceSquared(new Point(-3,1,0)), Constants.EPSILON);
	}

	@Test
	public void testGetDistance() {
		Line l = new Line(new Point(0,1,0), new Vector(2,0,0));
		Point proj = l.orthogonalProjection(new Point(4,1,20));
		assertEquals(proj.distance(new Point(4,1,20)), l.getDistance(new Point(4,1,20)), Constants.EPSILON);
		assertEquals(0, l.getDistance(new Point(0,1,0)), Constants.EPSILON);
		assertEquals(0, l.getDistance(new Point(-3,1,0)), Constants.EPSILON);
	}

	@Test
	public void testGetMaxDistanceSquared() {
		Line l = new Line(new Point(0,1,0), new Vector(2,0,0));
		PointList pl = new PointList();
		try{ 
			l.getMaxDistanceSquared(pl); 
			assertTrue(false); 
		}catch(Error err){/*Good. That method should throw an error when pl is empty*/}
		
		pl.add(new Point(-3,1,0));
		double m = l.getMaxDistanceSquared(pl);
		assertEquals(0, m, Constants.EPSILON);
		pl.add(new Point(0,0,0));
		m = l.getMaxDistanceSquared(pl);
		assertEquals(1, m, Constants.EPSILON);
		pl.add(new Point(4,1,20));
		m = l.getMaxDistanceSquared(pl);
		assertEquals(20d*20d, m, Constants.EPSILON);
	}

	@Test
	public void testGetMaxDistance() {
		Line l = new Line(new Point(0,1,0), new Vector(2,0,0));
		PointList pl = new PointList();
		try{ 
			l.getMaxDistance(pl); 
			assertTrue(false); 
		}catch(Error err){/*Good. That method should throw an error when pl is empty*/}
		
		pl.add(new Point(-3,1,0));
		double m = l.getMaxDistance(pl);
		assertEquals(0, m, Constants.EPSILON);
		pl.add(new Point(0,0,0));
		m = l.getMaxDistance(pl);
		assertEquals(1, m, Constants.EPSILON);
		pl.add(new Point(4,1,20));
		m = l.getMaxDistance(pl);
		assertEquals(20d, m, Constants.EPSILON);
	}

	@Test
	public void testGetSquaredDistance() {
		Line l1 = new Line(new Point(0,1,0), new Vector(1,0,0));
		assertEquals(0, l1.getSquaredDistance(l1),Constants.EPSILON);
		Line l2 = new Line(new Point(0,0,2), new Vector(0,1,0));
		assertEquals(4, l1.getSquaredDistance(l2),Constants.EPSILON);
		assertEquals(4, l2.getSquaredDistance(l1),Constants.EPSILON);
		l2 = new Line(new Point(0,0,2), new Vector(10,1,0));
		assertEquals(4, l1.getSquaredDistance(l2),Constants.EPSILON);
		assertEquals(4, l2.getSquaredDistance(l1),Constants.EPSILON);
	}
	
	@Test
	public void testRotate(){
		Line l = new Line(new Point(0,1,0), new Vector(1,0,0));
		Point p = new Point(2,0,0);
		assertTrue(l.rotate(p, Math.PI/2).equals(new Point(2,1,-1)));
		
		//Random non-trivial line and point. Rotate point and rotate it back. 
		l = new Line( new Point(0.2, 0.3, 0.4), new Vector(0.5,0.6,0.7));
		p = new Point(3.9,1.0,1.1);
		Point q = l.rotate(p, 1);//Rotate roughly a sixth of a circle
		assertTrue(l.rotate(q, -1).equals(p));
		assertTrue(l.rotate(q, Math.PI*2-1).equals(p));
	}

}
