package ProGAL.geom3d.predicates.tests;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.Test;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.Triangle;
import ProGAL.geom3d.predicates.ExactJavaPredicates;
import ProGAL.geom3d.predicates.Predicates;
import ProGAL.geom3d.predicates.Predicates.SphereConfig;
import ProGAL.geom3d.volumes.Tetrahedron;
import ProGAL.math.Constants;

public class ExactJavaPredicatesTest {

	Predicates pred;
	
	@Before
	public void setUp() throws Exception {
		this.pred = new ExactJavaPredicates();
	}
	
	@Test
	public void testCircumradiusPointPointPointPoint() {
		Point p0 = new Point(1,1,1);
		Point p1 = new Point(2,1,1);
		Point p2 = new Point(1,2,1);
		Point p3 = new Point(1,1,2);
		
		Tetrahedron tet = new Tetrahedron(p0,p1,p2,p3);
		assertEquals(tet.circumRadius(), pred.circumradius(p0, p1, p2, p3), Constants.EPSILON);
		assertEquals(tet.circumRadius(), pred.circumradius(p0, p1, p3, p2), Constants.EPSILON);
		assertEquals(tet.circumRadius(), pred.circumradius(p0, p2, p1, p3), Constants.EPSILON);
		assertEquals(tet.circumRadius(), pred.circumradius(p0, p3, p1, p2), Constants.EPSILON);
		assertEquals(tet.circumRadius(), pred.circumradius(p1, p0, p2, p3), Constants.EPSILON);
		assertEquals(tet.circumRadius(), pred.circumradius(p1, p0, p3, p2), Constants.EPSILON);
		assertEquals(tet.circumRadius(), pred.circumradius(p1, p2, p0, p3), Constants.EPSILON);
	}

	@Test
	public void testCircumradiusPointPointPoint() {
		Point p0 = new Point(1,1,1);
		Point p1 = new Point(2,1,1);
		Point p2 = new Point(1,2,1);
		
		Triangle t = new Triangle(p0,p1,p2);
		assertEquals(t.circumradius(), pred.circumradius(p0, p1, p2), Constants.EPSILON);
		assertEquals(t.circumradius(), pred.circumradius(p0, p2, p1), Constants.EPSILON);
		assertEquals(t.circumradius(), pred.circumradius(p1, p0, p2), Constants.EPSILON);
		assertEquals(t.circumradius(), pred.circumradius(p2, p0, p1), Constants.EPSILON);
		assertEquals(t.circumradius(), pred.circumradius(p1, p2, p0), Constants.EPSILON);
	}

	@Test
	public void testOrient() {
		Point p0 = new Point(1,1,1);
		Point p1 = new Point(2,1,1);
		Point p2 = new Point(1,2,1);
		Point p3 = new Point(1,1,2);
		Point p4 = new Point(2,2,1);

		/* From http://www.cs.berkeley.edu/~jrs/papers/robustr.pdf:
		 * orient(a,b,c,d) returns a positive value if d lies below the oriented plane passing through
		 * a, b and c. By oriented plane, I mean that a, b and c appear in counterclockwise order 
		 * when viewed from above the plane. (One can apply a left-hand rule: orient your left hand 
		 * with fingers curled to follow the circular sequence abc. If your thumb points toward d 
		 * orient returns a positive value). 
		 */
		assertTrue(pred.orient(p0, p1, p2, p3)<0);
		assertTrue(pred.orient(p0, p2, p1, p3)>0);
		assertTrue(pred.orient(p0, p1, p2, p4)==0);
	}

	@Test
	public void testInspherePointPointPointPointPoint() {
		Point p0 = new Point(1,1,1);
		Point p1 = new Point(2,1,1);
		Point p2 = new Point(1,2,1);
		Point p3 = new Point(1,1,2);
		Point p4 = new Point(2,2,2);
		
		assertEquals(SphereConfig.ON, pred.insphere(p0, p1, p2, p3, p0));
		assertEquals(SphereConfig.ON, pred.insphere(p0, p1, p2, p3, p1));
		assertEquals(SphereConfig.ON, pred.insphere(p0, p1, p2, p3, p2));
		assertEquals(SphereConfig.ON, pred.insphere(p0, p1, p2, p3, p3));
		assertEquals(SphereConfig.ON, pred.insphere(p0, p1, p2, p3, p4));

		p4 = new Tetrahedron(p0,p1,p2,p3).circumCenter();
		assertEquals(SphereConfig.INSIDE, pred.insphere(p0, p1, p2, p3, p4));
		p4 = new Point(1.9999,1.9999,1.9999);
		assertEquals(SphereConfig.INSIDE, pred.insphere(p0, p1, p2, p3, p4));

		p4 = new Point(0.9999, 0.9999, 0.9999);
		assertEquals(SphereConfig.OUTSIDE, pred.insphere(p0, p1, p2, p3, p4));
		p4 = new Point(2,2,2.01);
		assertEquals(SphereConfig.OUTSIDE, pred.insphere(p0, p1, p2, p3, p4));
		
		p1 = new Point(2,2,2);
		p2 = new Point(3,3,0);
		p3 = new Point(3,3,1);
		p4 = new Point(-1,-1,2);
		assertEquals(SphereConfig.COPLANAR, pred.insphere(p0, p1, p2, p3, p4));
		p4 = new Point(1,0,0);
		assertEquals(SphereConfig.COPLANAR, pred.insphere(p0, p1, p2, p3, p4));
	}

	@Test
	public void testInspherePointPointPointPoint() {
		Point p0 = new Point(1,1,1);
		Point p1 = new Point(2,1,1);
		Point p2 = new Point(1,2,1);
		
		Triangle t = new Triangle(p0,p1,p2);
		
		assertEquals(SphereConfig.INSIDE, pred.insphere(p0, p1, p2, t.circumcenter()));
		assertEquals(SphereConfig.INSIDE, pred.insphere(p1, p0, p2, t.circumcenter()));
		assertEquals(SphereConfig.INSIDE, pred.insphere(p2, p1, p0, t.circumcenter()));
		assertEquals(SphereConfig.INSIDE, pred.insphere(p0, p1, p2, new Point(1.0001,1,1)));
		assertEquals(SphereConfig.INSIDE, pred.insphere(p0, p2, p1, new Point(1.0001,1,1)));
		assertEquals(SphereConfig.INSIDE, pred.insphere(p1, p2, p0, new Point(1.0001,1,1)));
		assertEquals(SphereConfig.OUTSIDE, pred.insphere(p0, p1, p2, new Point(0.9999,1,1)));
		assertEquals(SphereConfig.OUTSIDE, pred.insphere(p0, p2, p1, new Point(0.9999,1,1)));
		assertEquals(SphereConfig.OUTSIDE, pred.insphere(p2, p1, p0, new Point(0.9999,1,1)));
		assertEquals(SphereConfig.OUTSIDE, pred.insphere(p2, p0, p1, new Point(0.9999,1,1)));
		assertEquals(SphereConfig.OUTSIDE, pred.insphere(p2, p0, p1, new Point(2.001,2,1)));
		assertEquals(SphereConfig.ON, pred.insphere(p2, p0, p1, p0));
		assertEquals(SphereConfig.ON, pred.insphere(p2, p0, p1, p1));
		assertEquals(SphereConfig.ON, pred.insphere(p2, p0, p1, p2));
		assertEquals(SphereConfig.ON, pred.insphere(p2, p0, p1, new Point(2,2,1)));
	}

	@Test
	public void testDiffsides() {
		fail("Not yet implemented");
	}

	@Test
	public void testInplane() {
		fail("Not yet implemented");
	}

	@Test
	public void testEdgeinsphere() {
		fail("Not yet implemented");
	}

	@Test
	public void testEdgecircumradius() {
		fail("Not yet implemented");
	}

}
