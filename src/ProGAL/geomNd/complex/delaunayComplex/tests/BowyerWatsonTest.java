package ProGAL.geomNd.complex.delaunayComplex.tests;

import static org.junit.Assert.*;

import org.junit.Test;

import ProGAL.geomNd.Point;
import ProGAL.geomNd.complex.Vertex;
import ProGAL.geomNd.complex.Tessel;
import ProGAL.geomNd.complex.delaunayComplex.BowyerWatson;

public class BowyerWatsonTest {

	@Test
	public void testInSphereTesselPointInt() {
		Vertex p0 = new Vertex(new Point(new double[]{1000,0}));
		Vertex p1 = new Vertex(new Point(new double[]{-500,-866}));
		Vertex p2 = new Vertex(new Point(new double[]{0,0}));
		Vertex p = new Vertex(new Point(new double[]{2,3}));
		Tessel t = new Tessel(new Vertex[]{p0,p1,p2});
		//boolean b = BW.inSphere(t, p, 2);
		assertFalse(BowyerWatson.inSphere(t, p));
		Vertex pp = new Vertex(Point.getMidpoint(p2, p1));
		assertTrue(BowyerWatson.inSphere(t, pp));
		assertTrue(BowyerWatson.inSphere(t, p0));
	}

}
