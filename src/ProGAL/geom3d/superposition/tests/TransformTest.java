package ProGAL.geom3d.superposition.tests;

import static org.junit.Assert.*;

import org.junit.Test;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.Vector;
import ProGAL.geom3d.superposition.Transform;
import ProGAL.math.Matrix;

public class TransformTest {

	
	
	@Test
	public void testTransformPoint() {
		Matrix m = new Matrix(new double[][]{{0,-1,0},{1,0,0},{0,0,1}});
		Vector pre = new Vector(-1,-2,0);
		Vector post = new Vector(3,0,0);
		Transform trans = new Transform(m, pre,post);

		Point p = new Point(4,1,0);
		assertTrue(  new Point(4,3,0).equals(trans.transform(p))  );
		assertFalse(  new Point(4,3,0).equals(p)  );
		
		trans = new Transform(m,post);
		p = new Point(4,1,0);
		assertTrue(  new Point(2,4,0).equals(trans.transform(p))  );
		assertFalse(  new Point(2,4,0).equals(p)  );
	}

	@Test
	public void testTransformInPoint() {
		Matrix m = new Matrix(new double[][]{{0,-1,0},{1,0,0},{0,0,1}});
		Vector pre = new Vector(-1,-2,0);
		Vector post = new Vector(3,0,0);
		Transform trans = new Transform(m, pre,post);

		Point p = new Point(4,1,0);
		assertTrue(  new Point(4,3,0).equals(trans.transformIn(p))  );
		assertTrue(  new Point(4,3,0).equals(p)  );

		trans = new Transform(m,post);
		p = new Point(4,1,0);
		assertTrue(  new Point(2,4,0).equals(trans.transformIn(p))  );
		assertTrue(  new Point(2,4,0).equals(p)  );
	}

}
