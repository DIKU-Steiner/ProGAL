package ProGAL;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({//Add all test-classes here
	ProGAL.geom3d.tests.LineTest.class,
	ProGAL.geom3d.tests.PlaneTest.class,
	ProGAL.geom3d.tests.PointTest.class,
	ProGAL.geom3d.tests.TriangleTest.class,
	ProGAL.math.test.CombinatoricsTest.class,
	ProGAL.math.test.PolynomialTest.class
	})
public class AllTests {}
