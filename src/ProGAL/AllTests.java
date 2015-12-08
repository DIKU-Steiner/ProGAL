package ProGAL;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({//Add all test-classes here
	ProGAL.geom3d.volumes.tests.LensTest.class,
	ProGAL.geom3d.tests.LineTest.class,
	ProGAL.geom3d.tests.PlaneTest.class,
	ProGAL.geom3d.tests.PointTest.class,
	ProGAL.geom3d.tests.TriangleTest.class,
	ProGAL.math.test.CombinatoricsTest.class,
	ProGAL.math.test.PolynomialTest.class,
	ProGAL.math.test.MatrixTest.class,
	ProGAL.proteins.belta.tests.SecondaryStructureTest.class,
	ProGAL.proteins.belta.tests.SheetTopologyTest.class,
	ProGAL.proteins.belta.tests.BetaTopologyTest.class
	})

/** 
 * A collection of all unit-tests in ProGAL.
 *  
 * @author R.Fonseca
 */
public class AllTests {}
