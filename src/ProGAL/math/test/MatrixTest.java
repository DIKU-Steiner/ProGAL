package ProGAL.math.test;

import static org.junit.Assert.*;

import org.junit.Test;

import ProGAL.math.Matrix;

public class MatrixTest {

	@Test
	public void testMatrixIntInt() {
		fail("Not yet implemented");
	}

	@Test
	public void testMatrixDoubleArrayArray() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetRow() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetColumn() {
		fail("Not yet implemented");
	}

	@Test
	public void testTranspose() {
		fail("Not yet implemented");
	}

	@Test
	public void testApplyToInVector() {
		fail("Not yet implemented");
	}

	@Test
	public void testApplyToPoint() {
		fail("Not yet implemented");
	}

	@Test
	public void testApplyToInPoint() {
		fail("Not yet implemented");
	}

	@Test
	public void testApplyToVector() {
		fail("Not yet implemented");
	}

	@Test
	public void testApplyToThis() {
		fail("Not yet implemented");
	}

	@Test
	public void testApplyToMatrix() {
		Matrix A = new Matrix(new double[][]{
				new double[]{1,0},
				new double[]{1,1},
				new double[]{0,1},
				new double[]{0,0},
				});
		Matrix B = new Matrix(new double[][]{
				new double[]{1,0,1},
				new double[]{1,1,1},
				});
		Matrix prod = A.applyTo(B);
		assertEquals(4, prod.getM());
		assertEquals(3, prod.getN());
		assertEquals(1.0, prod.get(0, 0), 0.000001);
		assertEquals(2.0, prod.get(1, 0), 0.000001);
		assertEquals(1.0, prod.get(2, 0), 0.000001);
		assertEquals(0.0, prod.get(3, 0), 0.000001);

		assertEquals(0.0, prod.get(0, 1), 0.000001);
		assertEquals(1.0, prod.get(1, 1), 0.000001);
		assertEquals(1.0, prod.get(2, 1), 0.000001);
		assertEquals(0.0, prod.get(3, 1), 0.000001);

		assertEquals(1.0, prod.get(0, 2), 0.000001);
		assertEquals(2.0, prod.get(1, 2), 0.000001);
		assertEquals(1.0, prod.get(2, 2), 0.000001);
		assertEquals(0.0, prod.get(3, 2), 0.000001);
		
	}

	@Test
	public void testAdd() {
		fail("Not yet implemented");
	}

	@Test
	public void testAddThis() {
		fail("Not yet implemented");
	}

	@Test
	public void testMultiply() {
		fail("Not yet implemented");
	}

	@Test
	public void testMultiplyThis() {
		fail("Not yet implemented");
	}

	@Test
	public void testDeterminant() {
		fail("Not yet implemented");
	}

	@Test
	public void testMinor() {
		fail("Not yet implemented");
	}

	@Test
	public void testInvert() {
		fail("Not yet implemented");
	}

	@Test
	public void testInvertThis() {
		fail("Not yet implemented");
	}

	@Test
	public void testReduce() {
		fail("Not yet implemented");
	}

	@Test
	public void testReduceThis() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetEigenvalueDecomposition() {
		fail("Not yet implemented");
	}

	@Test
	public void testCreateIdentityMatrix() {
		fail("Not yet implemented");
	}

	@Test
	public void testCreateColumnMatrix() {
		fail("Not yet implemented");
	}

	@Test
	public void testCreate4x4ColumnMatrix() {
		fail("Not yet implemented");
	}

	@Test
	public void testCreateRowMatrix() {
		fail("Not yet implemented");
	}

	@Test
	public void testCreateRotationMatrix() {
		fail("Not yet implemented");
	}

}
