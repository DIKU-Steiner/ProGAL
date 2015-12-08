package ProGAL.math.test;

import static org.junit.Assert.*;

import org.junit.Test;

import ProGAL.geomNd.Vector;
import ProGAL.math.Matrix;

public class MatrixTest {

	@Test
	public void testMatrixIntInt() {
		Matrix m = new Matrix(3, 4);
		assertEquals(3, m.getM());
		assertEquals(4, m.getN());
	}

	@Test
	public void testMatrixDoubleArrayArray() {
		Matrix m = new Matrix(new double[][]{
				new double[]{1, 2, 3, 4},
				new double[]{5, 6, 7, 8},
				new double[]{9,10,11,12},
		});
		assertEquals(3, m.getM());
		assertEquals(4, m.getN());
	}

	@Test
	public void testGetRow() {
		Matrix m = new Matrix(new double[][]{
				new double[]{1, 2, 3, 4},
				new double[]{5, 6, 7, 8},
				new double[]{9,10,11,12},
		});
		assertEquals(new Vector(new double[]{1,2,3,4}), m.getRow(0));
		assertEquals(new Vector(new double[]{5,6,7,8}), m.getRow(1));
		assertEquals(new Vector(new double[]{9,10,11,12}), m.getRow(2));
	}

	@Test
	public void testGetColumn() {
		Matrix m = new Matrix(new double[][]{
				new double[]{1, 2, 3, 4},
				new double[]{5, 6, 7, 8},
				new double[]{9,10,11,12},
		});
		assertEquals(new Vector(new double[]{1,5,9}), m.getColumn(0));
		assertEquals(new Vector(new double[]{2,6,10}), m.getColumn(1));
		assertEquals(new Vector(new double[]{3,7,11}), m.getColumn(2));
		assertEquals(new Vector(new double[]{4,8,12}), m.getColumn(3));
	}

	@Test
	public void testTranspose() {
		Matrix m = new Matrix(new double[][]{
				new double[]{1, 2, 3, 4},
				new double[]{5, 6, 7, 8},
				new double[]{9,10,11,12},
		});
		Matrix mt = m.transpose();
		Matrix t = new Matrix(new double[][]{
				new double[]{1,5,9},
				new double[]{2,6,10},
				new double[]{3,7,11},
				new double[]{4,8,12}
		});
		assertEquals(t,mt);
		assertEquals(m, mt.transpose());
	}

	@Test
	public void testMultiplyInVector() {
		Matrix m = new Matrix(new double[][]{
				new double[]{1, 2, 3, 4},
				new double[]{5, 6, 7, 8},
				new double[]{9,10,11,12},
		});
		Vector v = new Vector(new double[]{1,2,0,1});
		v = m.multiplyIn(v);
		
	}

	@Test
	public void testMultiplyPoint() {
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
		Matrix prod = A.multiply(B);
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
