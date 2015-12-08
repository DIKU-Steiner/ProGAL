package ProGAL.geom3d.superposition.tests;

import static org.junit.Assert.*;

import java.util.LinkedList;
import java.util.List;

import org.junit.Test;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.PointList;
import ProGAL.geom3d.Vector;
import ProGAL.geom3d.superposition.RMSD;
import ProGAL.geom3d.superposition.Transform;
import ProGAL.math.Matrix;

public class RMSDTest {

	@Test
	public void testOptimalSuperpositionListOfPointListOfPoint() {
		Matrix m = new Matrix(new double[][]{{0,-1,0},{1,0,0},{0,0,1}});
		Vector pre = new Vector(-1,-2,0);
		Vector post = new Vector(3,0,0);
		Transform trans = new Transform(m, pre,post);

		List<Point> orig = PointList.generatePointsInCube(5);
		List<Point> transformedPoints = trans.transform(orig);
		
		Transform backTrans = RMSD.optimalSuperposition(transformedPoints, orig);
		List<Point> backTransformedPoints = backTrans.transform(transformedPoints);
		
		for(int i=0;i<orig.size();i++){
			assertTrue(orig.get(i).distance(backTransformedPoints.get(i))<0.1);
		}
	
	}

	@Test
	public void testOptimalSuperpositionListOfPointListOfPointListOfInteger() {
		Matrix m = new Matrix(new double[][]{{0,-1,0},{1,0,0},{0,0,1}});
		Vector pre = new Vector(-1,-2,0);
		Vector post = new Vector(3,0,0);
		Transform trans = new Transform(m, pre,post);

		List<Point> orig = new PointList();
		orig.add(new Point(1,0,0));
		orig.add(new Point(1,2,0));
		orig.add(new Point(1,0,3));
		orig.add(new Point(4,0,0));
		orig.add(new Point(5,2,3));
		
		List<Point> transformedPoints = trans.transform(orig);

		//This ensures that the superposition gets bad if indices 2 and 3 are included
		orig.add(2, new Point(0,0,100));
		orig.add(3, new Point(0,0,-100));
		transformedPoints.add(2,new Point(0,100,0));
		transformedPoints.add(3,new Point(0,-100,0));
		List<Integer> indices = new LinkedList<Integer>();
		for(int i=0;i<orig.size();i++)
			if(i!=2 && i!=3) indices.add(i);
		
		Transform backTrans = RMSD.optimalSuperposition(transformedPoints, orig, indices);
		List<Point> backTransformedPoints = backTrans.transform(transformedPoints);
		System.out.println(indices);
		System.out.println(orig);
		System.out.println(backTransformedPoints);
		for(int i=0;i<orig.size();i++){
			if(i==2 || i==3) 
				assertFalse(orig.get(i).distance(backTransformedPoints.get(i))<0.1);
			else
				assertTrue(orig.get(i).distance(backTransformedPoints.get(i))<0.1);
		}
		
	}

	@Test
	public void testGetRMSD() {
		fail("Not yet implemented");
	}

}
