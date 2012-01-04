package ProGAL.geom3d.complex.delaunayComplex;

import java.util.LinkedList;
import java.util.List;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.PointWeighted;
import ProGAL.geom3d.predicates.*;
import ProGAL.math.Randomization;

public class DelaunayComplex extends RegularComplex {

	public DelaunayComplex(List<Point> points) {
		super(points, new ExactJavaPredicates());
	}
	
	public static void main(String[] args){
		int insertedPoints = 10000;
		LinkedList<Point> points= new LinkedList<Point>();
		for(int i=0;i<insertedPoints;i++){
			PointWeighted pw = new PointWeighted( 
					Randomization.randBetween(-10.0, 10.0),
					Randomization.randBetween(-10.0, 10.0),
					Randomization.randBetween(-10.0, 10.0) ,0);
			points.add(pw);
		}
		long start = System.nanoTime();
		DelaunayComplex rc = new DelaunayComplex(points);
		long end = System.nanoTime();
		System.out.printf("%d points took %.3fms\n",insertedPoints,(end-start)/1000000.0);
		
	}
}
