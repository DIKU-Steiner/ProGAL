package ProGAL.geom3d.superposition;

import java.util.Iterator;
import java.util.List;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.Vector;
import ProGAL.math.Matrix;


/**
 * @author R. Fonseca
 */
public class RMSD {
	public static Transform optimalSuperposition(List<Point> points, List<Point> superposeTo){
		assert(points.size()==superposeTo.size());
		
		int n = points.size();
		double[][] R = new double[3][3], U = new double[3][3];
		double[] mov_com = new double[3], mov_to_ref = new double[3];
		
		double[][] ref_xlist = new double[n][3], mov_xlist = new double[n][3];
		for(int i=0;i<n;i++){
			Point p = points.get(i);
			for(int d=0;d<3;d++) mov_xlist[i][d] = p.get(d);
			p = superposeTo.get(i);
			for(int d=0;d<3;d++) ref_xlist[i][d] = p.get(d);
		}

		double Eo = Superposition.setup_rotation(ref_xlist, mov_xlist, n, mov_com, mov_to_ref, R);
		@SuppressWarnings("unused")
		double residual = Superposition.calculate_rotation_matrix(R, U, Eo);

		Vector preTrans = new Vector(mov_com).multiply(-1);
		Vector postTrans = new Vector(mov_com).add(new Vector(mov_to_ref));
		return new Transform(new Matrix(U), preTrans, postTrans);
	}

	public static Transform optimalSuperposition(List<Point> points, List<Point> superposeTo, List<Integer> indices){
		assert(points.size()==superposeTo.size());
		
		int n = indices.size();
		double[][] R = new double[3][3], U = new double[3][3];
		double[] mov_com = new double[3], mov_to_ref = new double[3];
		
		double[][] ref_xlist = new double[n][3], mov_xlist = new double[n][3];
		for(int i=0;i<n;i++){
			Point p = points.get(indices.get(i));
			for(int d=0;d<3;d++) mov_xlist[i][d] = p.get(d);
			p = superposeTo.get(indices.get(i));
			for(int d=0;d<3;d++) ref_xlist[i][d] = p.get(d);
		}
		
		double Eo = Superposition.setup_rotation(ref_xlist, mov_xlist, n, mov_com, mov_to_ref, R);
		@SuppressWarnings("unused")
		double residual = Superposition.calculate_rotation_matrix(R, U, Eo);

		Vector preTrans = new Vector(mov_com).multiply(-1);
		Vector postTrans = new Vector(mov_com).add(new Vector(mov_to_ref));
		return new Transform(new Matrix(U), preTrans, postTrans);
	}
	
	public static double getRMSD(List<Point> points1, List<Point> points2){
		List<Point> superposed = optimalSuperposition(points1, points2).transform(points1);
		double rmsd = 0;
		int n=0;
		Iterator<Point> it1 = superposed.iterator();
		Iterator<Point> it2 = points2.iterator();
		while(it1.hasNext() && it2.hasNext()){
			Point p1 = it1.next();
			Point p2 = it2.next();
			rmsd+=p1.distanceSquared(p2);
			n++;
		}
		return Math.sqrt(rmsd/n);
	}

}
