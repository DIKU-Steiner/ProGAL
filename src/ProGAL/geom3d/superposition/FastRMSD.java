package ProGAL.geom3d.superposition;

import java.util.List;

import ProGAL.geom3d.Point;


/**
 * @author R. Fonseca
 */
public class FastRMSD {
	
	public static double getRMSD(List<Point> points1, List<Point> points2){
		
		double rmsd = 0;
		int n=Math.min(points1.size(), points2.size());
		for(int i=0;i<n;i++){
			for(int j=0;j<i;j++){
				double dist1 = points1.get(i).distance(points1.get(j));
				double dist2 = points2.get(i).distance(points2.get(j));
				rmsd+=(dist1-dist2)*(dist1-dist2);
			}
		}
		return Math.sqrt(rmsd/(n*(n-1)/2));
	}

}
