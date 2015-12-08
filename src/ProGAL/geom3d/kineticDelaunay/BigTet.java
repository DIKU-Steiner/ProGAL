package ProGAL.geom3d.kineticDelaunay;

import java.util.List;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.volumes.Tetrahedron;
import ProGAL.geomNd.Vector;

public class BigTet extends Tet{
	
	BigTet(List<Point> points){
		super(createCorners(points));		
	}
	
	public BigTet(Tetrahedron tetra) {
		super(tetra);
		
	}

	private static Vertex[] createCorners(List<Point> points){
		int d = 3;
		double[][] minMax = new double[d][2];
		Vector center = new Vector(d);
		double max = 0;
		for(int i=0;i<d;i++) {
			minMax[i][0] = Double.POSITIVE_INFINITY;
			minMax[i][1] = Double.NEGATIVE_INFINITY;
			for(Point p: points){
				if(p.get(i)<minMax[i][0]) minMax[i][0] = p.get(i);
				if(p.get(i)>minMax[i][1]) minMax[i][1] = p.get(i);
			}
			if(minMax[i][1]-minMax[i][0]>max) max = minMax[i][1]-minMax[i][0];
			center.set(i, (minMax[i][0]+minMax[i][1])/2.0);
		}
		max = Math.max(max, 0.3);
		
		//Regular d-dimensional simplex
		double angle = -1.0/d;
		Vertex[] ret = new Vertex[d+1];
		for(int s=0;s<d+1;s++) ret[s] = new BigVertex(new Point(0,0,0));
		double tmp = 0;
		for(int i=0;i<d;i++){
			double nCoord = Math.sqrt(1-tmp); 
			ret[i].set(i, nCoord);
			double rCoord = (angle - tmp)/nCoord;

			for(int s=i+1;s<d+1;s++)
				ret[s].set(i, rCoord);

			tmp += rCoord*rCoord;
		}
		for(int i=0;i<=d;i++){
			ret[i].multiplyThis(max*2);
			ret[i].addThis(center);
		}
		return ret;
	}
	
	private static class BigVertex extends Vertex{
		private static final long serialVersionUID = 1L;

		public BigVertex(Point p) {
			super(p);
		}

	}
}
