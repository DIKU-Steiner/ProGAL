package ProGAL.geom3d.tessellation.BowyerWatson;

import ProGAL.geom3d.PointWeighted;
import ProGAL.geom3d.volumes.Sphere;

public class Tet {
	public final Corner[] corners;
	public final Corner[] opposites;
	public Sphere circumSphere;
	
	/**
	 *  cornerSides[i] indicates the sign of orient(corners[i+1%4],corners[i+2%4],corners[i+3%4],corners[i]).
	 *  cornerSides[2] is, e.g., the sign of orient(corners[3],corners[0],corners[1],corners[2]).
	 */
	public final int[] cornerSides;
	
	Tet(PointWeighted p0, PointWeighted p1, PointWeighted p2, PointWeighted p3){
		corners = new Corner[]{
				new Corner(this,p0),
				new Corner(this,p1),
				new Corner(this,p2),
				new Corner(this,p3) };
		this.opposites = new Corner[4];
		
		if(orient3d(p1.getCoords(),p2.getCoords(),p3.getCoords(),p0.getCoords())>0)
			cornerSides = new int[]{1,-1,1,-1};
		else
			cornerSides = new int[]{-1,1,-1,1};
		
	}
	
	
	
	double orient3d(double[] pa, double[] pb, double[] pc, double[] pd)	{
		double adx = pa[0] - pd[0];
		double bdx = pb[0] - pd[0];
		double cdx = pc[0] - pd[0];
		double ady = pa[1] - pd[1];
		double bdy = pb[1] - pd[1];
		double cdy = pc[1] - pd[1];
		double adz = pa[2] - pd[2];
		double bdz = pb[2] - pd[2];
		double cdz = pc[2] - pd[2];

		double bdxcdy = bdx * cdy;
		double cdxbdy = cdx * bdy;

		double cdxady = cdx * ady;
		double adxcdy = adx * cdy;

		double adxbdy = adx * bdy;
		double bdxady = bdx * ady;

		return adz * (bdxcdy - cdxbdy)
				+ bdz * (cdxady - adxcdy)
				+ cdz * (adxbdy - bdxady);
	}
}
