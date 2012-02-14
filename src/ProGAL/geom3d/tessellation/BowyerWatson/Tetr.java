package ProGAL.geom3d.tessellation.BowyerWatson;

import ProGAL.geom3d.PointWeighted;
import ProGAL.geom3d.volumes.Tetrahedron;

public class Tetr extends Tetrahedron{
	public final PointWeighted[] corners = new PointWeighted[4];
	public final Tetr[] neighbors = new Tetr[4];

	/**
	 *  cornerSides[i] indicates the sign of orient(corners[i+1%4],corners[i+2%4],corners[i+3%4],corners[i]).
	 *  cornerSides[2] is, e.g., the sign of orient(corners[3],corners[0],corners[1],corners[2]).
	 */
	public final int[] cornerSides = new int[4];

	public Tetr(PointWeighted p0, PointWeighted p1, PointWeighted p2, PointWeighted p3){
		super(p0,p1,p2,p3);
		corners[0] = p0;
		corners[1] = p1;
		corners[2] = p2;
		corners[3] = p3;

		if(orient3d(p1.getCoords(),p2.getCoords(),p3.getCoords(),p0.getCoords())>0){
			cornerSides[0] = 1;
			cornerSides[1] = -1;
			cornerSides[2] = 1;
			cornerSides[3] = -1;
		}else{
			cornerSides[0] = -1;
			cornerSides[1] = 1;
			cornerSides[2] = -1;
			cornerSides[3] = 1;
		}
		
	}

	/** Return the index of the specified corner-point and -1 if p is not a corner. */
	public int cornerIdx(PointWeighted p){
		for(int i=0;i<4;i++) if(corners[i]==p) return i;
		return -1;
	}
	
	public String toString(){
		return String.format("Tetr[%s,%s,%s,%s]",corners[0],corners[1],corners[2],corners[3]);
	}

	private static double orient3d(double[] pa, double[] pb, double[] pc, double[] pd)	{
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
