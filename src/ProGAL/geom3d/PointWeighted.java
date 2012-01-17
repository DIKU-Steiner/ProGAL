package ProGAL.geom3d;

import ProGAL.geom3d.volumes.Sphere;

public class PointWeighted extends Point{
	/**
	 *  A weighted point in (x,y,z)-space with weight w represented with double precision. 
	 */

	private static final long serialVersionUID = 1L;   // what is this?
	double w;
	
	/** Construct a weighted point with the specified coordinates. */
	public PointWeighted(double x, double y, double z, double w) { 
		super(new double[]{x,y,z});
		this.w = w;
	}
	/** Construct a weighted point with the specified coordinates. */
	public PointWeighted(double[] coords, double w) { 
		super(coords);
		this.w = w;
	}
	/** Construct a weighted point that is a clone of p. */
	public PointWeighted(PointWeighted p) { 
		super(new double[]{p.coords[0],p.coords[1],p.coords[2]});
		this.w = p.w;
	}
	/** Construct a weighted point with the center of the sphere as its coordinates and with
	 * the square of its radius as weight. */
	public PointWeighted(Sphere s) { 
		super(s.getCenter());
		this.w = s.getRadiusSquared();
	}

	public double getWeight(){ return w; }

}
