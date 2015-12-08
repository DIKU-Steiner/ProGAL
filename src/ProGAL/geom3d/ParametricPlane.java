package ProGAL.geom3d;

import ProGAL.math.Matrix3x3;


public class ParametricPlane {
	public Point p;
	public Vector n, v1, v2;
	protected Matrix3x3 projInv;
	
	/** 
	 * 8 + 31 = 39HOps
	 */
	public ParametricPlane(Point p, Vector v1, Vector v2){
		this.p = p;
		this.n = v1.cross(v2).normalizeThis();
		this.v1 = v1;
		this.v2 = v2;
		this.projInv = Matrix3x3.createRowMatrix(v1, v2, n);
	}

	/** Construct a parametric plane with the specified point and normal. It is not specified 
	 * how the orthogonal basis is constructed except that the third basis-vector will be the 
	 * normal vector.
	 */
	public ParametricPlane(Point p, Vector normal){
		this.p = p;
		this.n = normal.normalize();
		this.v1 = new Vector(1.001,0.002,0).crossThis(normal).normalizeThis();
		this.v2 = this.n.cross(v1);
		this.projInv = Matrix3x3.createRowMatrix(v1, v2, n);
	}
	
	/** Projects the point v onto this plane and returns the parameters of 
	 * the projected point (scaling of v1, of v2 and finally the distance or 
	 * scaling along n). 9HOps
	 */
	public double[] projectPoint(Point v){
		Vector x = p.vectorTo(v);
		projInv.multiplyIn(x);
		return new double[]{x.x(), x.y(), x.z()};
	}
	
	/** Returns a parameter setting for the line describing the 
	 * intersection with this plane. Assumes there's a point 
	 * intersection.
	 */
	public double intersectionParameter(Line l){
		//From http://www.thepolygoners.com/tutorials/lineplane/lineplane.html
		return n.dot(l.p.vectorTo(p))/n.dot(l.dir);
	}
	
	public Point getP(){
		return p;
	}
	public Point getP(double[] pars){
		return p.add(v1.multiply(pars[0])).addThis(v2.multiply(pars[1]));
	}
}
