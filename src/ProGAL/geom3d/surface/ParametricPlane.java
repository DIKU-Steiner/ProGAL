package ProGAL.geom3d.surface;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.Vector;

public class ParametricPlane extends ParametricSurface {
	private Point p;
	private Vector[] bases;
	
	public ParametricPlane(Point p, Vector[] bases){
		this.p = p;
		this.bases = bases;
	}
	
	@Override
	public Point getPoint(double u, double v) {
		return p.add(bases[0].multiply(u)).addThis(bases[1].multiply(v));
	}
	
	public ParametricPlane clone(){
		return new ParametricPlane(p.clone(), new Vector[]{bases[0].clone(), bases[1].clone()});
	}


}
