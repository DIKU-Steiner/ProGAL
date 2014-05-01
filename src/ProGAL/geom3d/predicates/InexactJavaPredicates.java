package ProGAL.geom3d.predicates;

import ProGAL.geom3d.*;
import ProGAL.geom3d.volumes.Tetrahedron;
import ProGAL.math.Constants;

public class InexactJavaPredicates extends Predicates{
	public double circumradius(Point p0, Point p1, Point p2, Point p3){				return circumradius(new Tetrahedron(p0,p1,p2,p3)); 	}
	public double circumradius(Tetrahedron t){										return t.circumRadius();							}
	public double circumradius(Point p0, Point p1, Point p2){						return circumradius(new Triangle(p0,p1,p2));		}
	public double circumradius(Triangle t){											return t.circumradius(); 						}
	
	public double orient(Point p0, Point p1, Point p2, Point q){					
		Vector q0 = q.vectorTo(p0);
		Vector q1 = q.vectorTo(p1);
		Vector q2 = q.vectorTo(p2);
		return q0.dot(q1.crossThis(q2));//From [Ericsson 05, p. 33]
	}		
	
	public SphereConfig insphere(Point p0, Point p1, Point p2, Point p3, Point q){	return insphere(new Tetrahedron(p0,p1,p2,p3),q);	}
	public SphereConfig insphere(Tetrahedron t, Point q){
		Point tCenter = t.circumCenter();
		double diff = tCenter.distance(q)-tCenter.distance(t.getCorner(0));
		if(Double.isNaN(diff)) return SphereConfig.COPLANAR;
		if(Math.abs(diff)<Constants.EPSILON) return SphereConfig.ON;
		if(diff <= 0) return SphereConfig.INSIDE;
		return SphereConfig.OUTSIDE;
	}
	public SphereConfig insphere(Point p0, Point p1, Point p2, Point q){			return insphere(new Triangle(p0,p1,p2),q);			}
	public SphereConfig insphere(Triangle tri, Point q){
		Point tCenter = tri.circumcenter();
		double diff = tCenter.distance(q)-tCenter.distance(tri.getCorner(0));
		if(Math.abs(diff)<Constants.EPSILON) return SphereConfig.ON;
		if(diff<0) return SphereConfig.INSIDE;
		return SphereConfig.OUTSIDE;
		
	}

	public PlaneConfig diffsides(Point p0, Point p1, Point p2, Point q0, Point q1){
		Plane p = new Plane(p0,p1,p2);
		int a0 = p.above(q0);
		int a1 = p.above(q1);
		int prod = a0*a1;
		if(prod==0) 	return PlaneConfig.COPLANAR;
		else if(prod>0) return PlaneConfig.SAME;
		else 			return PlaneConfig.DIFF;
	}

	public boolean inplane(Point p0, Point p1, Point p2, Point p3){
		return orient(p0,p1,p2,p3)==0;
	}
	@Override
	public double edgecircumradius(LineSegment ls) {
		// TODO Auto-generated method stub
		return 0;
	}
	@Override
	public SphereConfig edgeinsphere(LineSegment e, Point q) {
		// TODO Auto-generated method stub
		return null;
	}

}