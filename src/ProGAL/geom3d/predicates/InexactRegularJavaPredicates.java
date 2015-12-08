package ProGAL.geom3d.predicates;

import ProGAL.geom3d.*;
import ProGAL.geom3d.predicates.Predicates.PlaneConfig;
import ProGAL.geom3d.predicates.Predicates.SphereConfig;
import ProGAL.geom3d.volumes.Tetrahedron;
import ProGAL.math.Constants;
import ProGAL.math.Matrix;

public class InexactRegularJavaPredicates extends Predicates{
	public double circumradius(Point p0, Point p1, Point p2, Point p3){				return circumradius(new Tetrahedron(p0,p1,p2,p3)); 	}
	public double circumradius(Tetrahedron t){										return t.circumRadius();							}
	public double circumradius(Point p0, Point p1, Point p2){						return circumradius(new Triangle(p0,p1,p2));		}
	public double circumradius(Triangle t){											return t.circumradius(); 						}
	
	public double orient(Point p0, Point p1, Point p2, Point q){					
		Vector q0 = q.vectorTo(p0);
		Vector q1 = q.vectorTo(p1);
		Vector q2 = q.vectorTo(p2);
		return q0.dot(q1.crossThis(q2));//From [Ericsson05, p. 33]
	}		
	
	public SphereConfig insphere(Point p0, Point p1, Point p2, Point p3, Point q){
		double orient = orient(p0,p1,p2,p3);
		if(Math.abs(orient)<Constants.EPSILON) return SphereConfig.COPLANAR;
		
		double w0 = ((PointWeighted)p0).getWeight();
		double w1 = ((PointWeighted)p1).getWeight();
		double w2 = ((PointWeighted)p2).getWeight();
		double w3 = ((PointWeighted)p3).getWeight();
		double wq = ((PointWeighted)q).getWeight();
		double l0 = p0.x()*p0.x() + p0.y()*p0.y() + p0.z()*p0.z();
		double l1 = p1.x()*p1.x() + p1.y()*p1.y() + p1.z()*p1.z();
		double l2 = p2.x()*p2.x() + p2.y()*p2.y() + p2.z()*p2.z();
		double l3 = p3.x()*p3.x() + p3.y()*p3.y() + p3.z()*p3.z();
		double lq = q.x()*q.x() + q.y()*q.y() + q.z()*q.z();
		
		Matrix m = new Matrix(new double[][]{
				{p0.x(), p0.y(), p0.z(), l0-w0, 1},
				{p1.x(), p1.y(), p1.z(), l1-w1, 1},
				{p2.x(), p2.y(), p2.z(), l2-w2, 1},
				{p3.x(), p3.y(), p3.z(), l3-w3, 1},
				{q.x(), q.y(), q.z(), lq-wq, 1}
		});
		double det = m.determinant();

		if(Math.abs(det)<Constants.EPSILON) return SphereConfig.ON;
		if(Math.signum(det)*Math.signum(orient)<0)
			return SphereConfig.OUTSIDE;
		else
			return SphereConfig.INSIDE;
	}
	public SphereConfig insphere(Tetrahedron t, Point q){
		return insphere(t.getCorner(0),t.getCorner(1),t.getCorner(2),t.getCorner(3),q);
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