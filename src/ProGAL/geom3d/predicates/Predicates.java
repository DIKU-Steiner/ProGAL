package ProGAL.geom3d.predicates;

import ProGAL.geom3d.*;
import ProGAL.geom3d.volumes.Tetrahedron;

public abstract class Predicates {
	public enum SphereConfig { INSIDE, OUTSIDE, ON, COPLANAR };
	public enum PlaneConfig  { SAME, DIFF, COPLANAR };

	public abstract double circumradius(Point p0, Point p1, Point p2, Point p3);
	public abstract double circumradius(Tetrahedron t);
	public abstract double circumradius(Point p0, Point p1, Point p2);
	public abstract double circumradius(Triangle t);
	public abstract double orient(Point p0, Point p1, Point p2, Point q);
	
	public abstract SphereConfig insphere(Point p0, Point p1, Point p2, Point p3, Point q);
	public abstract SphereConfig insphere(Tetrahedron t, Point q);
	public abstract SphereConfig insphere(Point p0, Point p1, Point p2, Point q);
	public abstract SphereConfig insphere(Triangle tri, Point q);

	public abstract PlaneConfig diffsides(Point p0, Point p1, Point p2, Point q0, Point q1);

	public abstract boolean inplane(Point p0, Point p1, Point p2, Point p3);

	public abstract SphereConfig edgeinsphere(LineSegment e, Point q);
	public abstract double edgecircumradius(LineSegment ls);

}