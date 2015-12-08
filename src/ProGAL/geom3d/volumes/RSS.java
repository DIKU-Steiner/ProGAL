package ProGAL.geom3d.volumes;

import java.util.Arrays;

import ProGAL.geom3d.ParametricPlane;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.PointList;
import ProGAL.geom3d.Rectangle;
import ProGAL.geom3d.Vector;
import ProGAL.math.Matrix3x3;

/** 
 * Implementation of a Rectangular Swept Sphere that supports overlap check and construction from 
 * a point-set and two RSS'.
 * @author P.Winter and R.Fonseca
 */
public class RSS implements Volume{
	public Rectangle rectangle;
	public double radius;

	public RSS(Point center, Vector[] bases, double radius){
		this.rectangle = new Rectangle(center, bases);
		this.radius = radius;
	}

	public boolean overlaps(Volume vol) {
		if(vol instanceof RSS) return overlaps((RSS)vol);
		throw new Error("Unimplemented ("+vol.getClass().getSimpleName()+")");
	}

	/**
	 * 674HOps 
	 */
	public boolean overlaps(RSS rss){
		double sqRads= (radius+rss.radius); sqRads = sqRads*sqRads; //1HOp
		double sqCenterDist = rectangle.center.distanceSquared(rss.rectangle.center);//3HOps
		if(sqCenterDist<=sqRads) return true;
		return rectangle.distance(rss.rectangle)<=radius+rss.radius; //670HOps
	}

	public double getVolume() {
		double boxVol = rectangle.bases[0].length()*rectangle.bases[1].length()*radius*8;
		double cylVol = 2*(rectangle.bases[0].length()+rectangle.bases[1].length())*(float)Math.PI*radius*radius;
		double sphereVol = (float)Math.PI*radius*radius*radius*4f/3;
		return boxVol+cylVol+sphereVol;
	}

	public Point getCenter() {
		return rectangle.center;
	}

	public String toString(){
		return String.format("RSS[center:%s,bases[%s,%s],radius:%f]", rectangle.center, rectangle.bases[0], rectangle.bases[1],radius);
	}

	public static RSS createBoundingRSS_covariance(PointList points){
		RSS ret;
		Matrix3x3 covMatr = points.getCovariance();
		Vector[] eigenVecs = covMatr.getEigenvectors();
		Vector tmp;
		if(eigenVecs[0].length()<eigenVecs[1].length()) {tmp = eigenVecs[0]; eigenVecs[0] = eigenVecs[1]; eigenVecs[1] = tmp; }
		if(eigenVecs[0].length()<eigenVecs[2].length()) {tmp = eigenVecs[0]; eigenVecs[0] = eigenVecs[2]; eigenVecs[2] = tmp; }
		if(eigenVecs[1].length()<eigenVecs[2].length()) {tmp = eigenVecs[1]; eigenVecs[1] = eigenVecs[2]; eigenVecs[2] = tmp; }
		
		eigenVecs[0].normalizeThis();
		eigenVecs[1].normalizeThis();
//		eigenVecs[2].normalizeThis();
		eigenVecs[2] = eigenVecs[0].cross(eigenVecs[1]);
//		scene.addShape(new Capsule3d(new Point(0,0,0),eigenVecs[0].toPoint(),0.05),Color.red);
//		scene.addShape(new Capsule3d(new Point(0,0,0),eigenVecs[1].toPoint(),0.05),Color.green);
//		scene.addShape(new Capsule3d(new Point(0,0,0),eigenVecs[2].toPoint(),0.05),Color.blue);
		ret = createBoundingRSSFromBases(eigenVecs, points);

		return ret;
	}

	private static RSS createBoundingRSSFromBases(Vector[] bases, PointList points){
		//Find radius along the third base
		double lowestDot = Double.POSITIVE_INFINITY, highestDot = Double.NEGATIVE_INFINITY;
		for(Point p: points){
			double dot = bases[2].dot(p.toVector());
			if(dot<lowestDot){	lowestDot = dot;	}
			if(dot>highestDot){	highestDot = dot;	}
		}
		double radius = (highestDot-lowestDot)/2;
		if(radius<0.0001) radius = 0.0001;
		ParametricPlane P = new ParametricPlane(bases[2].multiply((highestDot+lowestDot)/2).toPoint(), bases[0], bases[1]);

		//Slide half-cylinder caps along projection to the bases[0] and bases[1] planes.
		//dots contains the min and max along bases[0] and min and max along bases[1]
		double[] dots = {Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 
				Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY};
		for(Point p: points){
			double[] proj = P.projectPoint(p);
			double delta = Math.sqrt(radius*radius-proj[2]*proj[2]);
			if(radius*radius<proj[2]*proj[2]) delta = 0;
//			System.out.printf("Point: %s, projection: [%f , %f , %f], delta: %f .. %f %f\n", p,proj[0],proj[1],proj[2],delta, proj[0]+delta, proj[0]-delta);
			if(proj[0]+delta<dots[0]){	dots[0] = proj[0]+delta;	}
			if(proj[0]-delta>dots[1]){	dots[1] = proj[0]-delta;	}
			if(proj[1]+delta<dots[2]){	dots[2] = proj[1]+delta;	}
			if(proj[1]-delta>dots[3]){	dots[3] = proj[1]-delta;	}
		}
		//		System.out.printf("Final dots: %s\n",Arrays.toString(dots));
		//TODO: Fix corners


		double[] pars = { (dots[0]+dots[1])/2, (dots[2]+dots[3])/2 };
		double[] dim = { (dots[1]-dots[0])/2, (dots[3]-dots[2])/2 };
		Point center = new Point(P.getP(pars));
		Vector[] rssBases = {bases[0].multiply(dim[0]), bases[1].multiply(dim[1])};

		if(Double.isNaN(center.x())) throw new Error(" nana "+Arrays.toString(dots));
		
		return new RSS(center, rssBases, radius);
	}


	public static RSS createBoundingRSS_covariance(RSS s1, RSS s2) {
		PointList points = new PointList();
		for(Point p: s1.rectangle.getCorners()) points.add(p); 
		for(Point p: s2.rectangle.getCorners()) points.add(p);
		RSS ret = createBoundingRSS_covariance(points);
		ret.radius+=Math.max(s1.radius, s2.radius);
		return ret;
	}

	public RSS clone(){
		Point p = rectangle.center.clone();
		Vector[] bases = {rectangle.bases[0].clone(), rectangle.bases[1].clone()};
		return new RSS(p,bases,radius);
	}

}
