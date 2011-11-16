package ProGAL.geom3d.volumes;

import java.util.ArrayList;
import java.util.List;

import ProGAL.geom2d.Circle;
import ProGAL.geom3d.*;


/** 
 * An infinitely extended cylinder defined using a line and a radius. 
 */
public class InfCylinder {
	protected Line line;
	protected double radius;

	public InfCylinder(Line l, double r){
		line = l;
		radius = r;
	}

	public LSS capWithHalfSpheres(PointList enclosedPoints){
		double lowerT = Float.POSITIVE_INFINITY, upperT=Float.NEGATIVE_INFINITY;
		for (int i = 0; i < enclosedPoints.size(); i++) {
			Point p = enclosedPoints.get(i);
			Sphere s = new Sphere(p, radius+0.00001f);
			double[] intersections = s.intersectionParameters(line);

			if(intersections.length<2){
				double intersection = line.orthogonalProjectionParameter(p);
				if(intersection>upperT) upperT = intersection;
				if(intersection<lowerT) lowerT = intersection;
			}else{
				if(intersections[0]>upperT) upperT = intersections[0];
				if(intersections[1]<lowerT) lowerT = intersections[1];
			}
		}
		return new LSS(line.getPoint(lowerT), line.getPoint(upperT), radius);
	}

	public Cylinder capWithDiscs(PointList enclosedPoints){
		double lowerT = Float.POSITIVE_INFINITY, upperT=Float.NEGATIVE_INFINITY;
		//System.out.println("capWithDiscs(..) .. "+line);
		for(Point p: enclosedPoints){
			Plane plane = new Plane(p, line.getDir());
			double intersection = plane.getIntersectionParameter(line);
			if(intersection>upperT) upperT = intersection;
			if(intersection<lowerT) lowerT = intersection;
		}

		return new Cylinder(line.getPoint(lowerT), line.getPoint(upperT), radius);
	}


	/** Assumes dir is normalized */
	public static InfCylinder createMinRadCylinderFromDirection(PointList points, Vector dir){
		Plane p = new Plane(dir);
		List<ProGAL.geom2d.Point> points2d = new ArrayList<ProGAL.geom2d.Point>();
		Vector x = new Vector(dir.x(), dir.y(), dir.z()+1);
		if(dir.x()==0 && dir.x()==0) x.setX(x.x()+1);
		x = x.cross(dir).scaleToLength(1);

		Vector y = x.cross(dir).scaleToLength(1);
		for (int i = 0; i < points.size(); i++) {
			Point po = points.get(i);
			Point projected = p.projectPoint(po);
			ProGAL.geom2d.Point p2d = new ProGAL.geom2d.Point(x.dot(projected.toVector()), y.dot(projected.toVector()));
			points2d.add(p2d);
			System.out.println(p2d);
		}
		Circle mec = Circle.minimumEnclosingCircle_Welzl(points2d);//minimumEnclosingCircle_bruteforce(points2d);
		
		Point linePoint = x.multiply(mec.center().x()).addThis(y.multiply(mec.center().y())).toPoint();
		return new InfCylinder(new Line(linePoint, dir.clone()), mec.getRadius()+mec.getRadius()*0.001);
	}

	public String toString(){
		return String.format("InfCylinder[%s,%.2f]", line.toString(), radius);
	}

	/* 74HOps */
	public LSS capWithHalfSpheres(LSS c1, LSS c2) {
		double lowerT = Double.POSITIVE_INFINITY, upperT=Double.NEGATIVE_INFINITY;
		Point[] centers = {c1.segment.getA(), c1.segment.getB(), c2.segment.getA(), c2.segment.getB()};
		double[] rads = {c1.rad, c1.rad, c2.rad, c2.rad};
		for(int i=0;i<4;i++){
			Point p = centers[i];
			double rad = rads[i];
			
			Sphere s = new Sphere(p, radius-rad+0.00001);
			double[] intersections = s.intersectionParameters(line);//10HOps
			
			if(intersections.length<2){
				double intersection = line.orthogonalProjectionParameter(p);//7HOps
				if(intersection+rad>upperT) upperT = intersection+rad;
				if(intersection-rad<lowerT) lowerT = intersection-rad;
			}else{
				double upperIntersection = intersections[0];
				double lowerIntersection = intersections[1];
				if(upperIntersection>upperT) 	upperT = upperIntersection;
				if(lowerIntersection<lowerT) 	lowerT = lowerIntersection;
			}
		}//4*17=68HOps
		return new LSS(line.getPoint(lowerT), line.getPoint(upperT), radius);//6HOps
	}
}
