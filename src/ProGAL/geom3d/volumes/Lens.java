package ProGAL.geom3d.volumes;

import ProGAL.geom3d.Circle;
import ProGAL.geom3d.ParametricPlane;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.Vector;

/** 
 * A lens is the intersection between two spheres. 
 * 
 * @todo Add getters and setters for the focal points and radii
 * @author R.Fonseca
 */
public class Lens implements Volume {
	private final Sphere s0, s1;
	private Circle equator;
	private ParametricPlane plane;
	private double d0, d1, r;

	/**
	 * Construct a lens from the two spheres. A IllegalArgumentException is thrown if the spheres 
	 * have empty or single point intersection. Note that the two spheres given as arguments are 
	 * not stored, and subsequent changes to the lens should be performed through appropriate setters.  
	 */
	public Lens(Sphere s0, Sphere s1){
		this.s0 = s0.clone();
		this.s1 = s1.clone();
		double d = s0.getCenter().distance(s1.getCenter());
		d0 = (d*d-s1.getRadiusSquared()+s0.getRadiusSquared())/(2*d);
		d1 = d-d0;
		if(d0<=0 || d1<=0) throw new IllegalArgumentException("Lens spheres are not allowed to contain eachother");
		r = Math.sqrt(s0.getRadiusSquared()-d0*d0);
		if(Double.isNaN(r)) throw new IllegalArgumentException("Lens is undefined unless the lens-spheres intersect");
		Vector normal = s0.getCenter().vectorTo(s1.getCenter()).divideThis(d);
		Point center = s0.getCenter().add(normal.multiply(d0));
		equator = new Circle(center, r, normal);

		Vector x = new Vector(1,0.001,0.0001).crossThis(normal).normalizeThis();
		Vector y = normal.cross(x);
		plane = new ParametricPlane(center, x,y);
//		System.out.printf("d %.3f, d0 %.3f, d1 %.3f, r %.3f\n", d,d0,d1,r);
	}
	
	/** Radius of the equator of the lens */
	public double getRadius(){ 
		return r;
	}
	/** The distance from a focal point to the equator of the lens */
	public double getFocalDistance(int i){ 
		if(i==0) return d0; else return d1;
	}
	public double getSphereRadius(int i){
		if(i==0) return s0.getRadius(); else return s1.getRadius(); 
	}
	public Point getSphereCenter(int i){
		if(i==0) return s0.getCenter().clone(); else return s1.getCenter().clone();
	}
	
	@Override
	public Point getCenter() {
		double dSq = s0.getCenter().distanceSquared(s1.getCenter());
		double d0p = (dSq-s1.getRadius()*s1.getRadius()+s0.getRadius()*s0.getRadius())/(2*dSq);
		return s0.getCenter().add(s0.getCenter().vectorTo(s1.getCenter()).multiplyThis(d0p));
	}
	
	
	@Override
	public boolean overlaps(Volume vol) {
		throw new RuntimeException("Not yet implemented");
	}

	/** The volume of the lens */
	public double getVolume() {
		//http://mathworld.wolfram.com/Sphere-SphereIntersection.html
		double R = s0.getRadius();
		double r = s1.getRadius();
		double dSq = s0.getCenter().distanceSquared(s1.getCenter());
		double d = Math.sqrt(dSq);
		return Math.PI*(R+r-d)*(R+r-d)*(dSq+(2*d-3*r)*r+(2*d+6*r-3*R)*R)/12*d;
	}
	
	public Lens clone(){
		return new Lens(s0.clone(), s1.clone());
	}
	
	public double distance(Point p){
		Vector cp = equator.getCenter().vectorTo(p);
		
		if(cp.dot(equator.getNormal())>0){//Its on the s1 side
			Vector a0p = s0.center.vectorTo(p);double a0pDist = a0p.length();
			if(a0pDist<s0.radius) return 0;

			double a0Theta = Math.atan(r/d0);
			double a0Angle = a0p.angle(equator.getNormal());
			if(a0Angle<=a0Theta) return a0pDist-s0.radius;
			return discDistance(p);
			
		}else{//Its on the s0 side
			Vector a1p = s1.center.vectorTo(p);double a1pDist = a1p.length();
			double a1Theta = Math.atan(r/d1);
			double a1Angle = Math.PI-a1p.angle(equator.getNormal());
			if(a1Angle<=a1Theta) return a1pDist-s1.radius; 
			return discDistance(p);
		}
		
	}
	
	private double discDistance(Point p){
		//Follows: http://www.geometrictools.com/Documentation/DistanceCircle3Disk3.pdf
		double[] xyz = plane.projectPoint(p);
		double rSq = equator.getRadius()*equator.getRadius();
		double xySq = xyz[0]*xyz[0]+xyz[1]*xyz[1];
		if(xySq<=rSq) return xyz[2];
		return Math.sqrt(xySq + xyz[2]*xyz[2]+rSq-2*equator.getRadius()*Math.sqrt(xySq));
	}
//	private Point projectPointToDisc(Point p){
//		ProGAL.geom2d.Vector proj = new ProGAL.geom2d.Vector(plane.projectPoint(p));
//		double len = proj.length();
//		if(len>getRadius()) proj.multiplyThis(getRadius()/len);
//		
//		return plane.getP(proj.getCoords());
//	}

	private Point getCirclePoint(double s){
		return equator.getCenter().add(plane.v1.multiply(Math.cos(s)*getRadius()).addThis(plane.v2.multiply(Math.sin(s)*getRadius())));
	}
	
	public double distance(Lens l){
//		J3DScene scene = J3DScene.createJ3DSceneInFrame();
//		if(scene!=null){
//			scene.setAxisEnabled(true);
//			scene.addShape(s0, new Color(100,100,100,100));
//			scene.addShape(s1, new Color(100,100,100,100));
//			scene.addShape(l.s0, new Color(100,250,100,100));
//			scene.addShape(l.s1, new Color(100,250,100,100));
//			scene.addShape(new Sphere(s0.center,0.03), Color.BLACK);
//			scene.addShape(new Sphere(s1.center,0.03), Color.BLACK);
//			scene.addShape(new Sphere(l.s0.center,0.03), Color.BLACK);
//			scene.addShape(new Sphere(l.s1.center,0.03), Color.BLACK);
//		}
		Vector a0b0 = s0.center.vectorTo(l.s0.center).normalizeThis();
		Vector a0b1 = s0.center.vectorTo(l.s1.center).normalizeThis();
		Vector a1b0 = s1.center.vectorTo(l.s0.center).normalizeThis();
		Vector a1b1 = s1.center.vectorTo(l.s1.center).normalizeThis();

		//TODO: Replace angle with acos(..dot..)
		double a0b0Angle = 			a0b0.angle(equator.getNormal());
		double a0b1Angle = 			a0b1.angle(equator.getNormal());
		double a1b0Angle = Math.PI- a1b0.angle(equator.getNormal());
		double a1b1Angle = Math.PI-	a1b1.angle(equator.getNormal());
		double b0a0Angle = Math.PI- a0b0.angle(l.equator.getNormal());
		double b0a1Angle = Math.PI- a1b0.angle(l.equator.getNormal());
		double b1a0Angle = 			a0b1.angle(l.equator.getNormal());
		double b1a1Angle = 			a1b1.angle(l.equator.getNormal());

		double a0Theta = Math.atan(r/d0);
		double a1Theta = Math.atan(r/d1);
		double b0Theta = Math.atan(l.r/l.d0);
		double b1Theta = Math.atan(l.r/l.d1);
		double a0b0Dist = Math.max(0,s0.center.distance(l.s0.center)-s0.radius-l.s0.radius);
		double a0b1Dist = Math.max(0,s0.center.distance(l.s1.center)-s0.radius-l.s1.radius);
		double a1b0Dist = Math.max(0,s1.center.distance(l.s0.center)-s1.radius-l.s0.radius);
		double a1b1Dist = Math.max(0,s1.center.distance(l.s1.center)-s1.radius-l.s1.radius);
		if(a0b0Angle<=a0Theta && b0a0Angle<=b0Theta ) return a0b0Dist;
		if(a0b1Angle<=a0Theta && b1a0Angle<=b1Theta ) return a0b1Dist;
		if(a1b0Angle<=a1Theta && b0a1Angle<=b0Theta ) return a1b0Dist;
		if(a1b1Angle<=a1Theta && b1a1Angle<=b1Theta ) return a1b1Dist;

		//The shortest distance must lie between the discs or between a disc and a sphere. 
		
		Point disc1Point=null, disc2Point=null;
		double delta = Math.PI, scale = 0.5, deltaRed = 0.5;
		double sDist, tDist, tmpPlusDist=Double.POSITIVE_INFINITY, tmpMinusDist = 0;
		double s = 0, t=0;
		//		Point p = getPoint(s);
		//		pDist = c.discDistance(p);
		while(delta>0.001){
			Point tmpPlus = getCirclePoint(s+delta/100);
			tmpPlusDist = l.discDistance(tmpPlus);//c.discDistance(tmpPlus);
			Point tmpMinus = getCirclePoint(s-delta/100);
			tmpMinusDist = l.discDistance(tmpMinus);
			if(tmpPlusDist<tmpMinusDist){
				disc1Point = tmpPlus;
				s += delta*scale;
			}else{
				disc1Point = tmpMinus;
				s -= delta*scale;
			}
			delta*=deltaRed;
		}
		sDist = Math.min(tmpPlusDist,tmpMinusDist);
		if(sDist<0.0001) return 0;

		delta = Math.PI;
		while(delta>0.001){
			Point tmpPlus = l.getCirclePoint(t+delta/100);
			tmpPlusDist = discDistance(tmpPlus);
			Point tmpMinus = l.getCirclePoint(t-delta/100);
			tmpMinusDist = discDistance(tmpMinus);
			if(tmpPlusDist<tmpMinusDist){
				disc2Point = tmpPlus;
				t += delta*scale;
			}else{
				disc2Point = tmpMinus;
				t -= delta*scale;
			}
			delta*=deltaRed;
			
		}
		tDist = Math.min(tmpPlusDist,tmpMinusDist);
		if(tDist<0.0001) return 0;
//		if(scene!=null) scene.addShape(new LSS(disc1Point, disc2Point, 0.05), Color.BLUE, 5);

		//Now just to check if the distance could be futher shortened by moving along the sphere
		
		Vector pApB = disc1Point.vectorTo(disc2Point);
		double pApBAngle = pApB.angle(equator.getNormal());
		if(			pApBAngle<a0Theta) return s0.center.distance(disc2Point)-s0.radius;
		if(Math.PI-	pApBAngle<a1Theta) return s1.center.distance(disc2Point)-s1.radius;

		pApBAngle = pApB.angle(l.equator.getNormal());
		if(			pApBAngle<b0Theta) return l.s0.center.distance(disc1Point)-l.s0.radius;
		if(Math.PI- pApBAngle<b1Theta) return l.s1.center.distance(disc1Point)-l.s1.radius;
		
		return Math.min(sDist, tDist);
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Sphere a0 = new Sphere(new Point(0,1,0), Math.sqrt(2));
		Sphere a1 = new Sphere(new Point(2,1,0), Math.sqrt(2));
		Lens A = new Lens(a0,a1);

		Sphere b0,b1;
		Lens B;
		
		//Case 2: Distance should be from the equator of B to the sphere of A (a0). 
		b0 = new Sphere(new Point(2,4,0), Math.sqrt(2));
		b1 = new Sphere(new Point(6,4,0), Math.sqrt(10));
		B = new Lens(b0,b1);
		System.out.println(A.distance(B));
	}

}
