package ProGAL.geom3d.volumes;

import java.awt.Color;
import java.text.DecimalFormat;
import java.util.List;

import ProGAL.geom3d.Circle;
import ProGAL.geom3d.Line;
import ProGAL.geom3d.LineSegment;
import ProGAL.geom3d.Plane;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.Vector;
import ProGAL.geom3d.kineticDelaunay.Tri;
import ProGAL.geom3d.kineticDelaunay.Vertex;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.math.Constants;
import ProGAL.math.Matrix;
import ProGAL.math.Polynomial;

/**
 * Let the radius from the center of the hole to the center of the torus tube be R, and the radius of 
 * the tube be r. Then the equation in Cartesian coordinates for a torus azimuthally symmetric about 
 * the z-axis is 
 * 	              (R - sqrt(x^2+y^2))^2 + z^2 = r^2
 * and the parametric equations are
 *       x = (R + rcos(v)) cos(u)
 *       y = (R + rsin(v)) sin(u)
 *       z = rsin(v)
 * for u, v in [0, 2PI[
 * 
 * The three different classes of standard tori correspond to the three possible relative sizes of r 
 * and R. When R > r, the surface will be the familiar ring torus. The case R = r corresponds to the 
 * horn torus, which in effect is a torus with no "hole". The case R < r describes the self-intersecting 
 * spindle torus. When R = 0, the torus degenerates to the sphere.
 */

public class Torus  {
	protected Point center;
	protected Vector normal;
	protected double R;
	protected double r;


	public Torus(Point center, Vector normal, double majorRadius, double minorRadius) {
		this.center = center;
		this.normal = normal;
		this.R = majorRadius;
		this.r = minorRadius;
	}
	
	/** Returns TRUE if point q is in the interior of the torus */
	public boolean contains(Point q) {
		double RR = R*R;
		Vector vp = new Vector(center, q);
		return Math.pow(q.distanceSquared(center) + RR - r*r,2) - 4*RR*(vp.cross(normal).getLengthSquared()) < 0.0;
	}

	public Circle getMainCircle() {
		return new Circle(center, R, normal);
	}
	
	public Plane getMainPlane() {
		return new Plane(center, normal);
	}
	
	public Circle getPoloidalCircle() {
		return null;
	}
	
	public Circle getSweepingCircle() {
		Vector v = normal.getOrthonormal();
		return new Circle(center.add(v.scaleToLength(R)), r, normal.cross(v));
	}
	
	public Sphere getSweepingSphere() {
		Vector v = normal.getOrthonormal();
		return new Sphere(center.add(v.scaleToLength(R)), r); 
	}

	public double getSurfaceArea() {
		return 4*Math.PI*Math.PI*R*r;
	}
	public Circle getToroidalCircle() {
		return new Circle(center, R, normal);
	}
	
	public double getVolume() {
		return 2*Math.PI*Math.PI*R*r*r;
	}
	/** Two circles of major radius R, tilted to the center plane
		by slopes of r/R and -r/R and offset from the center by minor radius distance r */
	public Circle[] getVillarceauCircles () {
		return null;
	}
	
	public double getMajorRadius(){ return R; }
	
	public double getMinorRadius(){ return r; }
	
	public Vector getNormal(){ return normal; }
	
	public Point getCenter() {
		return center;
	}
	
	// Daisy
	/** Find up to four intersection points with circle C(p,n).
	    See http://ubm.opus.hbz-nrw.de/volltexte/2009/2037/pdf/diss.pdf **/
	public Point[] getIntersectionCircle(Circle C) {
//		DecimalFormat newFormat = new DecimalFormat("#.#########");
		/*Vector newCircleNormal = C.getNormal();
		Point newCircleCenter = C.getCenter();*/
		Vector e3 = new Vector(0,0,1);
		Vector rotAxis = normal.cross(e3);
		if (!(rotAxis.equals(new Vector(0,0,0)))) {
			rotAxis.normalizeThis();
		}
		double angle = normal.angle(e3);
		Matrix rotMatrix = Matrix.createRotationMatrix(angle, rotAxis);
		Vector translate = new Vector(-center.x(),-center.y(),-center.z());
//		Vector newTorusNormal = rotMatrix.multiply(normal).normalize();
		Vector newCircleNormal = rotMatrix.multiply(C.getNormal()).normalizeThis();
		Point CC0 = C.getCenter().add(translate);
		Point newCircleCenter = (Point)rotMatrix.multiplyIn(CC0);
		
		Point[] intersections = new Point[4];
		
		double ax, ay, az, bx, by, bz, nz;
		Vector a_, b_;
		
		double nx = newCircleNormal.get(0);
		double ny = newCircleNormal.get(1);
		if (nx*nx+ny*ny == 0) {
			ax = 1;
			ay = 0;
			az = 0;
			a_ = new Vector(ax, ay, az);
			
			bx = 0;
			by = 1;
			bz = 0;
			b_ = new Vector(bx, by, bz);
			
			nx = 0;
			ny = 0;
			nz = 1;
		} else {
			nz = newCircleNormal.get(2);
			a_ = newCircleNormal.getOrthogonal().normalizeThis();
			b_ = a_.cross(newCircleNormal).normalizeThis();
		}
		double[] m1 = {a_.get(0), b_.get(0), nx};
		double[] m2 = {a_.get(1), b_.get(1), ny};
		double[] m3 = {a_.get(2), b_.get(2), nz};
		double[][] mlist = {m1, m2, m3};
//		Matrix M = new Matrix(mlist); 
		//Point p = M.invert().multiply(center.subtract(pPoint));
//		Point p = new Point(0,0,0);
		double alph = newCircleCenter.dot(newCircleCenter);
		double gam = newCircleCenter.dot(a_);
		double del = newCircleCenter.dot(b_);
		double eps = alph + C.getRadius()*C.getRadius()+R*R-r*r;
		double a = 4*gam*gam - 4*R*R*(a_.get(0)*a_.get(0)+a_.get(1)*a_.get(1));
		double b = 4*del*del - 4*R*R*(b_.get(0)*b_.get(0)+b_.get(1)*b_.get(1));
		double f = 4*gam*del - 4*R*R*(a_.get(0)*b_.get(0)+a_.get(1)*b_.get(1));
		double l = 2*gam*eps - 4*R*R*(newCircleCenter.x()*a_.get(0) + newCircleCenter.y()*a_.get(1));
		double m = 2*del*eps - 4*R*R*(newCircleCenter.x()*b_.get(0) + newCircleCenter.y()*b_.get(1));
		double d =   eps*eps - 4*R*R*(newCircleCenter.x()*newCircleCenter.x() + newCircleCenter.y()*newCircleCenter.y());
		double[] er0 = {a, f, l};
		double[] er1 = {f, b, m};
		double[] er2 = {l, m, d};
		double[][] emat =  {er0, er1, er2};
		Matrix E = new Matrix(emat);
		
/*		Vector a1 = M.multiply(new Vector(1,0,0));
		double cPHI = a1.dot(M.multiply(new Vector(1,0,0)));
		double sPHI = -a1.dot(M.multiply(new Vector(0,1,0)));
		double[] ar0 = {cPHI, -sPHI, p.x()};
		double[] ar1 = {sPHI, cPHI, p.y()};
		double[] ar2 = {0, 0, 1};
		double[][] amat = {ar0, ar1, ar2};
		Matrix A = new Matrix(amat);
		System.out.println("A = "+A.toString());
		
		Matrix ENew = ((A.invert().transpose()).multiply(E)).multiply(A.invert());
		System.out.println("E' = "+ENew.toString());*/
		Matrix ENew = E;
		double[] pars = new double[5];
		pars[4] =    ENew.get(0,0)*C.getRadius()*C.getRadius() - 2*ENew.get(2,0)*C.getRadius() 				 +   ENew.get(2,2);
		//pars[4] =  Double.valueOf(newFormat.format(pars[4]));
		pars[3] = -4*ENew.get(1,0)*C.getRadius()*C.getRadius() + 4*ENew.get(1,2)*C.getRadius();
		//pars[3] =  Double.valueOf(newFormat.format(pars[3]));
		pars[2] = -2*ENew.get(0,0)*C.getRadius()*C.getRadius() + 4*ENew.get(1,1)*C.getRadius()*C.getRadius() + 2*ENew.get(2,2);
		//pars[2] =  Double.valueOf(newFormat.format(pars[2]));
		pars[1] =  4*ENew.get(0,1)*C.getRadius()*C.getRadius() + 4*ENew.get(1,2)*C.getRadius();
		//pars[1] =  Double.valueOf(newFormat.format(pars[1]));
		pars[0] =    ENew.get(0,0)*C.getRadius()*C.getRadius() + 2*ENew.get(0,2)*C.getRadius() 				 +   ENew.get(2,2);
		//pars[0] =  Double.valueOf(newFormat.format(pars[0]));
		for (int z = 0 ; z < 5 ; z++) {
			if (Math.abs(pars[z]) <= Constants.EPSILON) pars[z] = 0.0;
		}
		Polynomial func;
		int tmpDeg = 4;
        while ((tmpDeg != 0) && (pars[tmpDeg] == 0.0)) { tmpDeg--; }
        if (tmpDeg != 4) {
        	double[] coeffs = new double[tmpDeg+1];
        	for (int i = 0; i <= tmpDeg ; i++) {
        		coeffs[i] = pars[i];
        	}
        	func = new Polynomial(coeffs);
        } else {
        	func = new Polynomial(pars);
        }
//        System.out.println("func = "+func.toString());
//        System.out.println("func = "+pars[4]+"x^4 + "+pars[3]+"x^3 + "+pars[2]+"x^2 + "+pars[1]+"x + "+pars[0]);
		Double[] roots = Polynomial.solveQuartic(func.coeff);
		if (roots==null) return null;
/*		for (Double r : roots) {
			if (r==null) break;
			System.out.println("root = "+r);
		}*/
//		Point[] ps = new Point[4];
		
		if (roots.length == 0) return null;
		for (int k = 0 ; k<roots.length ; k++) {
//				if (roots0[k]==null) continue;
			double root = roots[k];
			double[] parameters = new double[4];
			parameters[0] = (newCircleCenter.x()-C.getRadius()*a_.get(0))*root*root + 2*C.getRadius()*b_.get(0)*root + newCircleCenter.x()+C.getRadius()*a_.get(0);
			parameters[1] = (newCircleCenter.y()-C.getRadius()*a_.get(1))*root*root + 2*C.getRadius()*b_.get(1)*root + newCircleCenter.y()+C.getRadius()*a_.get(1);
			parameters[2] = (newCircleCenter.z()-C.getRadius()*a_.get(2))*root*root + 2*C.getRadius()*b_.get(2)*root + newCircleCenter.z()+C.getRadius()*a_.get(2);
			parameters[3] = root*root+1;
//			ps[k] = new Point((parameters[0]/parameters[3]), (parameters[1]/parameters[3]), (parameters[2]/parameters[3]));
			Point pk = new Point((parameters[0]/parameters[3]), (parameters[1]/parameters[3]), (parameters[2]/parameters[3]));
			//ps[k].toScene(scene, 0.02, java.awt.Color.CYAN);
//			Point rotPoint = rotMatrix.invert().multiplyIn(ps[k]);
			Point rotPoint = (Point)rotMatrix.invert().multiplyIn(pk);
			intersections[k] = rotPoint.subtractThis(translate);
		}
		return intersections;
	}
	
	public void toScene(J3DScene scene, Color clr, int vStep, int uStep) {
		Point p = new Point(0,0,0);
		Vector n = new Vector(0,0,1);
		double vAngle = 2*Math.PI/vStep;
		double uAngle = 2*Math.PI/uStep;
		double rcos;
		for (int v = 0; v < vStep; v++) {
			rcos = r*Math.cos(v*vAngle);
			p.setZ(r*Math.sin(v*vAngle));
			Circle c = new Circle(p, R+rcos, n);
			c.toScene(scene, 0.002, 32, clr);
		}
		p.setX(R);
		p.setY(0);
		p.setZ(0);
		Vector m = new Vector(0,1,0);
		for (int u = 0; u < uStep; u++) {
			p.rotation(n, uAngle, new Point(0,0,0));//TODO: Origo tilf�jet. Check om korrekt
			m.rotation(n, uAngle);
			Circle c = new Circle(p,r,m);
			c.toScene(scene,  0.002, 32, clr);
		}
	}
	
	public void toScene(J3DScene scene, Color clr) {
		int iMax = 50;
		double delta = 2*Math.PI/iMax;
		
		Vector v = normal.getOrthonormal();
		Point p = center.add(v.scaleToLength(R));
		for (int i = 0; i < iMax; i++) {

			new Circle(p, r, v.cross(normal)).toScene(scene, 0.01, 16, clr);
//			new Sphere(p, r).toScene(scene, clr);
//			new Circle(p, r, v.cross(normal)).toScene(scene, 0.01, 32);
			p.rotation(normal, delta, center);
		}
				
	}
	
	public void toSceneSkeleton(J3DScene scene, Color clr, int vStep, int uStep) {
		Point p = new Point(0,0,0);
		double vAngle = 2*Math.PI/vStep;
		double uAngle = 2*Math.PI/uStep;
		double rcos;
		for (int v = 0; v < vStep; v++) {
			for (int u = 0; u < uStep; u++) {
				rcos = r*Math.cos(v*vAngle);
				p.setX((R + rcos)*Math.cos(u*uAngle));
				p.setY((R + rcos)*Math.sin(u*uAngle));
				p.setZ(r*Math.sin(v*vAngle));
				p.toScene(scene, 0.02, clr);
			}
		}
	}

	
	public static void main(String[] args) {
		
		System.out.print(new Vector(1,1,0).dot(new Vector(1,1,0)));
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		Point Tcenter = new Point(0, 0, 0);
		Vector Tnormal = new Vector(0,0,1).normalize();
//		Vector CNormal = new Vector(0, 1, 1).normalize();
		double TR = 0.5;
		double Tr = 0.3;
		Torus T = new Torus(Tcenter, Tnormal, TR, Tr);
//		Point CCenter = new Point(1, 0.2, 0);
//		double CRadius = 0.5;
//		Circle Circle = new Circle(CCenter, CRadius, CNormal);
		double alphaVal = 0.7;
		Point A = new Point(-0.07125192717748097, -0.35838495022033634, -0.21524016943473911);
		Point B = new Point(-0.5587823647379834, -0.4774072677661424, -0.8968986134203439);
//		Point A = new Point(1,0,0);
//		Point B = new Point(0.5, 0.5,0);
		Plane biP = Point.getBisector(A, B);
		Point torusCenter = biP.projectPoint(A);
		//Find major radius using pythagoras:
		double R = Math.pow(Math.abs(alphaVal*alphaVal-A.distanceSquared(torusCenter)), 0.5);
		
		//Torus
		Vector torusNormal = A.vectorTo(torusCenter).normalize();
		Circle Tcircle = new Circle(torusCenter, R, torusNormal);
		Torus torus = new Torus(torusCenter, torusNormal, R, alphaVal);
		//Circle1
		Circle circle = new Circle(new Point(0.00,0.00,-0.7455782353977968), 0.668358614580686, new Vector(0, 0, 1).normalize());
//		Circle circle = new Circle(new Point(0.50,-0.75,-0.01), 1.0, new Vector(0, 0, 1).normalize());
		
		
		Point C = new Point(-0.07125192717748097, -0.35838495022033634, -0.21524016943473911);
		Point D = new Point(-0.5587823647379834, -0.4774072677661424, -0.8968986134203439);
//		Point A = new Point(1,0,0);
//		Point B = new Point(0.5, 0.5,0);
		Plane biP2 = Point.getBisector(C, D);
		Point torusCenter2 = biP2.projectPoint(C);
		//Find major radius using pythagoras:
		double R2 = Math.pow(Math.abs(alphaVal*alphaVal-C.distanceSquared(torusCenter2)), 0.5);
		
		//Torus
		Vector torusNormal2 = C.vectorTo(torusCenter2).normalize();
		Circle Tcircle2 = new Circle(torusCenter2, R2, torusNormal2);
		Torus torus2 = new Torus(torusCenter2, torusNormal2, R2, alphaVal);
		//Circle1
		Circle circle2 = new Circle(new Point(0.00,0.00,-0.7455782353977968), 0.668358614580686, new Vector(0, 0, 1).normalize());
		

		Color transBlue = new Color(0, 0, 255, 2);
		circle.toScene(scene, 0.01, Color.blue);
		circle2.toScene(scene, 0.01, Color.BLACK);
		torus.toScene(scene, transBlue);
/*		Point center = new Point(0, 3, 1);
		Vector normal = new Vector(0, 1, 0);
		Vector circleNormal = new Vector(0, 0, 1).normalize();
		double R = 0.5;
		double r = 0.2;
		Torus torus = new Torus(center, normal, R, r);

	//	Point p = new Point(0.45, -0.30, 0.01);
	//	if (torus.contains(p)) System.out.println("inside"); else System.out.println("outside");
		torus.toScene(scene, Color.blue, 36, 72);
//TODO Put in test class
//		Point circleCenter = new Point(0, 3, 1.2);
//		double circleRadius = 0.5;
//		Circle C = new Circle(circleCenter, circleRadius, circleNormal);
//		Vector e3 = new Vector(0,0,1);
//		Vector rotAxis = normal.cross(e3);
//		double angle = normal.angle(e3);
//		Matrix rotMatrix = Matrix.createRotationMatrix(angle, rotAxis);
//		System.out.println("Rotation matrix = "+rotMatrix.toString());
//		Vector translate = new Vector(-center.x(),-center.y(),-center.z());
//		Vector newTorusNormal = rotMatrix.multiply(normal);
//		System.out.println("Torus normal = "+newTorusNormal.normalize().toString(6));
//		Vector newCircleNormal = rotMatrix.multiply(circleNormal);
//		Point newCircleCenter = rotMatrix.multiply(circleCenter.add(translate));
//		Torus transTorus = new Torus(center.add(translate), newTorusNormal.normalize(), R, r);
//		Circle transCircle = new Circle(newCircleCenter, circleRadius, newCircleNormal.normalize());
//		System.out.println("new torus : "+transTorus.center.toString()+" "+transTorus.normal.toString());
//		System.out.println("new circle :"+transCircle.toString());
//		
//		System.out.println("TransTorus center = "+transTorus.center);
//		System.out.println("T center          = "+T.center);
//		System.out.println("TransTorus normal = "+transTorus.normal);
//		System.out.println("T normal          = "+T.normal);
//		System.out.println("TransCircle center = "+transCircle.getCenter());
//		System.out.println("Circle center          = "+Circle.getCenter());
//		System.out.println("TransCircle normal = "+transCircle.getNormal());
//		System.out.println("Circle normal          = "+Circle.getNormal());*/
//		
//		Point[] intersections = torus.getIntersectionCircle(circle);
//		double error = 0;
//		if (intersections==null) {
//			System.out.println("No intersections!");
//		} else {
//			for (Point intersect : intersections) {
//				if (intersect == null) continue;
////				System.out.println("Intersectionns at : "+intersect);
//				intersect.toScene(scene, 0.02, java.awt.Color.ORANGE);
//				double radius = (new Tri(new Vertex(A), new Vertex(B), new Vertex(intersect))).getCircumRadius();
//				System.out.println("Distance from intersection point to TCircle = "+radius);
//				error += alphaVal-radius;
//			}
//		}
//		System.out.println("Average error = "+(Math.abs(error)/4));
//		intersections = torus2.getIntersectionCircle(circle2);
//		error = 0;
//		if (intersections==null) {
//			System.out.println("No intersections2!");
//		} else {
//			for (Point intersect : intersections) {
//				if (intersect == null) continue;
////				System.out.println("Intersectionns at : "+intersect);
//				intersect.toScene(scene, 0.02, java.awt.Color.GREEN);
//				double radius = (new Tri(new Vertex(C), new Vertex(D), new Vertex(intersect))).getCircumRadius();
//				System.out.println("Distance from intersection point to TCircle2 = "+radius);
//				error += alphaVal-radius;
//			}
//		}
//		System.out.println("Average error2 = "+(Math.abs(error)/4));
///*		double[] pars0 = new double[]{-200.99, 1};
//		Polynomial tmp = new Polynomial(pars0);
//		Polynomial poly0 = tmp.times(tmp);
//		System.out.println("1. poly = "+poly0.toString());
//		double[] roots0 = poly0.solveSecondDegree();
//		double[] pars1 = new double[]{0.980295, 1.9802, 1};
//		Polynomial poly1 = new Polynomial(pars1);
//		System.out.println("2. poly = "+poly1.toString());
//		double[] roots1 = poly1.solveSecondDegree();*/
///*		double[] roots = {-0.65465, 0.65465, -1.5275, 1.5275};
//		int j = 0;
//		Point[] ps = new Point[4];
//		Point[] psTrans = new Point[4];
//		Color[] colors = {java.awt.Color.GREEN, java.awt.Color.RED, java.awt.Color.YELLOW, java.awt.Color.GRAY};
//		for (int k = 0 ; k<roots.length ; k++) {
//			double root = roots[k];
//			double[] parameters = new double[4];
//			parameters[0] = (newCircleCenter.x()-circleRadius*0)*root*root+2*circleRadius*(1)*root+newCircleCenter.x()+circleRadius*0;
//			System.out.println("parameters0 = "+parameters[0]);
//			parameters[1] = (newCircleCenter.y()-circleRadius*-0)*root*root+2*circleRadius*0*root+newCircleCenter.y()+circleRadius*-0;
//			System.out.println("parameters1 = "+parameters[1]);
//			parameters[2] = (newCircleCenter.z()-circleRadius*-1)*root*root+2*circleRadius*(0)*root+newCircleCenter.z()+circleRadius*-1;
//			System.out.println("parameters2 = "+parameters[2]);
//			parameters[3] = root*root+1;
//			System.out.println("parameters3 = "+parameters[3]);
//			ps[j] = new Point((parameters[0]/parameters[3]), (parameters[1]/parameters[3]), (parameters[2]/parameters[3]));
//			psTrans[j] = new Point(rotMatrix.invert().multiply(ps[j]).subtract(translate));
//			scene.addShape(new Sphere(psTrans[j], 0.01), colors[j]);
//			
//			System.out.println("Intersection point = "+psTrans[j]);
//			j += 1;
//		}*/
///*		Point test1 = new Point(circleCenter.x()+(-0.41), circleCenter.y()+0.28618176042508375, 0.0);
//		Point test2 = new Point(circleCenter.x()+(-0.41), circleCenter.y()+(-0.28618176042508375), 0.0);
//		
//		Point test3 = new Point(circleCenter.x()+(-0.01), circleCenter.y()+(0.4998999899979995), 0.0);
//		Point test4 = new Point(circleCenter.x()+(-0.01), circleCenter.y()+(-0.4998999899979995), 0.0);
//*/	
///*		Point test1 = new Point(circleCenter.x()+(-0.48987), circleCenter.y()+0.10013, circleCenter.z());
//		Point test2 = new Point(circleCenter.x()+(-0.10013), circleCenter.y()+(0.48987), circleCenter.z());
//		Point test3 = new Point(circleCenter.x()+(-0.20311), circleCenter.y()+(-0.45688765347730725), 0.0);
//		Point test4 = new Point(circleCenter.x()+(-0.45689), circleCenter.y()+(-0.2031047215108501), 0.0);
//		Point[] orgTest = {test1, test2, test3, test4};
//		for (int i = 0 ; i<4 ; i++) {
//			System.out.println("Org. intersecction at : "+orgTest[i]);
//		}
//		test1.toScene(scene, 0.01, java.awt.Color.YELLOW);
//		test2.toScene(scene, 0.01, java.awt.Color.YELLOW);
//		test3.toScene(scene, 0.01, java.awt.Color.YELLOW);
//		test4.toScene(scene, 0.01, java.awt.Color.YELLOW);*/
///*		j = 0;
//		for (int k = 0 ; k<roots1.length ; k++) {
////			if (roots0[k]==null) continue;
//			double root = roots1[k];
//			double[] parameters = new double[4];
//			parameters[0] = (0.5-0.5*1)*root*root+2*0.5*0*root+0.5+0.5*1;
//			parameters[1] = (0.5-0.5*0)*root*root+2*0.5*0*root+0.5+0.5*0;
//			parameters[2] = (0-0.5*0)*root*root+2*0.5*0*root+0+0.5*0;
//			parameters[3] = root*root+1;
//			ps[j] = new Point(parameters[0]/parameters[3], parameters[1]/parameters[3], parameters[2]/parameters[3]);
////			scene.addShape(new Sphere(ps[j], 0.06), java.awt.Color.GREEN);
//			
//			System.out.println("Intersection point = "+ps[j]);
//			j += 1;
//		}*/
	}

}
