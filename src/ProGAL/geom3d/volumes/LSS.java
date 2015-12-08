package ProGAL.geom3d.volumes;

import ProGAL.geom3d.*;
import ProGAL.math.Constants;
import ProGAL.math.Matrix3x3;

/**
 * A line-segment swept sphere (also known as a line-swept-sphere, capsule or sometimes 
 * 'cigar') class. The LSS is represented by a line-segment and a radius, and is a cylinder 
 * capped with hemispheres. 
 * 
 * Distance calculations (and thereby collision checks) can be performed very fast, but 
 * finding the minimum capsule bounding a set of points can be somewhat time-consuming and 
 * no well-documented methods exist for doing this. For a heuristic see [Ericsson 05]. 
 */
public class LSS implements Volume{
	public LineSegment segment;
	public double rad;

	/**
	 * Construct a capsule using two endpoints (center of hemispheres) and a radius (used 
	 * both for hemispheres and cylinder shape). 
	 */
	public LSS(Point p1, Point p2, double r){
		this(new LineSegment(p1,p2), r);
	}

	/**
	 * Construct a capsule using a line-segment and a radius (used both for hemispheres and 
	 * cylinder shape). 
	 */
	public LSS(LineSegment segment, double r){
		this.segment = segment;
		this.rad = r;
	}

	public static LSS createBoundingLSS(PointList points){
		return createBoundingLSS_covariance(points);
	}

	public static LSS createBoundingLSS_covariance(PointList points){
		if(points.size()<=0)	throw new Error("Cannot create capsule enclosing 0 points");
		if(points.size()==1)	return new LSS(points.get(0).clone(), points.get(0).clone(), 0);
		if(points.size()==2)	return new LSS(points.get(0).clone(), points.get(1).clone(), 0);

		Matrix3x3 covMatr = points.getCovariance();
		covMatr.toConsole(3);
		Vector[] eigenVecs = covMatr.getEigenvectors();
//		eigenVecs[0].toConsole(3);
//		eigenVecs[1].toConsole(3);
//		eigenVecs[2].toConsole(3);
		if(eigenVecs[0]==null) eigenVecs[0] = eigenVecs[1].cross(eigenVecs[2]);
		if(eigenVecs[1]==null) eigenVecs[1] = eigenVecs[2].cross(eigenVecs[0]);
		if(eigenVecs[2]==null) eigenVecs[2] = eigenVecs[0].cross(eigenVecs[1]);
		
		Vector dir = eigenVecs[0];
		if ((eigenVecs[1] != null) && (eigenVecs[1].length()>dir.length())) dir = eigenVecs[1];
		if ((eigenVecs[2] != null) && (eigenVecs[2].length()>dir.length())) dir = eigenVecs[2];

		InfCylinder iCyl = InfCylinder.createMinRadCylinderFromDirection(points, dir.normalizeThis());
		LSS ret = iCyl.capWithHalfSpheres(points);
		return ret;
	}

	/* 284HOps */
	public static LSS createBoundingLSS_MaxDist(LSS v1, LSS v2) {
		double[] rads = {v1.rad, v1.rad, v2.rad, v2.rad};
		Point[] points = {v1.segment.getA(), v1.segment.getB(), v2.segment.getA(), v2.segment.getB()};

		int m1 = 0, m2 = 1;
		double best = v1.segment.getLength() + v1.rad+v1.rad;//4HOps
		double dist = v2.segment.getLength() + v2.rad+v2.rad;//4HOps
		if (dist > best) { best = dist;	m1 = 2;	m2 = 3; }
		
		double sumOfRads = v1.rad + v2.rad;
		for (int i=0; i<2; i++) {
			for (int j=2; j<4; j++) {
				dist = points[i].distance(points[j]) + sumOfRads;//4HOps
				if (dist > best) { best = dist;	m1 = i;	m2 = j; }
			}
		}

		Vector dir = points[m1].vectorTo(points[m2]).scaleToLength(1);//8HOps
		int exclude = 0;
		if(m1>1 && m2>1) exclude = 2;//A circle enclosing only 0,1 and 3 must be created
		else {
			if(rads[m1]>rads[m2]) exclude = m2;//A circle enclosing all but m2 must be created
			else exclude = m1;//A circle enclosing all but m1 must be created;
		}
		//32HOps so far.
		
		InfCylinder iCyl = createCylinderFromDirAndThreeSpheres(dir,rads,points,exclude);//126HOps
		LSS ret = iCyl.capWithHalfSpheres(v1, v2); //74HOps
		return ret;
	}
	
	/* 126 HOps */
	private static final InfCylinder createCylinderFromDirAndThreeSpheres(Vector dir, double[] rads, Point[] points, int exclude){
		Plane p = new Plane(new Point(0,0,0),dir);
		
		//Rand is guaranteed not to be parallel with dir
		Vector rand = new Vector(1,(dir.x()==-1||dir.x()==1)?1:0,0);
		
		Vector x = rand.cross(dir).scaleToLength(1);	//10HOps
		Vector y = x.cross(dir);						//6HOps
		ProGAL.geom2d.Circle[] cArr = new ProGAL.geom2d.Circle[3];
		int c=0;
		for(int i=0;i<3;i++){
			if(exclude==c) c++;
			Vector proj = p.projectPoint(points[c]).toVector();					//6HOps
			cArr[i] = new ProGAL.geom2d.Circle(new ProGAL.geom2d.Point(x.dot(proj), y.dot(proj)), rads[c]);	//6HOps
			c++;
		}
		//10+6+12*3=52HOps so far

		ProGAL.geom2d.Circle mec = new ProGAL.geom2d.Circle( cArr[0],cArr[1],cArr[2]);//68HOps
		Point linePoint = x.multiply(mec.center().x()).add(y.multiply(mec.center().y())).toPoint();//6HOps
		return new InfCylinder(new Line(linePoint,dir), mec.getRadius());
	}
	
	private static double clamp(double s){
		if(s<0) return 0;
		if(s>1) return 1;
		return s;
	}

	public double distanceToPoint(Point point){
		Vector d = segment.getAToB();
		double t = clamp( -(point.vectorTo(segment.getA()).dot(d))/(d.dot(d)) );
		return ( segment.getA().add(d.multiplyThis(t)).subtractThis(point.toVector()) ).toVector().length();
	}
	public boolean overlaps(LSS capsule){
		double minDist = closestSegmentPoint(capsule); 
		return minDist<=(rad+capsule.rad);
	}

	public double closestSegmentPoint(LSS capsule){
		Point startPoint1 = segment.getA();
		Point startPoint2 = capsule.segment.getA();

		Vector dir1 = segment.getAToB();
		Vector dir2 = capsule.segment.getAToB();
		double a = dir1.getLengthSquared();					//|S1| squared		.. 3HOp
		double e = dir2.getLengthSquared();					//|S2| squared		.. 3HOp
		
		if(a<Constants.EPSILON && e<Constants.EPSILON )
			return startPoint1.distance(startPoint2);
		if(a<Constants.EPSILON) return closestSegmentPoint(startPoint2, capsule.segment.getB(), startPoint1);
		if(e<Constants.EPSILON) return closestSegmentPoint(startPoint1, segment.getB(), startPoint2);
		
		Vector r = startPoint2.vectorTo(startPoint1);
		double f = dir2.dot(r);//               .. 3HOp
		double c = dir1.dot(r);//               .. 3HOp
		double b = dir1.dot(dir2);//			.. 3HOp
		double denom = a*e-b*b;//               .. 2HOp

		//If segments not parallel, compute closest point on L1 and L2
		//and clamp to S1
		double s, t;
		if(denom!=0.0f)  s = clamp( (b*f-c*e)/denom );//      .. 3HOp
		else             s = 0.0f;

		//Compute point on L2 closest to S1(S)
		double tnom = b*s+f;//                     .. 1HOp

		//If t in [0,1] done. Else clamp t and recompute and clamp s
		//.. 1 HOp
		if(tnom<0.0f){
			t = 0.0f;
			s = clamp(-c/a);
		}else if(tnom>e){
			t = 1.0f;
			s = clamp( (b-c)/a );
		}else{
			t = tnom/e;
		}

		Point c1 = startPoint1.add(dir1.multiplyThis(s));//      vec-scalar mult  .. 3HOp
		Point c2 = startPoint2.add(dir2.multiplyThis(t));//      vec-scalar mult  .. 3HOp
		return c1.distance(c2);	                //          .. 3HOp
	}
	
	private static double closestSegmentPoint(Point p11, Point p12, Point p2){
		Line l = new Line(p11, p11.vectorTo(p12));
		double t = l.orthogonalProjectionParameter(p2);
		t = clamp(t)*p11.distance(p12);
		return l.getPoint(t).distance(p2);
	}


	public boolean overlaps(Volume vol) {
		if(vol instanceof LSS) return overlaps((LSS)vol);
		throw new Error("Unimplemented");
	}

	public boolean contains(Point p){
		Line l = new Line(segment);
		double t = l.orthogonalProjectionParameter(p);
		if(t>1) t=1; else if(t<0) t=0;
		return l.getPoint(t).distance(p)<=this.rad;
	}


	public double getVolume() {
		double sphereVols = (4d/3d)*Math.PI*rad*rad*rad;
		double cylVol = Math.PI*rad*rad*segment.getLength();
		return sphereVols+cylVol;
	}

	
	public LSS clone(){
		return new LSS(segment.clone(), rad);
	}
	public Point getCenter() {
		return segment.getMidPoint();
	}

	public String toString(){
		return String.format("LSS[ls=%s,r=%f]", segment, rad);
	}

}
