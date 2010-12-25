package ProGAL.geom3d.predicates;

import ProGAL.geom3d.LineSegment;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.Triangle;
import ProGAL.geom3d.volumes.Tetrahedron;

public class ExactNativePredicates extends Predicates{

	private native void exactinit();
	private native double insphere(double x1, double y1, double z1, double x2, double y2, double z2,double x3, double y3, double z3, double x4, double y4, double z4, double x5, double y5, double z5);
	private native double orient3d(double x1, double y1, double z1, double x2, double y2, double z2,double x3, double y3, double z3, double x4, double y4, double z4);
	private native double tetcircumradius(double x1, double y1, double z1, double x2, double y2, double z2,double x3, double y3, double z3, double x4, double y4, double z4);
	private native double tricircumradius3d(double x1, double y1, double z1, double x2, double y2, double z2,double x3, double y3, double z3);
	private native double inspheretri(double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3, double q1, double q2, double q3);


	public ExactNativePredicates() {
		super();

		String os = System.getProperty("os.name");
		String bit = System.getProperty("sun.arch.data.model");
		String lib = null;

		if(os.contains("Linux") && bit.contains("32")){
			lib="ProGALPredicates_linux32";
			
		}else if(os.contains("Linux") && bit.contains("64")){
			lib="ProGALPredicates_linux64";
			
		}else if(os.contains("Mac")){
			lib="ProGALPredicates_mac";
			
		}else if(os.contains("Windows") && bit.contains("32")){
			lib = "libProGALPredicates_winxp32";
			
		}else{
			
			throw new Error(String.format("ProGAL has no primitive support for %s_%s yet",os,bit));
		}

		System.loadLibrary(lib);
		
		exactinit();		
	}
	
	@Override
	public double circumradius(Point p0, Point p1, Point p2, Point p3){
		double x0 = p0.getX(); double y0 = p0.getY(); double z0 = p0.getZ();
		double x1 = p1.getX(); double y1 = p1.getY(); double z1 = p1.getZ();
		double x2 = p2.getX(); double y2 = p2.getY(); double z2 = p2.getZ();
		double x3 = p3.getX(); double y3 = p3.getY(); double z3 = p3.getZ();

		return tetcircumradius(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);
	}

	@Override
	public double circumradius(Tetrahedron t){
		return circumradius(t.getPoint(0),t.getPoint(1),t.getPoint(2),t.getPoint(3));
	}
	
	@Override
	public double circumradius(Point p0, Point p1, Point p2){
		double radius = tricircumradius3d(	p0.getX(),p0.getY(),p0.getZ(),
				p1.getX(),p1.getY(),p1.getZ(),
				p2.getX(),p2.getY(),p2.getZ());

		return radius;
	}

	@Override
	public double circumradius(Triangle tri){
		double radius = tricircumradius3d(	tri.getPoint(0).getX(),tri.getPoint(0).getY(),tri.getPoint(0).getZ(),
				tri.getPoint(1).getX(),tri.getPoint(1).getY(),tri.getPoint(1).getZ(),
				tri.getPoint(2).getX(),tri.getPoint(2).getY(),tri.getPoint(2).getZ() );

		return radius;
	}


	@Override
	public double orient(Point p0, Point p1, Point p2, Point q){
		double orient = orient3d(	p0.getX(),p0.getY(),p0.getZ(),
				p1.getX(),p1.getY(),p1.getZ(),
				p2.getX(),p2.getY(),p2.getZ(),
				q.getX(),q.getY(),q.getZ());

		return orient;
	}

	@Override
	public SphereConfig insphere(Point p0, Point p1, Point p2, Point p3, Point q){
		double orient = orient(p0,p1,p2,p3);
		if(orient==0){
			return SphereConfig.COPLANAR;
		}

		double result = insphere(	p0.getX(),p0.getY(),p0.getZ(),
				p1.getX(),p1.getY(),p1.getZ(),
				p2.getX(),p2.getY(),p2.getZ(),
				p3.getX(),p3.getY(),p3.getZ(),
				q.getX(),q.getY(),q.getZ());

		if(orient>0){
			if(result>0){
				return SphereConfig.INSIDE;
			}
			if(result<0){
				return SphereConfig.OUTSIDE;
			}
			if(result==0){
				return SphereConfig.ON;
			}
		}
		if(orient<0){
			if(result>0){
				return SphereConfig.OUTSIDE;
			}
			if(result<0){
				return SphereConfig.INSIDE;
			}
			if(result==0){
				return SphereConfig.ON;
			} 
		}
		return null;
	}

	@Override
	public SphereConfig insphere(Tetrahedron t, Point q){
		return insphere(t.getPoint(0),t.getPoint(1),t.getPoint(2),t.getPoint(3),q);
	}

	@Override
	public SphereConfig insphere(Point p0, Point p1, Point p2, Point q){
		double result = inspheretri(	p0.getX(),p0.getY(),p0.getZ(),
				p1.getX(),p1.getY(),p1.getZ(),
				p2.getX(),p2.getY(),p2.getZ(),
				q.getX(),q.getY(),q.getZ()	);

		if(result>0) 		return SphereConfig.INSIDE;
		else if(result<0) 	return SphereConfig.OUTSIDE;		
		else				return SphereConfig.ON;
	}

	@Override
	public SphereConfig insphere(Triangle tri, Point q){
		return insphere(tri.getPoint(0),tri.getPoint(1),tri.getPoint(2),q);
	}

	@Override
	public PlaneConfig diffsides(Point p0, Point p1, Point p2, Point q0, Point q1){
		double a,b;

		a=orient(p0,p1,p2,q0);
		b=orient(p0,p1,p2,q1);

		if((a>0 && b<0) || (a<0 && b>0)) return PlaneConfig.DIFF;
		if((a>0 && b>0) || (a<0 && b<0)) return PlaneConfig.SAME;
		if(a==0 || b==0) return PlaneConfig.COPLANAR;

		//never happens
		return null;		
	}

	@Override
	public boolean inplane(Point p0, Point p1, Point p2, Point p3) {
		if(orient(p0,p1,p2,p3)==0) return true;
		else return false;
	}	

	/**TODO: This isnt exact */
	public SphereConfig edgeinsphere(LineSegment ls, Point q){
//		double x,y,z,r;		
//		double circumradius = edgecircumradius(e);
//		x= e.getPoint(0).getX()-q.getX();
//		y= e.getPoint(0).getY()-q.getY();
//		z= e.getPoint(0).getZ()-q.getZ();
//		r = Math.sqrt(Math.pow(x,2)+Math.pow(y,2)+Math.pow(z,2));
//		if(r<circumradius) return SphereConfig.INSIDE;
//		else if(r>circumradius) return SphereConfig.OUTSIDE;
//		else return SphereConfig.ON;
		double d_sq = ls.getMidPoint().getDistanceSquared(q);
		double r_sq = ls.getLengthSquared()/4;
		if(d_sq==r_sq) 	return SphereConfig.ON;
		if(d_sq<r_sq)	return SphereConfig.INSIDE;
		else			return SphereConfig.OUTSIDE;
	}

	/**TODO: This isnt exact */
	public double edgecircumradius(LineSegment ls){
//		double x,y,z,x0,x1,y0,y1,z0,z1;		
//		x0= e.getPoint(0).getX(); y0= e.getPoint(0).getY(); z0= e.getPoint(0).getZ();
//		x1= e.getPoint(1).getX(); y1= e.getPoint(1).getY(); z1= e.getPoint(1).getZ();
//		x= (x0+x1)/2; y= (y0+y1)/2; z= (z0+z1)/2;
//
//		return Math.sqrt( Math.pow(x-x0, 2) + Math.pow(y-y0, 2) + Math.pow(z-z0, 2)  );
		return ls.getLength()/2;
	}

}