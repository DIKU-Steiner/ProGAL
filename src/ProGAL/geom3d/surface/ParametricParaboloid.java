package ProGAL.geom3d.surface;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.Vector;
import ProGAL.math.Matrix;

//import static java.lang.Math.cos;
//import static java.lang.Math.sin;
 
public class ParametricParaboloid extends ParametricSurface{
	private double a, b, c;
	private Vector displacement;
	private Matrix rotation;

	public ParametricParaboloid(double a, double b, double c){
		this(a,b,c,new Vector(0,0,0));
	}

	public ParametricParaboloid(double a, double b, double c, Vector displacement){
		this.a = a;
		this.b = b;
		this.c = c;
		this.displacement = displacement;
		this.rotation = Matrix.createIdentityMatrix(3);
	}
	
	public Vector getDisplacement(){ return displacement.clone(); }
	public void setDisplacement(Vector v){
		displacement.set(v);
	}
	public double getA(){ return a; }
	public void setA(double v){ a = v;}
	public double getB(){ return b; }
	public void setB(double v){ b = v;}
	public double getC(){ return c; }
	public void setC(double v){ c = v;}
	public Matrix getRotation(){ return rotation.clone(); }
	public void setRotation(Matrix rot){ 
		for(int r=0;r<rotation.getM();r++){
			for(int c=0;c<rotation.getN();c++){
				rotation.set(r, c, rot.get(r, c));
			}
		}
		
	}
	

	public Point getPoint(double u, double v){
		return (Point)rotation.multiplyIn(new Point(a*(u+v),b*(v-u),c*v*u)).addThis(displacement);
//		return rotation.applyToIn(new Point(cos(c)*u-sin(c)*v,sin(c)*u+cos(c)*v,(v*v/(b*b)-u*u/(a*a)))).addThis(displacement);
	}
	public Point getPoint(ProGAL.geom2d.Point par){
		return getPoint(par.get(0), par.get(1));
	}

	public ParametricParaboloid clone(){
		ParametricParaboloid ret = new ParametricParaboloid(a, b, c, displacement.clone());
		ret.setRotation(rotation);
		return ret;
	}
	
	public String toString(){
		return String.format("ParametricParaboloid[a:%.3f,b:%.3f,c:%.3f,dis:%s,rot:\n%s",a,b,c,displacement.toString(3),rotation.toString(3));
	}
	
}
