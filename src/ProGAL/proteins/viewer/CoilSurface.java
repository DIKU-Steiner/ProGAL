package ProGAL.proteins.viewer;

import java.awt.Color;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.Vector;
import ProGAL.geom3d.surface.ParametricSurface;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.LSS;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.math.Constants;

public class CoilSurface extends ParametricSurface{
	private Point[] cAlphas;
	private Vector[] tangents;
	private double bondWidth = 0.3;
	private double tension = -0;

	CoilSurface(List<Point> cAlphaPoints){
		this.cAlphas = new Point[cAlphaPoints.size()];
		this.tangents = new Vector[cAlphaPoints.size()];
		int c = 0;
		for(Point p: cAlphaPoints) 	this.cAlphas[c++] = p;
		
		tangents[0] = cAlphas[0].vectorTo(cAlphas[1]).multiplyThis((1-tension)*0.5);
		for(c=1;c<cAlphas.length-1;c++){
			tangents[c] = cAlphas[c-1].vectorTo(cAlphas[c+1]).multiplyThis((1-tension)*0.5);
		}
		tangents[cAlphas.length-1] = cAlphas[cAlphas.length-2].vectorTo(cAlphas[cAlphas.length-1]).multiplyThis((1-tension)*0.5);
	}


	@Override
	public Point getPoint(double u, double v) {
		int frame = Math.min((int)Math.floor(u), cAlphas.length-2);
		double localT = u-frame;
		Vector p;
		p = cAlphas[frame].toVector().multiplyThis(h00(localT));
		p.addThis(tangents[frame].multiply(h10(localT)));
		p.addThis(cAlphas[frame+1].toVector().multiply(h01(localT)));
		p.addThis(tangents[frame+1].multiply(h11(localT)));

		Vector dir;
		dir = tangents[frame].multiply(h00(localT)).addThis(tangents[frame+1].multiply(h01(localT))).normalizeThis();
		Vector out = new Vector(-dir.y(),dir.x()+0.1, dir.z()).crossThis(dir).normalizeThis();
		Vector up = dir.cross(out);
		p.addThis(out.multiplyThis(bondWidth*Math.cos(v*2*Math.PI)));
		p.addThis(up.multiplyThis(bondWidth*Math.sin(v*2*Math.PI)));

		return p.toPoint();
	}


	public Vector getNormal(double u, double v){
		int frame = Math.min((int)Math.floor(u),cAlphas.length-2);
		double localT = u-frame;

		Vector dir;
		dir = tangents[frame].multiply(h00(localT)).addThis(tangents[frame+1].multiply(h01(localT))).normalizeThis();
		Vector out = new Vector(-dir.y(),dir.x()+0.1, dir.z()).crossThis(dir).normalizeThis();
		Vector up = dir.cross(out);

		return out.multiplyThis(Math.cos(v*2*Math.PI)).addThis(up.multiplyThis(Math.sin(v*2*Math.PI)));
	}

	private static double h00(double t){ return (1+2*t)*(1-t)*(1-t); }
	private static double h10(double t){ return t*(1-t)*(1-t); }
	private static double h01(double t){ return t*t*(3-2*t); }
	private static double h11(double t){ return t*t*(t-1); }

	public Point getCenter(){ return cAlphas[cAlphas.length/2]; }

	@Override
	public ParametricSurface clone() {
		return null;
	}


	/**
	 * @param args
	 */
	public static void main(String[] args) {
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		LinkedList<Point> cAlphas= new LinkedList<Point>();
		for(int i=0;i<10;i++){
//			double z = i*(2*Math.PI/3.6);
//			Point p = new Point(1.5*Math.cos(z),1.5*Math.sin(z), z);
			Point p = new Point(i*2.3-5*2.4, i%2==0?1:0, 0.1*i*i);
			cAlphas.add(p);
		}
		for(Point p: cAlphas){
			scene.addShape(new Sphere(p,0.5));
		}
		CoilSurface surf = new CoilSurface(cAlphas);
		scene.addSurface(surf,new Color(200,200,20), 0,9,200, 0,1,20);

//		for(double u=6;u<8;u+=0.2){
//			Point p = surf.getPoint(u,0);
//			System.out.println(p);
//			scene.addShape(new Sphere(p,(u-6)*0.2));
//		}
//				scene.addShape(normalShape(surf,1.45,0.25));
//				scene.addShape(normalShape(surf,2,0));

	}
	private static LSS normalShape(CoilSurface surf, double u, double v){
		Point p = surf.getPoint(u,v);
		Vector n = surf.getNormal(u, v);
		return new LSS(p,p.add(n), 0.1);
	}
}
