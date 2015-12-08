package ProGAL.proteins.viewer;

import java.awt.Color;
import java.util.LinkedList;
import java.util.List;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.Vector;
import ProGAL.geom3d.surface.ParametricSurface;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.LSS;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.math.Constants;

public class HelixSurface extends ParametricSurface{
	private Point[] cAlphas;
	private Vector[] tangents;
	private Vector[] normals;
	private Vector[] ups;
	private Vector axis;
	private double bondWidth = 0.2;
	private double bondHeight = 1.3;
	private double tension = -0.25;

	HelixSurface(List<Point> cAlphas){
		this.cAlphas = new Point[cAlphas.size()];
		this.tangents = new Vector[cAlphas.size()];
		this.normals = new Vector[cAlphas.size()];
		this.ups = new Vector[cAlphas.size()];
		this.axis = cAlphas.get(0).vectorTo(cAlphas.get(cAlphas.size()-1)).normalizeThis();

		int c = 0;
		for(Point p: cAlphas) 	this.cAlphas[c++] = p;

		for(c=0;c<cAlphas.size();c++){
			if(c==0 || c>=cAlphas.size()-1){
				if(c==0) {
					tangents[c] = cAlphas.get(c).vectorTo(cAlphas.get(c+1));
					normals[c] = tangents[c].cross(cAlphas.get(c+1).vectorTo(cAlphas.get(c+2))).crossThis(tangents[c]).normalizeThis();
				}else{
					tangents[c] = cAlphas.get(c-1).vectorTo(cAlphas.get(c));
					normals[c] = tangents[c].cross(cAlphas.get(c-2).vectorTo(cAlphas.get(c))).crossThis(tangents[c]).normalizeThis();
				}
			}else{
				tangents[c] = cAlphas.get(c-1).vectorTo(cAlphas.get(c+1)).multiplyThis(0.5*(1-tension));
				normals[c] = tangents[c-1].cross(tangents[c]).crossThis(tangents[c]).normalizeThis();
			}
		}
		ups[0] = axis;
		ups[1] = axis;
		ups[ups.length-1] = axis;
		ups[ups.length-2] = axis;
		for(c=2;c<cAlphas.size()-2;c++){
			ups[c] = this.cAlphas[c-2].vectorTo(this.cAlphas[c+2]).normalizeThis();
		}
		
	}


	@Override
	public Point getPoint(double u, double v) {
		int frame = (int)Math.floor(u);
		double localT = u-frame;
		Vector p;
		if(frame>=cAlphas.length-1){
			p = cAlphas[frame-1].toVector().multiplyThis(h01(localT));
			p.addThis(tangents[frame-1].multiply(h11(localT)));
			p.addThis(cAlphas[frame].toVector().multiply(h00(localT)));
			p.addThis(tangents[frame].multiply(h10(localT)));
		}else{
			p = cAlphas[frame].toVector().multiplyThis(h00(localT));
			p.addThis(tangents[frame].multiply(h10(localT)));
			p.addThis(cAlphas[frame+1].toVector().multiply(h01(localT)));
			p.addThis(tangents[frame+1].multiply(h11(localT)));
		}

		Vector dir;
		if(frame==cAlphas.length-1){
			dir = tangents[frame-1].multiply(h00(localT)).addThis(tangents[frame].multiply(h01(localT))).normalizeThis();
		}else{
			dir = tangents[frame].multiply(h00(localT)).addThis(tangents[frame+1].multiply(h01(localT))).normalizeThis();
		}
		Vector up;
		if(frame==cAlphas.length-1){
			up = ups[frame-1].multiply(h00(localT)).addThis(ups[frame].multiply(h01(localT))).normalizeThis();
		}else{
			up = ups[frame].multiply(h00(localT)).addThis(ups[frame+1].multiply(h01(localT))).normalizeThis();
		}
//		Vector up = axis;
//		if(frame>0 && frame<cAlphas.length-3){
//			up = cAlphas[frame-1].vectorTo(cAlphas[frame+3]).normalizeThis();
//		}
		Vector out = up.cross(dir).normalizeThis();

		p.addThis(out.multiplyThis(Math.cos(v*2*Math.PI)*bondWidth));
		if(frame==0 || frame==cAlphas.length-1){
			p.addThis(up.multiply(Math.sin(v*2*Math.PI)*(localT*bondHeight+(1-localT)*bondWidth)));
		}else if(frame==cAlphas.length-2){
			p.addThis(up.multiply(Math.sin(v*2*Math.PI)*((1-localT)*bondHeight+(localT)*bondWidth)));
		}else{
			p.addThis(up.multiply(Math.sin(v*2*Math.PI)*bondHeight));
		}

		return p.toPoint();
	}


	public Vector getNormal(double u, double v){
		int frame = (int)Math.floor(u);
		double localT = u-frame;

		Vector dir;
		if(frame==cAlphas.length-1)
			dir = tangents[frame-1].multiply(h01(localT)).addThis(tangents[frame].multiply(h00(localT))).normalizeThis();
		else
			dir = tangents[frame].multiply(h00(localT)).addThis(tangents[frame+1].multiply(h01(localT))).normalizeThis();
		Vector up = axis;
		Vector in = up.cross(dir).normalizeThis();

		return in.multiplyThis(Math.cos(v*Math.PI*2)).addThis(up.multiply(Math.sin(v*Math.PI*2)));
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
			double z = i*(2*Math.PI/3.6);
			Point p = new Point(1.5*Math.cos(z),1.5*Math.sin(z), z);
			cAlphas.add(p);
		}

		for(Point p: cAlphas){
//			scene.addShape(new Sphere(p,0.5));
		}
		HelixSurface surf = new HelixSurface(cAlphas);
		scene.addSurface(surf,new Color(200,20,20), 0,9,100, 0,1,10);

		//		scene.addShape(normalShape(surf,1.45,0.25));
		//		scene.addShape(normalShape(surf,2,0));

	}
	private static LSS normalShape(HelixSurface surf, double u, double v){
		Point p = surf.getPoint(u,v);
		Vector n = surf.getNormal(u, v);
		return new LSS(p,p.add(n), 0.1);
	}
}
