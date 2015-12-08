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

public class StrandSurface extends ParametricSurface{
	private Point[] cAlphas;
	private Point[] controls;
	private Vector[] tangents;
	private Vector[] normals;
	private double[] tks;
	private double bondWidth = 2.5;
	private double bondHeight = 0.3;

	StrandSurface(List<Point> cAlphas){
		this.cAlphas = new Point[cAlphas.size()];
		this.controls = new Point[cAlphas.size()-1];
		this.tks = new double[cAlphas.size()-1];
		this.tangents = new Vector[cAlphas.size()-1];
		this.normals = new Vector[cAlphas.size()-1];

		int c = 0;
		for(Point p: cAlphas) 	this.cAlphas[c++] = p;

		tks[0] = 0;
		for(c=1;c<cAlphas.size()-2;c++){
				controls[c] = Point.getMidpoint(cAlphas.get(c),cAlphas.get(c+1));
				if(c==1) 	tks[c] = 1.5;
				else		tks[c] = tks[c-1]+1;
		}
		tks[cAlphas.size()-2] = tks[cAlphas.size()-3]+1.5;
		
		controls[0] = cAlphas.get(0);
		controls[controls.length-1] = cAlphas.get(cAlphas.size()-1);
		
		tangents[0] = controls[0].vectorTo(controls[1]).multiplyThis(0.5/tks[1]);
		for(c=1;c<controls.length-1;c++){
			tangents[c] = controls[c-1].vectorTo(controls[c+1]).divideThis(tks[c+1]-tks[c-1]);
		}
		tangents[controls.length-1] = controls[controls.length-2].vectorTo(controls[controls.length-1]).multiplyThis(0.5/(tks[controls.length-1]-tks[controls.length-2]));
		
		for(c=1;c<controls.length-1;c++){
			normals[c] = cAlphas.get(c-1).vectorTo(cAlphas.get(c)).crossThis(cAlphas.get(c).vectorTo(cAlphas.get(c+1))).normalizeThis();
			if(c>1 && normals[c].dot(normals[c-1])<0) normals[c].multiplyThis(-1);
		}
		normals[0] = normals[1];
		normals[normals.length-1] = normals[normals.length-2];
	}


	@Override
	public Point getPoint(double u, double v) {
		int frame = Math.min((int)Math.floor(u), controls.length-2);
		double localT = u-frame;
		Vector p;
//		if(frame>=controls.length-2){
//			p = controls[frame-1].toVector().multiplyThis(h00(localT));
//			p.addThis(tangents[frame-1].multiply(h10(localT)));
//			p.addThis(controls[frame].toVector().multiply(h01(localT)));
//			p.addThis(tangents[frame].multiply(h11(localT)));
//		}else{
			p = controls[frame].toVector().multiplyThis(h00(localT));
			p.addThis(tangents[frame].multiply(h10(localT)));
			p.addThis(controls[frame+1].toVector().multiply(h01(localT)));
			p.addThis(tangents[frame+1].multiply(h11(localT)));
//		}

		Vector dir;
//		if(frame==controls.length-2){
//			dir = tangents[frame-1].multiply(h00(localT)).addThis(tangents[frame].multiply(h01(localT))).normalizeThis();
//		}else{
			dir = tangents[frame].multiply(h00(localT)).addThis(tangents[frame+1].multiply(h01(localT))).normalizeThis();
//		}
		Vector up;
//		if(frame==controls.length-2){
//			up = normals[frame-1].multiply(h00(localT)).addThis(normals[frame].multiply(h01(localT))).normalizeThis();
//		}else{
			up = normals[frame].multiply(h00(localT)).addThis(normals[frame+1].multiply(h01(localT))).normalizeThis();
//		}
		Vector out = up.cross(dir).normalizeThis();

		double arrowStartU = controls.length-1.5;
		if(u<arrowStartU){
			double norm = 2*bondWidth+2*bondHeight;
			if(v*norm<bondWidth/2) {
				p.addThis(out.multiplyThis(v*norm));
				p.addThis(up.multiplyThis(bondHeight/2));
			}else if(v*norm<bondWidth/2+bondHeight){
				p.addThis(out.multiplyThis(bondWidth/2));
				p.addThis(up.multiplyThis(bondHeight/2-(v*norm-bondWidth/2)));
			}else if(v*norm<3*bondWidth/2+bondHeight){
				p.addThis(out.multiplyThis(bondWidth/2-(v*norm-bondWidth/2-bondHeight)));
				p.addThis(up.multiplyThis(-bondHeight/2));
			}else if(v*norm<3*bondWidth/2+2*bondHeight){
				p.addThis(out.multiplyThis(-bondWidth/2));
				p.addThis(up.multiplyThis(-bondHeight/2+(v*norm-bondWidth/2-bondHeight-bondWidth)));
			}else{
				p.addThis(out.multiplyThis(-bondWidth/2+(v*norm-bondWidth/2-bondHeight-bondWidth-bondHeight)));
				p.addThis(up.multiplyThis(bondHeight/2));
			}
		}else{
			double outFact = 4*(controls.length-1-u);
			double norm = 2*bondWidth+2*bondHeight;
			if(v*norm<bondWidth/2) {
				p.addThis(out.multiplyThis(v*norm*outFact));
				p.addThis(up.multiplyThis(bondHeight/2));
			}else if(v*norm<bondWidth/2+bondHeight){
				p.addThis(out.multiplyThis(outFact*bondWidth/2));
				p.addThis(up.multiplyThis(bondHeight/2-(v*norm-bondWidth/2)));
			}else if(v*norm<3*bondWidth/2+bondHeight){
				p.addThis(out.multiplyThis(outFact*(bondWidth/2-(v*norm-bondWidth/2-bondHeight))));
				p.addThis(up.multiplyThis(-bondHeight/2));
			}else if(v*norm<3*bondWidth/2+2*bondHeight){
				p.addThis(out.multiplyThis(outFact*(-bondWidth)/2));
				p.addThis(up.multiplyThis(-bondHeight/2+(v*norm-bondWidth/2-bondHeight-bondWidth)));
			}else{
				p.addThis(out.multiplyThis(outFact*(-bondWidth/2+(v*norm-bondWidth/2-bondHeight-bondWidth-bondHeight))));
				p.addThis(up.multiplyThis(bondHeight/2));
			}
			p.addThis(dir.multiplyThis(-0.5*(controls.length-1-u)));
		}

		return p.toPoint();
	}


	public Vector getNormal(double u, double v){
		int frame = Math.min((int)Math.floor(u),controls.length-1);
		double localT = u-frame;

		Vector dir;
		if(frame==controls.length-1)
			dir = tangents[frame-1].multiply(h01(localT)).addThis(tangents[frame].multiply(h00(localT))).normalizeThis();
		else
			dir = tangents[frame].multiply(h00(localT)).addThis(tangents[frame+1].multiply(h01(localT))).normalizeThis();
		Vector up;
		if(frame==controls.length-1){
			up = normals[frame-1].multiply(h00(localT)).addThis(normals[frame].multiply(h01(localT))).normalizeThis();
		}else{
			up = normals[frame].multiply(h00(localT)).addThis(normals[frame+1].multiply(h01(localT))).normalizeThis();
		}
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
//			double z = i*(2*Math.PI/3.6);
//			Point p = new Point(1.5*Math.cos(z),1.5*Math.sin(z), z);
			Point p = new Point(i*2.3-5*2.4, i%2==0?1:0, 0.1*i*i);
			cAlphas.add(p);
		}
		for(Point p: cAlphas){
			scene.addShape(new Sphere(p,0.5));
		}
		StrandSurface surf = new StrandSurface(cAlphas);
		scene.addSurface(surf,new Color(200,200,20), 0,8,200, 0,1,20);

		for(double u=6;u<8;u+=0.2){
			Point p = surf.getPoint(u,0);
			System.out.println(p);
//			scene.addShape(new Sphere(p,(u-6)*0.2));
		}
//				scene.addShape(normalShape(surf,1.45,0.25));
//				scene.addShape(normalShape(surf,2,0));

	}
	private static LSS normalShape(StrandSurface surf, double u, double v){
		Point p = surf.getPoint(u,v);
		Vector n = surf.getNormal(u, v);
		return new LSS(p,p.add(n), 0.1);
	}
}
