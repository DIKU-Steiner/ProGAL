package ProGAL.geomNd;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import ProGAL.geom2d.viewer.J2DScene;
import ProGAL.geom3d.viewer.J3DScene;

public class Simplex {
	private Point[] corners;
	
	public Simplex(Point[] corners){
		this.corners = corners;
	}
	
	public static Simplex regularSimplex(int dimension){
		int corners = dimension+1;
		double angle = -1.0/dimension;
		Point[] ret = new Point[corners];

		for(int s=0;s<corners;s++) ret[s] = new Point(new double[dimension]);

		double tmp = 0;
		for(int i=0;i<Math.min(corners,dimension);i++){
			double nCoord = Math.sqrt(1-tmp); 
			ret[i].set(i, nCoord);
			double rCoord = (angle - tmp)/nCoord;

			for(int s=i+1;s<corners;s++)
				ret[s].set(i, rCoord);

			tmp += rCoord*rCoord;
		}
		return new Simplex(ret);
	}
	
	public Simplex extend(int apex){
		Point[] newPointSet = new Point[corners.length];
		int c = 0;
		Vector faceAvg = new Vector(corners[0].getDimensions());
		for(int p=0; p<corners.length; p++){
			if(p==apex) continue;
			faceAvg.addThis(corners[p].toVector());
			newPointSet[c++] = corners[p];
		}
		faceAvg.divideThis(corners.length-1);
		
		Point newCorner = faceAvg.toPoint().addThis( corners[apex].vectorTo(faceAvg.toPoint()) );
		newPointSet[c++] = newCorner;
		return new Simplex(newPointSet);
	}
	
	public String toString(){
		return String.format("Simplex%s",Arrays.toString(corners));
	}

	public static void main(String[] args){
		cross_n16_d3(args);
	}
	
	public static void sausage_n16_d3(String[] args) {
		Set<Point> points = new HashSet<Point>();
		Simplex i = Simplex.regularSimplex(3);
		System.out.println(i); points.add(i.corners[0]);points.add(i.corners[1]);points.add(i.corners[2]);points.add(i.corners[3]);
		Simplex s = i;

		for(int it=0;it<13;it++){
			s = s.extend(0);
			System.out.println(s); points.add(s.corners[0]);points.add(s.corners[1]);points.add(s.corners[2]);points.add(i.corners[3]);
		}
		
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		int c=1;
		for(Point p: points){
			scene.addShape(new ProGAL.geom3d.volumes.Sphere(new ProGAL.geom3d.Point(p.getCoords()), 0.1), java.awt.Color.BLACK, 32);
			System.out.printf("DDD %d %.15f %.15f %.15f\n",c++, p.coords[0], p.coords[1], p.coords[2]);
		}
		
	}
	public static void cross_n16_d3(String[] args) {
		Set<Point> points = new HashSet<Point>();
		Simplex i = Simplex.regularSimplex(3);
		points.add(i.corners[0]);points.add(i.corners[1]);points.add(i.corners[2]);points.add(i.corners[3]);
		Simplex s = i;

		for(int it=0;it<10;it++){
			s = s.extend(0);
			points.add(s.corners[0]);points.add(s.corners[1]);points.add(s.corners[2]);points.add(i.corners[3]);
		}
		
//		s = i;
//		s = s.extend(1);
//		points.add(s.corners[0]);points.add(s.corners[1]);points.add(s.corners[2]);points.add(i.corners[3]);
//		s = s.extend(1);
//		points.add(s.corners[0]);points.add(s.corners[1]);points.add(s.corners[2]);points.add(i.corners[3]);
//	
//		for(int it=0;it<13;it++){
//			s = s.extend(0);
//			System.out.println(s); points.add(s.corners[0]);points.add(s.corners[1]);points.add(s.corners[2]);points.add(i.corners[3]);
//		}
	
		s = i;
		s = s.extend(2);
		points.add(s.corners[0]);points.add(s.corners[1]);points.add(s.corners[2]);points.add(i.corners[3]);
		s = s.extend(2);
		points.add(s.corners[0]);points.add(s.corners[1]);points.add(s.corners[2]);points.add(i.corners[3]);
	
		for(int it=0;it<10;it++){
			s = s.extend(0);
			System.out.println(s); points.add(s.corners[0]);points.add(s.corners[1]);points.add(s.corners[2]);points.add(i.corners[3]);
		}
		
		s = i;
		s = s.extend(3);
		points.add(s.corners[0]);points.add(s.corners[1]);points.add(s.corners[2]);points.add(i.corners[3]);
		s = s.extend(1);
		points.add(s.corners[0]);points.add(s.corners[1]);points.add(s.corners[2]);points.add(i.corners[3]);
	
		for(int it=0;it<10;it++){
			s = s.extend(0);
			System.out.println(s); points.add(s.corners[0]);points.add(s.corners[1]);points.add(s.corners[2]);points.add(i.corners[3]);
		}
		
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		int c=1;
		for(Point p: points){
			scene.addShape(new ProGAL.geom3d.volumes.Sphere(new ProGAL.geom3d.Point(p.getCoords()), 0.1), java.awt.Color.BLACK, 32);
			System.out.printf("DDD %d %.15f %.15f %.15f\n",c++, p.coords[0], p.coords[1], p.coords[2]);
		}
		
	}
	public static void bar_n16_d2(String[] args) {
		J2DScene scene = J2DScene.createJ2DSceneInFrame();
		Set<Point> points = new HashSet<Point>();
		Simplex i = Simplex.regularSimplex(2);
		System.out.println(i); points.add(i.corners[0]);points.add(i.corners[1]);points.add(i.corners[2]);
		Simplex s = i;

		for(int it=0;it<13;it++){
			s = s.extend(0);
			System.out.println(s); points.add(s.corners[0]);points.add(s.corners[1]);points.add(s.corners[2]);
		}
		
		int c=1;
		for(Point p: points){
			scene.addShape(new ProGAL.geom2d.Circle(new ProGAL.geom2d.Point(p.getCoords()), 0.1), java.awt.Color.BLACK, 0, true);
			System.out.printf("DD %d %.15f %.15f\n",c++, p.coords[0], p.coords[1]);
		}
	}

	public static void cross_n18_d2(String[] args) {
		J2DScene scene = J2DScene.createJ2DSceneInFrame();
		Set<Point> points = new HashSet<Point>();
		Simplex i = Simplex.regularSimplex(2);
		System.out.println(i); points.add(i.corners[0]);points.add(i.corners[1]);points.add(i.corners[2]);
		Simplex s = i;

		for(int it=0;it<5;it++){
			s = s.extend(0);
			System.out.println(s); points.add(s.corners[0]);points.add(s.corners[1]);points.add(s.corners[2]);
		}
		s = i.extend(2);
		System.out.println(s);points.add(s.corners[0]);points.add(s.corners[1]);points.add(s.corners[2]);
		for(int it=0;it<4;it++){
			s = s.extend(0);
			System.out.println(s);points.add(s.corners[0]);points.add(s.corners[1]);points.add(s.corners[2]);
		}

		s = i.extend(1);
		System.out.println(s);points.add(s.corners[0]);points.add(s.corners[1]);points.add(s.corners[2]);
		s = s.extend(1);
		System.out.println(s);points.add(s.corners[0]);points.add(s.corners[1]);points.add(s.corners[2]);
		for(int it=0;it<3;it++){
			s = s.extend(0);
			System.out.println(s);points.add(s.corners[0]);points.add(s.corners[1]);points.add(s.corners[2]);
		}
		
		int c=1;
		for(Point p: points){
			scene.addShape(new ProGAL.geom2d.Circle(new ProGAL.geom2d.Point(p.getCoords()), 0.1), java.awt.Color.BLACK, 0, true);
			System.out.printf("DD %d %.15f %.15f\n",c++, p.coords[0], p.coords[1]);
		}
	}

}
