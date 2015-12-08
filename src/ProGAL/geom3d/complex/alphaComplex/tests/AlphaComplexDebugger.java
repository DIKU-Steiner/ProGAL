package ProGAL.geom3d.complex.alphaComplex.tests;

import java.awt.Color;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;

import ProGAL.geom3d.LineSegment;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.PointList;
import ProGAL.geom3d.Simplex;
import ProGAL.geom3d.Triangle;
import ProGAL.geom3d.complex.alphaComplex.AlphaFiltration;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.LSS;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.geom3d.volumes.Tetrahedron;
import ProGAL.math.Randomization;

public class AlphaComplexDebugger {
	public static void main(String[] args){
	
		debug1();
		
	}
	
	
	
	static J3DScene scene;
	static AlphaFiltration af;
	static int i=0;
	static int[][] betti;
	static void debug1(){
		Randomization.seed(0);
		PointList pl = PointList.generatePointsOnSphere(5);
//		PointList pl = new PointList();
//		pl.add(new Point(0,0,0));
//		pl.add(new Point(1,0,0));
//		pl.add(new Point(0,1,0));
//		pl.add(new Point(1,1,0));
//		pl.add(new Point(0,0,1));
//		pl.add(new Point(1,0,1));
//		pl.add(new Point(0,1,1));
//		pl.add(new Point(1,1,1));
		
//		pl.add(new Point(0,0,0));
//		pl.add(new Point(1,0,0));
//		pl.add(new Point(1,1,0));
//		pl.add(new Point(0,1,0));
//		pl.add(new Point(0.5,0.5,2));
//		pl.add(new Point(0.5,0.5,-2));
		af = new AlphaFiltration(pl);
		betti = af.getBettiNumbers();
		
		scene = J3DScene.createJ3DSceneInFrame();
		for(Point p: af.getVertices()){
			scene.addShape(new Sphere(p,0.03), new Color(0,0,0,100));
		}
		
//		for(Triangle t: af.getTriangles()){
//			scene.addShape(t.clone(), new Color(200,100,100,10));
//		}
		
//		Triangle tri23 = (Triangle)af.getSimplices().get(32);
//		System.out.println(af.getAttached(tri23));
//		System.out.println(af.getInAlpha(tri23));
//		System.out.println(  tri23.circumradius()  );
//		scene.addShape(new Sphere(tri23.circumcenter(), tri23.circumradius()), new Color(100,100,100,100));
//		scene.addShape(tri23, Color.blue);
		
		scene.getCanvas().addKeyListener(new KeyAdapter(){
			public void keyReleased(KeyEvent arg0) {
				if(arg0.getKeyChar()!='n') return;
				if(i==af.getSimplices().size()) return;
				
				Simplex s = af.getSimplices().get(i);
				System.out.printf("%d %d %d %d %d %d %.6f %s\n",i,betti[0][i],betti[1][i],betti[2][i],betti[3][i],betti[4][i], af.getInAlpha(s), s.toString());
				switch(s.getDimension()){
				case 0: scene.addShape(new Sphere((Point)s, 0.03), Color.BLACK);break;
				case 1: scene.addShape(new LSS((LineSegment)s, 0.005), Color.GRAY);break;
				case 2: scene.addShape((Triangle)s, new Color(100,100,250,200));break;
				case 3: scene.addShape((Tetrahedron)s, Color.BLACK);break;
				}
				
				i++;
			}
		});
	}
	
}
