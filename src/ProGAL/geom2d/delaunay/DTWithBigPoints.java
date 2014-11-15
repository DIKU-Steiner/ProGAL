package ProGAL.geom2d.delaunay;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import ProGAL.geom2d.LineSegment;
import ProGAL.geom2d.Point;
import ProGAL.geom2d.viewer.J2DScene;
import ProGAL.geom3d.predicates.ExactJavaPredicates;
import ProGAL.math.Randomization;

/** 
 * 
 * @author ras
 */
public class DTWithBigPoints {
	private final ExactJavaPredicates pred = new ExactJavaPredicates();
	final List<Vertex> vertices;
	final List<Triangle> triangles;
	final Vertex[] bigPoints;

	public DTWithBigPoints(List<Point> points){
		this.vertices = new ArrayList<Vertex>(points.size());
		this.triangles = new ArrayList<Triangle>(points.size());
		bigPoints = new Vertex[3];
		bigPoints[0] = new Vertex(new Point(-3000,-3000)); 
		bigPoints[1] = new Vertex(new Point( 3000, 0));
		bigPoints[2] = new Vertex(new Point( 0, 3000));
		Triangle bigTri = new Triangle(bigPoints[0], bigPoints[1], bigPoints[2]);
		triangles.add(bigTri);
		for(Point p: points) 
			addPoint(p);
	}

	public List<Triangle> getTriangles(){ return new ArrayList<Triangle>(triangles); }
	public List<int[]> getEdges(){
		ArrayList<int[]> ret = new ArrayList<int[]>();
		for(Triangle tri: triangles){
			for(int i=0;i<3;i++){
				int[] e = new int[]{tri.corners[i].id-3, tri.corners[(i+1)%3].id-3};
				if(e[0]<0 || e[1]<0) continue;
				if(e[0]>e[1]){
					int tmp = e[0];
					e[0] = e[1];
					e[1] = tmp;
				}
				if( !ret.contains(e) )
					ret.add(e);
			}
			
		}
		return ret;
	}
	
	public Triangle walk(Point p){
		Triangle t = triangles.get(triangles.size()-1);
		while(true){
			double a1 = Point.area(t.corners[0], t.corners[1], p); orientPredCounter++;
			double a2 = Point.area(t.corners[1], t.corners[2], p); orientPredCounter++;
			if(a1<0 && a2<0) //No need to check the third side
				if(a2<a1) t = t.neighbors[0]; else t = t.neighbors[2];
			else{
				double a3 = Point.area(t.corners[2], t.corners[0], p); orientPredCounter++;
				if(Math.min(a1, Math.min(a2, a3))<0){ //Does t contain p yet?
					if(a3<a1 && a3<a2) 	t = t.neighbors[1];
					else if(a2<a1)		t = t.neighbors[0];
					else				t = t.neighbors[2];
				}else break;
			}
			
			//This is the straightforward way. The above is faster. 
			//if(Point.rightTurn(t.corners[0], t.corners[1], p)) t = t.neighbors[2];
			//else if(Point.rightTurn(t.corners[1], t.corners[2], p)) t = t.neighbors[0];
			//else if(Point.rightTurn(t.corners[2], t.corners[0], p)) t = t.neighbors[1];
			//else break;
		}
		return t;
	}
	public void addPoint(Point point){
		Vertex p = new Vertex(point);
		vertices.add((Vertex)p);
		Triangle t = walk(p);

		Triangle[] ts = splitTriangle(p, t);
		
		legalizeEdge(ts[0], 0);
		legalizeEdge(ts[1], 0);
		legalizeEdge(ts[2], 0);
	}

	private Triangle[] splitTriangle(Vertex p, Triangle t) {
		triangles.remove(t);

		Vertex[] ps = {t.corners[0], t.corners[1], t.corners[2]};
		Triangle[] ns = {t.neighbors[0], t.neighbors[1], t.neighbors[2]};
		Triangle[] ts = {new Triangle(p, ps[1], ps[2]), new Triangle(p, ps[2], ps[0]), new Triangle(p, ps[0], ps[1])};
		for(int i=0;i<3;i++){
			ts[i].neighbors[0] = ns[i]; 
			if(ns[i]!=null) ns[i].neighbors[ (ns[i].indexOf(ps[(i+1)%3])+1)%3 ] = ts[i];
			ts[i].neighbors[1] = ts[(i+1)%3]; 
			ts[i].neighbors[2] = ts[(i+2)%3];
			ps[i].first = ts[(i+1)%3];
			triangles.add(ts[i]);
		}
		p.first = ts[0];
		
		return ts;
	}

	

	void legalizeEdge(Triangle t, int e){
		if(t.neighbors[e]==null) return;

		Triangle u = t.neighbors[e];
		int f = (u.indexOf(t.corners[(e+1)%3])+1)%3;
		double inc = pred.incircle(t.corners[0].getCoords(),t.corners[1].getCoords(),t.corners[2].getCoords(), u.corners[f].getCoords());
		boolean illegal = inc>0;
		predCounter++;
		if(illegal){
			flip(t,e,u,f);
			//			scene.repaint();
			legalizeEdge(t, e);
			legalizeEdge(u, (f+2)%3);
		}

	}

	static void flip(Triangle t, int e, Triangle u, int f){
		flipCounter++;
		Vertex p0 = t.corners[e];
		Vertex p1 = t.corners[(e+1)%3];
		Vertex p2 = u.corners[f];
		Vertex p3 = t.corners[(e+2)%3];
		Triangle n01 = t.neighbors[(e+2)%3];
		Triangle n12 = u.neighbors[(f+1)%3];
		Triangle n23 = u.neighbors[(f+2)%3];
		Triangle n34 = t.neighbors[(e+1)%3];

		t.corners[(e+2)%3] = p2;
		u.corners[(f+2)%3] = p0;
		t.setCorner(p2, (e+2)%3);
		u.setCorner(p0, (f+2)%3);

		t.neighbors[e] = n12;
		t.neighbors[(e+1)%3] = u;
		t.neighbors[(e+2)%3] = n01;
		u.neighbors[f] = n34;
		u.neighbors[(f+1)%3] = t;
		u.neighbors[(f+2)%3] = n23;

		if(n12!=null) n12.neighbors[(n12.indexOf(p1)+1)%3] = t;
		if(n34!=null) n34.neighbors[(n34.indexOf(p3)+1)%3] = u;
		if(p1.first==u) p1.first = t;
		if(p3.first==t) p3.first = u;
	}


	public boolean inCH(Triangle t){
		if(t==null) return false;
		for(Point bp: bigPoints){
			for(int c=0;c<3;c++) if(t.corners[c]==bp) return false;
		}
		return true;
	}
	public boolean onCH(Triangle t){
		return inCH(t) && (!inCH(t.neighbors[0]) || !inCH(t.neighbors[1]) || !inCH(t.neighbors[2]));
	}

	public static int flipCounter = 0, predCounter = 0, orientPredCounter = 0;

	public static void main(String[] args){
		Randomization.seed(2);
		List<Point> points = new ArrayList<Point>();
		for(int i=0;i<1000;i++)
			points.add(new Point(Randomization.randBetween(-1.0, 1.0), Randomization.randBetween(-1.0, 1.0)));
		long start = System.nanoTime();
		DTWithBigPoints dt = new DTWithBigPoints(points);
		long end = System.nanoTime();
		System.out.printf("Took %.2fms. Used %d flips, %d predicates, %d orients\n", (end-start)/1000000.0, flipCounter, predCounter,  orientPredCounter);

		J2DScene scene = J2DScene.createJ2DSceneInFrame();
		scene.addShape(new LineSegment(new Point(0,0), new Point(0,1)), Color.black);
		scene.addShape(new LineSegment(new Point(0,0), new Point(1,0)), Color.black);
		//		for(Point p: points) {
		//			scene.addShape(new Circle(p,0.03), Color.BLACK);
		//			scene.addShape(new TextShape(points.indexOf(p)+"", p, 0.3), Color.BLACK);
		//		}
		for(Triangle t: dt.triangles){
			boolean fill = dt.onCH(t);
			scene.addShape(t, Color.GRAY, 0.0001, fill);
		}
		scene.autoZoom();
	}


}
