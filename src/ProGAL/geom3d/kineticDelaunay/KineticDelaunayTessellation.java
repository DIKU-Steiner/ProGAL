package ProGAL.geom3d.kineticDelaunay;

import java.awt.Color;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Stack;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.PointList;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.LSS;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.geomNd.Vector;
import ProGAL.math.Randomization;


public class KineticDelaunayTessellation {
	private List<Tet> tets = new LinkedList<Tet>();
	private Tet lastTet;
	
	public KineticDelaunayTessellation(List<Point> points){
		lastTet = new BigTet(points);
		
		
		for(Point p: points) insertPoint(p);		
	}
	

	public void insertPoint(Point p){
		Vertex v = new Vertex(p);
		System.out.println("Inserting "+v);
		Tet c = walk(p);
		
		//Corresponds to findNTes
		List<Tet> newTets = new LinkedList<Tet>();
		HashSet<Tet> processed = new HashSet<Tet>();
		processed.add(null);
		Stack<Tet> fringe = new Stack<Tet>();
		fringe.add(c);

		while(!fringe.isEmpty()){
			c = fringe.pop();
			if(processed.contains(c)) continue;
			for(int f=0;f<4;f++){
				Tet neigh = c.neighbors[f];
				if(neigh==null || !neigh.insideCircumsphere(p)){
					//-3 -2  0  1 .. c
					//-2 -1  0  1 .. neigh .. neigh.apex(c) should be 0
					//Create new cell
					Vertex[] corners = new Vertex[4];
					corners[3] = v;
					for(int i=1;i<4;i++) corners[i-1] = c.corners[(f+i)%4];
					Tet newTet = new Tet(corners);
					newTet.neighbors[3] = neigh;
					if(neigh!=null) neigh.neighbors[neigh.apex(c)] = newTet;
					newTets.add(newTet);
				}else if(!processed.contains(neigh)){
					fringe.add(neigh);
				}
			}
			processed.add(c);
			tets.remove(c);
		}
		
		lastTet = newTets.get(0); 
		restoreNeighborhood(newTets);
		tets.addAll(newTets);
	}

	private void restoreNeighborhood(List<Tet> newCells){
		for(Tet c1: newCells){
			cellLoop: for(Tet c2: newCells){
				if(c1==c2) break;
				
				//Check if c1 and c2 share a face based on corner vertices
				//Both will contain the last inserted vertex (highest index)
				//i and j run over the arrays, I and J record the location of differences
				int i=0,j=0, I=-1, J=-1;
				while(i<=3&&j<=3){
					if(c1.corners[i]==c2.corners[j]){
						i++;j++;
					}else if(i<3 && c1.corners[i+1]==c2.corners[j]){
						if(I>=0) continue cellLoop;
						I=i;i++;
					}else {
						if(J>=0) continue cellLoop;
						J=j;j++;
					}
				}
				c1.neighbors[I] = c2;
				c2.neighbors[J] = c1;
			}
		}
	}
	

	private Tet walk(Point p){
		Tet t = lastTet;
		mainWalk: while(true){
			for(int f=0;f<4;f++){
				if(!t.insideFace(f, p)){
					t = t.neighbors[f];
					lastTet = t;
					continue mainWalk;
				}
			}
			return t;
		}
	}
	
	public List<Tet> getTetrahedra(){ return tets; }
	
	
	public static void main(String[] args){
//		List<Point> points = new java.util.LinkedList<Point>();
//		points.add(new Point(0,0,0));
//		points.add(new Point(1,0,0));
//		points.add(new Point(0,1,0));
//		points.add(new Point(0,0,1));
//		points.add(new Point(1.1,1.1,1.1));

//		Vertex v0 = new Vertex(new Point(0,0,0));
//		Vertex v1 = new Vertex(new Point(1,0,0));
//		Vertex v2 = new Vertex(new Point(0,1,0));
//		Vertex v3 = new Vertex(new Point(0,0,1));
//		Tet t = new Tet(new Vertex[]{v0,v1,v2,v3});
//		System.out.println(t.inSphere(new Point(0,0,0)));
//		System.out.println();
//		System.out.println(t.inSphere(new Point(-0.1,0,0)));
//		System.out.println();
//		System.out.println(t.inSphere(new Point(0.1,0,0)));
		
		Randomization.seed(2);
		List<Point> points = PointList.generatePointsInCube(10);
		KineticDelaunayTessellation dt = new KineticDelaunayTessellation(points);
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		System.out.println(dt.getTetrahedra().size());
		for(Tet t: dt.getTetrahedra()){
			System.out.println(t);
			for(int i=0;i<4;i++){
				for(int j=i+1;j<4;j++){
					scene.addShape(new LSS(t.corners[i], t.corners[j], 0.01), Color.BLACK, 3);
				}
					
			}
		}
	}
	
}
