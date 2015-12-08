package ProGAL.geom3d.complex.mss;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.PriorityQueue;

import ProGAL.dataStructures.UnionFind;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.PointList;
import ProGAL.geom3d.complex.*;
import ProGAL.geom3d.complex.delaunayComplex.DelaunayComplex;
import ProGAL.geom3d.viewer.J3DScene;

public class MinimumSpanningSurface {
	public final DelaunayComplex delaunayComplex;
	public final List<CTriangle> triangles; 
	
	public MinimumSpanningSurface(DelaunayComplex dc){
		this.delaunayComplex = dc;
		this.triangles = new ArrayList<CTriangle>();
		buildSurface();
	}
	
	private void buildSurface(){
		UnionFind<CTetrahedron> ds = new UnionFind<CTetrahedron>();
		for(CTetrahedron tet: delaunayComplex.getAllTetrahedra()){
			if(!tet.containsBigPoint()) continue;
			CTetrahedron n;
			n = tet.getNeighbour(0); if(n!=null && n.containsBigPoint() && ds.find(tet)!=ds.find(n) ) ds.union(tet, n);
			n = tet.getNeighbour(1); if(n!=null && n.containsBigPoint() && ds.find(tet)!=ds.find(n) ) ds.union(tet, n);
			n = tet.getNeighbour(2); if(n!=null && n.containsBigPoint() && ds.find(tet)!=ds.find(n) ) ds.union(tet, n);
			n = tet.getNeighbour(3); if(n!=null && n.containsBigPoint() && ds.find(tet)!=ds.find(n) ) ds.union(tet, n);
		}
		
		
		PriorityQueue<CTriangle> backwardQueue = new PriorityQueue<CTriangle>(
				1000, new Comparator<CTriangle>(){
					public int compare(CTriangle t1, CTriangle t2) {
						return -Double.compare(t1.getArea(), t2.getArea());
					}
				} );
		
		backwardQueue.addAll(delaunayComplex.getTriangles());
		
		while(!backwardQueue.isEmpty()){
			CTriangle tri = backwardQueue.poll();
			CTetrahedron t1 = tri.getAdjacentTetrahedron(0);
			CTetrahedron t2 = tri.getAdjacentTetrahedron(1);
			if(ds.find(t1)==ds.find(t2))
				triangles.add(tri);
			else ds.union(t1, t2);
		}
		
	}
	
//	public static List<CTriangle> getMSTInducedTriangles(DelaunayComplex dc){
//		MinimumSpanningTree
//	}
	
	public static void main(String[] args) {
		List<Point> points = PointList.generatePointsInCube(100);
		DelaunayComplex dc = new DelaunayComplex(points);
		MinimumSpanningSurface mss = new MinimumSpanningSurface(dc);

		
		
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		for(Point p: points){
			p.toScene(scene, 0.02, Color.BLACK);
		}
		for(CTriangle tri: mss.triangles){
			scene.addShape(tri, new Color(100,100,200));
		}
	}

}
