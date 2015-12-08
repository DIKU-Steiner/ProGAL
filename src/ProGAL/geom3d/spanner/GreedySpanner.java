package ProGAL.geom3d.spanner;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import ProGAL.dataStructures.Pair;
import ProGAL.dataStructures.WeightedGraph;
import ProGAL.dataStructures.shortestPath.Dijkstra;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.complex.CEdge;
import ProGAL.geom3d.complex.CTetrahedron;
import ProGAL.geom3d.complex.CVertex;
import ProGAL.geom3d.complex.delaunayComplex.DelaunayComplex;
import ProGAL.geom3d.complex.tessellation.DelaunayTessellation;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.LSS;
import ProGAL.geom3d.volumes.Tetrahedron;
import ProGAL.math.Randomization;

/**
 * The greedy spanner is the t-spanner that results from repeatedly adding an edge between the nearest 
 * pair of points. Since the number of pairs is quadratic, naive implementations run in O(n^3log(n)) time. 
 * A factory method is therefore supplied that computes the greedy spanner only on edges in the Delaunay 
 * tessellation. 
 * @author R.Fonseca
 *
 */
public class GreedySpanner extends WeightedGraph{
	private final Point[] points;
	private final Double[][] dMap; //Caching of distances
	
	public GreedySpanner(Point[] points, double t){
		this(points,t, createAllEdges(points.length));
	}
	
	private GreedySpanner(Point[] points, double t, List<Pair<Integer,Integer>> allEdges){
		super(points.length);
		this.points = points;
		this.dMap = new Double[points.length][points.length];
		super.weightFunction = new WeightFunction() {
			@Override
			public double w(Pair<Integer, Integer> p) {
				if(dMap[p.fst][p.snd]==null){
					dMap[p.fst][p.snd] = dMap[p.snd][p.fst] = GreedySpanner.this.points[p.fst].distance(GreedySpanner.this.points[p.snd]);  
				}
				return dMap[p.fst][p.snd];
			}
		};
		
		constructFromEdges(allEdges, t);
	}
	
	private static List<Pair<Integer,Integer>> createAllEdges(int vertices){
		List<Pair<Integer,Integer>> allEdges = new ArrayList<Pair<Integer,Integer>>(vertices*vertices);
		for(int i=0;i<vertices;i++){
			for(int j=i+1;j<vertices;j++){
				allEdges.add(new Pair<Integer,Integer>(i,j));
			}
		}
		return allEdges;
	}
	
	
	private void constructFromEdges(List<Pair<Integer,Integer>> allEdges, double t){
		Collections.sort(allEdges, new Comparator<Pair<Integer,Integer>>(){
			@Override
			public int compare(Pair<Integer, Integer> arg0,
					Pair<Integer, Integer> arg1) {
				return Double.compare(weightFunction.w(arg0), weightFunction.w(arg1));
			}});
		
		for(Pair<Integer,Integer> pair: allEdges){
			Dijkstra sp = new Dijkstra(this,pair.fst);
			double dG = sp.getDistance(pair.snd);
			double d = weightFunction.w(pair);
//			System.out.printf("Considering %s .. %.3f %.3f\n",pair, dG,d);
//			System.out.println(this);
			if(d*t<dG){
				addEdge(pair.fst, pair.snd);
//				System.out.println("> Adding edge");
			}
		}
	}

	public static GreedySpanner createSpannerFromDelaunay(Point[] points, double t){
		List<Point> pointList = new ArrayList<Point>(points.length);
		for(Point p: points) pointList.add(p);
		
		List<Pair<Integer,Integer>> edges = new ArrayList<Pair<Integer,Integer>>();
		DelaunayComplex dt = new DelaunayComplex(pointList);
		for(CEdge e: dt.getEdges()){
			int c0 = ((CVertex)e.getA()).idx;
			int c1 = ((CVertex)e.getB()).idx;
			edges.add(new Pair<Integer,Integer>(c0,c1));
			
		}
//		for(CTetrahedron tet: dt.getTetrahedra()){
//			int c0 = ((CVertex)tet.getCorner(0)).idx;
//			int c1 = ((CVertex)tet.getCorner(1)).idx;
//			int c2 = ((CVertex)tet.getCorner(2)).idx;
//			int c3 = ((CVertex)tet.getCorner(3)).idx;
//			edges.add(new Pair<Integer,Integer>(Math.min(c0, c1), Math.max(c0, c1)));
//			edges.add(new Pair<Integer,Integer>(Math.min(c0, c2), Math.max(c0, c2)));
//			edges.add(new Pair<Integer,Integer>(Math.min(c0, c3), Math.max(c0, c3)));
//			edges.add(new Pair<Integer,Integer>(Math.min(c1, c2), Math.max(c1, c2)));
//			edges.add(new Pair<Integer,Integer>(Math.min(c2, c3), Math.max(c2, c3)));
//			edges.add(new Pair<Integer,Integer>(Math.min(c3, c1), Math.max(c3, c1)));
//		}
		
		return new GreedySpanner(points, t, new ArrayList<Pair<Integer,Integer>>(edges));
		
	}
	
	public static void main(String[] args) {
		Randomization.seed(8);
		Point[] points = new Point[1000];
		for(int i=0;i<points.length;i++) points[i] = new Point(Point.getRandomPoint(3, 0.0, 1.0));
		GreedySpanner spanner = GreedySpanner.createSpannerFromDelaunay(points, 6);//new GreedySpanner(points, 4);
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
//		int i=0;
		for(Point p: points){
			p.toScene(scene, 0.0022, java.awt.Color.BLACK);
//			scene.addText((i++)+"", p);
		}
		for(Edge e: spanner.getEdges()){
			scene.addShape(new LSS(points[e.fst], points[e.snd], 0.002));
		}
		scene.centerCamera();
	}

}
