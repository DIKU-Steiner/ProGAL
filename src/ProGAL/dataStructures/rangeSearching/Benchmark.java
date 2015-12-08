package ProGAL.dataStructures.rangeSearching;

import java.util.LinkedList;
import java.util.List;

import ProGAL.dataStructures.rangeSearching.rangeTree.RangeTree;
import ProGAL.geomNd.Point;
import ProGAL.geomNd.PointList;

public class Benchmark {

	private static void timeQuery(RangeSearchDataStructure ds, List<Point[]> queries){
		ds.query(queries.get(0)[0], queries.get(0)[1]);
		System.gc();
		long start = System.nanoTime();
		int check = 0;
		for(Point[] q: queries){
			List<Point> ret = ds.query(q[0], q[1]);
			check+=ret.size();
		}
		long end = System.nanoTime();
		System.out.printf("Done with %20s (check:%d) took %10.3fms. \n",ds.getClass().getSimpleName(), check, (end-start)/1000000.0);
	}
	
	private static List<Point[]> generateRandomQueries(int n, int d){
		List<Point[]> ret = new LinkedList<Point[]>();
		for(int i=0;i<n;i++){
			List<Point> highLow = PointList.generatePointsInCube(2, d);
			Point[] lowHigh = {highLow.get(0), highLow.get(1)};
			for(int c=0;c<d;c++){
				if(lowHigh[0].get(c)>lowHigh[1].get(c)){
					double tmp = lowHigh[0].get(c);
					lowHigh[0].set(c, lowHigh[1].get(c));
					lowHigh[1].set(c, tmp);
				}
			}
			ret.add(lowHigh);
		}
		return ret;
	} 
	private static void queryExperiment(int d, int n){
		
		List<Point> points = PointList.generatePointsInCube(n, d);
		GridMap gm = new GridMap(0.1, d);
		LinearRangeSearch lrs = new LinearRangeSearch(d);
		RangeTree rt = new RangeTree(points);
		for(Point p: points){
			gm.addPoint(p);
			lrs.addPoint(p);
		}
		List<Point[]> queries = generateRandomQueries(1000, d);

		timeQuery(gm, queries);
		timeQuery(lrs, queries);
		timeQuery(rt, queries);
		timeQuery(gm, queries);
		timeQuery(lrs, queries);
		timeQuery(rt, queries);
		timeQuery(gm, queries);
		timeQuery(lrs, queries);
		timeQuery(rt, queries);
		timeQuery(gm, queries);
		timeQuery(lrs, queries);
		timeQuery(rt, queries);
		
		
	}
	public static void main(String[] args) {
		queryExperiment(5, 10);
	}

}
