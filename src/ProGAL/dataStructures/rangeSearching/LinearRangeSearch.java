package ProGAL.dataStructures.rangeSearching;

import java.util.LinkedList;
import java.util.List;

import ProGAL.geomNd.Point;

public class LinearRangeSearch implements RangeSearchDataStructure{
	private final List<Point> points = new LinkedList<Point>(); 
	private final int dim;
	
	public LinearRangeSearch(int dimensions){
		this.dim = dimensions;
	}
	
	public void addPoint(Point p) {
		points.add(p);
	}

	@Override
	public List<Point> query(Point low, Point high) {
		List<Point> ret = new LinkedList<Point>();
		
		pointloop: for(Point p: points ){
			for(int d=0;d<dim;d++){
				if(p.get(d)<low.get(d) || p.get(d)>high.get(d))	continue pointloop;
			}
			ret.add(p);
		}
		return ret;
	}

	public static void main(String[] args) {
		LinearRangeSearch gm = new LinearRangeSearch(3);
		gm.addPoint(new Point(new double[]{0.7,0.7,0.7}));
		gm.addPoint(new Point(new double[]{0.3,0.3,0.3}));
		gm.addPoint(new Point(new double[]{-0.3,-0.3,-0.3}));
		gm.addPoint(new Point(new double[]{2.2,2.2,3.2}));
		List<Point> points = gm.query(new Point(new double[]{0.5,0.5,0.5}), new Point(new double[]{2.5,2.5,3.5}));
		System.out.println(points);
	}
}
