package ProGAL.dataStructures.rangeSearching.rangeTree;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import ProGAL.dataStructures.rangeSearching.RangeSearchDataStructure;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geomNd.Point;

public class RangeTree implements RangeSearchDataStructure {
	
	private final RangeTreeNode root;
	private final int dimensions;
	
	/**
	 * Initialize the range tree by presorting according to each dimension.
	 */
	public RangeTree(List<Point> points) {
		this.dimensions = points.get(0).getDimensions();
		
		List< List<Point> > sortedpoints = new ArrayList< List<Point> >();
		
		// Presort the points according to each dimension.
		for(int dimension = 0; dimension < dimensions; dimension++) {
			sortedpoints.add(dimension, new ArrayList<Point>(points));
			Collections.sort(sortedpoints.get(dimension), new PointDimensionComparator(dimension));
		}
		
		// Initialize the root of the 0th dimension.
		if(this.dimensions == 1) {
			this.root = new RangeTreeNode1d(sortedpoints);
		} else if(this.dimensions == 2) {
			this.root = new RangeTreeNode2d(sortedpoints, 0);
		} else {
			this.root = new RangeTreeNodeNd(sortedpoints, 0);
		}
	}
	
	/**
	 * Query the range tree given a low and a high value for each dimension.
	 * The values in high does not necessary has to be higher than the values in low.
	 * It is the interval which matters.
	 */
	public List<Point> query(double[] low, double[] high) {
		// Make sure low <= high for every dimension.
//		for(int dimension = 0; dimension < this.dimensions; dimension++) {
//			if(low[dimension] > high[dimension]) {
//				double temp = low[dimension];
//				low[dimension] = high[dimension];
//				high[dimension] = temp;
//			}
//		}
		
		return this.root.query(low, high);
	}
	
	/**
	 * Query the range given two points.
	 * As with the query above, the coordinates perform an interval.
	 * p1 < p2 is not a necessity.
	 */
	public List<Point> query(Point p1, Point p2) {
		return this.query(p1.getCoords(), p2.getCoords());
	}
	
	
	/**
	 * Test the range tree.
	 */
	public static void main(String[] args) {
		// Generate our point list.
		List<Point> list = new ArrayList<Point>();
		double[] coords;
		
		for(int i = 0; i < 9; i++) {
			coords = new double[3];
			coords[0] = i + ((i % 2) * 7);
			coords[1] = i + 3*(1 - 2*(i % 2));
			coords[2] = i + 5*(1 - 2*(i % 2));
			
			list.add(new Point(coords));
		}
		
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		for(Point p: list){
			ProGAL.geom3d.Point p3 = new ProGAL.geom3d.Point(p.getCoords());
			scene.addShape(new ProGAL.geom3d.volumes.Sphere(p3, 0.1));
		}
		
		
		/*
		coords = new double[3];
		coords[0] = 4;
		coords[1] = 5.0;
		coords[2] = 5.0;
		list.set(1, new Point(coords));
		coords = new double[3];
		coords[0] = 4;
		coords[1] = 5.0;
		coords[2] = 5.0;
		list.set(2, new Point(coords));
		coords = new double[3];
		coords[0] = 4;
		coords[1] = 5.0;
		coords[2] = 5.0;
		list.set(3, new Point(coords));
		coords = new double[3];
		coords[0] = 4;
		coords[1] = 5.0;
		coords[2] = 5.0;
		list.set(4, new Point(coords));
		coords = new double[3];
		coords[0] = 4;
		coords[1] = 5.0;
		coords[2] = 5.0;
		list.set(5, new Point(coords));
		*/

		RangeTree rangetree = new RangeTree(list);
		
		double[] c1 = new double[3];
		c1[0] = 8;
		c1[1] = -1;
		c1[2] = 1;
		double[] c2 = new double[3];
		c2[0] = 20;
		c2[1] = 11;
		c2[2] = 20;
		System.out.println(rangetree.query(c1, c2));
		
		ProGAL.geom3d.Point center = new ProGAL.geom3d.Point((c1[0]+c2[0])/2, (c1[1]+c2[1])/2, (c1[2]+c2[2])/2);
		ProGAL.geom3d.Vector[] bases = {
				new ProGAL.geom3d.Vector(1,0,0),
				new ProGAL.geom3d.Vector(0,1,0),
				new ProGAL.geom3d.Vector(0,0,1)
		};
		double[] extents = { (c2[0]-c1[0])/2, (c2[1]-c1[1])/2, (c2[2]-c1[2])/2 };
		ProGAL.geom3d.volumes.OBB query = new ProGAL.geom3d.volumes.OBB(center, bases, extents);
		scene.addShape(query, new java.awt.Color(200,50,50,100));
	}

}
