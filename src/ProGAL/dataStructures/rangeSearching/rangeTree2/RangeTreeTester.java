package ProGAL.dataStructures.rangeSearching.rangeTree2;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import ProGAL.Benchmark;
import ProGAL.geomNd.Point;

/**
 * Class used for testing the Range Tree. 
 * @author Søren Lynnerup 11.01.2012
 */
@SuppressWarnings("unused")
public class RangeTreeTester {

	/**
	 * Generates a Point at a given dimension and with coordinates in the given intervals.
	 */
	private static Point generatePoint(int dimensions, double[] mins, double[] maxs) {
		double[] coords = new double[dimensions];

		Random generator = new Random();

		// Generate a random number in the interval [mins_d, maxs_d) for each dimension.
		for(int d = 0; d < dimensions; d++) {
			double value = mins[d] + generator.nextDouble() * (maxs[d] - mins[d]);
			coords[d] = (double) value;
		}

		return new Point(coords);
	}

	/**
	 * Generates a list of points in the given ranges.
	 */
	private static List<Point> generatePoints(int dimensions, int number_of_points, double[] mins, double[] maxs) {
		List<Point> points = new ArrayList<Point>();

		for(int i = 0; i < number_of_points; i++) {
			points.add(generatePoint(dimensions, mins, maxs));
		}

		return points;
	}

	/**
	 * Generate a list points. We have an outside and inside interval: [  outside   [   inside   ]   outside   ]
	 * Half of points_outside points are generated on each side of the inside interval.
	 */
	private static List<Point> generateIntervals(int dimensions, int points_inside, int points_outside, double[] inside_mins, double[] inside_maxs, double[] outside_mins, double[] outside_maxs) {
		List<Point> points = new ArrayList<Point>();

		// Generate the left outside interval.
		points.addAll(generatePoints(dimensions, points_outside/2, outside_mins, inside_mins));

		// Generate the inside interval.
		points.addAll(generatePoints(dimensions, points_inside, inside_mins, inside_maxs));

		// Generate the right outside interval.
		points.addAll(generatePoints(dimensions, points_outside/2, inside_maxs, outside_maxs));

		// Randomly permutate the points.
		Collections.shuffle(points);

		return points;
	}

	/**
	 * Generate a list points in the given interval.
	 */
	private static List<Point> generateInterval(int dimensions, int numb_of_points, double[] mins, double[] maxs) {
		List<Point> points = new ArrayList<Point>();

		// Generate the points..
		points.addAll(generatePoints(dimensions, numb_of_points, mins, maxs));

		// Randomly permutate the points.
		Collections.shuffle(points);

		return points;
	}

	private static void testCorrectness(int tests, int dimensions, int points) {
		Random generator = new Random();
		boolean allcorrect = true;

		// Define the interval.
		double[] mins = new double[dimensions];
		double[] maxs = new double[dimensions];
		Arrays.fill(mins, 0);
		Arrays.fill(maxs, 100);

		// Define the query windows.
		double[] querymins = new double[dimensions];
		double[] querymaxs = new double[dimensions];

		// Print stuff.
		System.out.println("Testing correctness");
		System.out.println("-------------------");
		System.out.println("Tests: " + tests + "  Dimensions: " + dimensions + "  Points: " + points);

		for(int i = 0; i < tests; i++) {
			System.out.print("Test " + i + ": ");

			// Generate points randomly in the interval [0, 100]
			List<Point> list = generateInterval(dimensions, points, mins, maxs);

			// Generate the query window for each dimension.
			for(int d = 0; d < dimensions; d++) {
				// Generate a interval.
				double min = generator.nextInt(100);
				double max = generator.nextInt(100);

				if(min < max) {
					querymins[d] = min;
					querymaxs[d] = max;
				} else {
					querymins[d] = max;
					querymaxs[d] = min;
				}
			}

			// Build tree.
			RangeTree rangetree = new RangeTree(list, true);

			// Query tree.
			List<Point> query1 = rangetree.query(querymins, querymaxs);

			// Build tree.
			NaivRangeArray rangearray = new NaivRangeArray(list);

			// Query array.
			List<Point> query2 = rangearray.query(querymins, querymaxs);

			// Check if the queries have returned the same points.
			// The two queries are equal if the have the same length and the same elements.
			if(query1.containsAll(query2) && query1.size() == query2.size()) {
				System.out.print("OK! (" + query1.size() + " returned)\n");
				allcorrect = allcorrect && true;
			} else {
				System.out.print("FAIL!\n");
				allcorrect = false;
			}
		}

		// Print more stuff.
		System.out.println("-------------------");

		if(allcorrect) {
			System.out.println("All tests: OK!");
		} else {
			System.out.println("All tests: FAIL!");
		}
	}


	private static class QueryBenchmark extends Benchmark{
		private final int tests, dimensions, points; 
		private QueryBenchmark(int tests, int dimensions, int points) {
			super();
			this.tests = tests;
			this.dimensions = dimensions;
			this.points = points;
		}

		
		public void runBenchmark(){
			Random generator = new Random();

			// Print stuff.
			System.out.println("Benchmarking");
			System.out.println("-------------------");
			System.out.println("Tests: " + tests + "  Dimensions: " + dimensions + "  Points: " + points);
			System.out.print("Generating query windows: ");

			// Define the query windows.
			List<double[]> querymins = new ArrayList<double[]>();
			List<double[]> querymaxs = new ArrayList<double[]>();

			// Generate query windows.
			for(int i = 0; i < tests; i++) {
				double[] qmins = new double[dimensions];
				double[] qmaxs = new double[dimensions];

				// Generate the query window for each dimension.
				for(int d = 0; d < dimensions; d++) {
					// Generate a interval.
					double min = generator.nextInt(100);
					double max = generator.nextInt(100);

					// Specify a interval.
					//min = -2;
					//max = -1;

					if(min < max) {
						qmins[d] = min;
						qmaxs[d] = max;
					} else {
						qmins[d] = max;
						qmaxs[d] = min;
					}
				}

				querymins.add(qmins);
				querymaxs.add(qmaxs);
			}

			System.out.print("OK!\n");

			// Generate points.
			double[] pmins = new double[dimensions];
			double[] pmaxs = new double[dimensions];
			Arrays.fill(pmins, 0);
			Arrays.fill(pmaxs, 100);
			List<Point> pointlist = generateInterval(dimensions, points, pmins, pmaxs);

			// Build the Range Tree without Fractional Cascading.
			new RangeTree(pointlist, false);
			double t1 = getCPUTime();//System.nanoTime();
			RangeTree rangetree = new RangeTree(pointlist, false);
			double t2 = getCPUTime();//System.nanoTime();
			System.out.println("RTbuild: " + (t2 - t1) + "ms");

			// Perform a query "tests" number of times.
			for(int i = 0; i < tests; i++) {
				List<Point> query1 = rangetree.query(querymins.get(i), querymaxs.get(i));
			}
			t1 = getCPUTime();//System.nanoTime();
			for(int i = 0; i < tests; i++) {
				List<Point> query1 = rangetree.query(querymins.get(i), querymaxs.get(i));
			}
			t2 = getCPUTime();//System.nanoTime();
			System.out.println("RTquery: " + (t2 - t1) + "ms (" + ((t2 - t1) / tests) + "ms avrg)");

			// Build the Range Tree using Fractional Cascading.
			new RangeTree(pointlist, true);
			t1 = getCPUTime();//System.nanoTime();
			rangetree = new RangeTree(pointlist, true);
			t2 = getCPUTime();//System.nanoTime();
			System.out.println("RTFCbuild: " + (t2 - t1) + "ms");

			// Perform a query "tests" number of times.
			for(int i = 0; i < tests; i++) 
				rangetree.query(querymins.get(i), querymaxs.get(i));
			t1 = getCPUTime();//System.nanoTime();
			for(int i = 0; i < tests; i++) {
				List<Point> query1 = rangetree.query(querymins.get(i), querymaxs.get(i));
			}
			t2 = getCPUTime();//System.nanoTime();
			System.out.println("RTFCquery: " + (t2 - t1) + "ms (" + ((t2 - t1) / tests) + "ms avrg)");


			// Build the Range Tree using Fractional Cascading.
			new NaivRangeArray(pointlist);
			t1 = getCPUTime();//System.nanoTime();
			NaivRangeArray rangearray = new NaivRangeArray(pointlist);
			t2 = getCPUTime();//System.nanoTime();
			System.out.println("Arraybuild: " + (t2 - t1) + "ms");

			// Perform a query "tests" number of times.
			for(int i = 0; i < tests; i++) {
				List<Point> query1 = rangearray.query(querymins.get(i), querymaxs.get(i));
			}
			t1 = getCPUTime();//System.nanoTime();
			for(int i = 0; i < tests; i++) {
				List<Point> query1 = rangearray.query(querymins.get(i), querymaxs.get(i));
			}
			t2 = getCPUTime();//System.nanoTime();
			System.out.println("Arrayquery: " + (t2 - t1) + "ms (" + ((t2 - t1) / tests) + "ms avrg)");
		}
	}

	/**
	 * 
	 */
	public static void main(String[] args) {
		int points = 1000;
		int dimensions = 3;
		int tests = 10000;

		//testCorrectness(tests, dimensions, points);

		new QueryBenchmark(tests, dimensions, 1000).runBenchmark();
//		new QueryBenchmark(test, dimensions, 10000).runBenchmark();
		//		benchmarkQuery(tests, dimensions, 100000);
	}

}
