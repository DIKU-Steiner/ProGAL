package ProGAL.statistics;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

/** 
 * A class representing a one-dimensional histogram with irregular width bins. A min and max value 
 * can be specified as well as the individual bin widths. If a value higher than the max is added then it 
 * is added to the rightmost bin (vice versa for values below the min. value). 
 * @author R.Fonseca
 */
public class Histogram {
	private final List<Double> values = new LinkedList<Double>();
	private double minValue, binWidths[];
	private int bins;

	public Histogram(double minValue, double maxValue, int bins){
		this.minValue = minValue;
		this.bins = bins;
		this.binWidths = new double[bins];
		double width = (maxValue-minValue)/bins;
		for(int b=0;b<bins;b++) binWidths[b]=width;
	}
	
	public Histogram(double minValue, double[] binWidths){
		this.minValue = minValue;
		this.binWidths = binWidths;
		this.bins = binWidths.length;
	}
	
	public Histogram(double[] binWidths){		this(0,binWidths);		}

	public void addValue(double v){
		values.add(v);
	}
	public void addValues(List<Double> vals){
		values.addAll(vals);
	}
	
	public int getBin(double value){
		double bSup = minValue;
		for(int b=0;b<bins;b++){
			bSup+=binWidths[b];
			if(value<bSup) return b;
		}
		return bins-1;
	}

	public int[] histogramArray(){
		int[] ret = new int[bins];
		for(Double v: values){
			int bin = getBin(v);
			ret[bin]++;
		}
		return ret;
	}
	public double[] histogramNormalizedArray(){
		double inc = 1.0/values.size();
		double[] ret = new double[bins];
		for(Double v: values){
			int bin = getBin(v);
			ret[bin]+=inc;
		}
		return ret;
	}
	public double[] binCenters(){
		double[] ret = new double[bins];
		double bInf = minValue;
		for(int b=0;b<bins;b++){
			ret[b] = bInf+(binWidths[b]/2);
			bInf+=binWidths[b];
		}
		return ret;
	}
	public double[] binMinima(){
		double[] ret = new double[bins+1];
		double bInf = minValue;
		for(int b=0;b<bins;b++){
			ret[b] = bInf;
			bInf+=binWidths[b];
		}
		return ret;
	}
	public double[] binWidths(){
		double[] ret = new double[bins];
		for(int b=0;b<bins;b++){
			ret[b] = binWidths[b];
		}
		return ret;
	}
	
	public String toGnuplotString(int dec){
		StringBuilder sb = new StringBuilder();
		sb.append("#bin, bin-center, hist, hist-norm, bin-width\n");
		int[] hist = histogramArray();
		double[] histNorm = histogramNormalizedArray();
		double[] binCenters = binCenters();
		for(int b=0;b<bins;b++){
			sb.append(String.format("%d %."+dec+"f %d %."+dec+"f %."+dec+"f\n",b,binCenters[b],hist[b], histNorm[b], binWidths[b]));
		}
		return sb.toString();
	}
	public String toGnuplotString(){ return toGnuplotString(4); }
	
	
	public String toString(){
		StringBuilder sb = new StringBuilder();
		sb.append(String.format("IrregularHistogram[min:%.2f,widths:%s]\n",minValue, Arrays.toString(binWidths)));
		//TODO: Add an ascii normalized histogram
		sb.append(String.format("> %s\n",Arrays.toString(histogramArray())));
		return sb.toString();
	}
	
	
	/** 
	 * Construct an IrregularHistogram with min-value and bin-widths such that each bin contains roughly 
	 * the same number of values if the values are added. If all values are different then the occupancy 
	 * of one bin will, at most, be one different from the occupancy of any other bin. Note that the 
	 * returned histogram is empty.  
	 */
	public static Histogram createBalancedHistogram(List<Double> values, int bins){
		List<Double> sortedVals = new ArrayList<Double>(values);
		Collections.sort(sortedVals);
		double[] binWidths = new double[bins];
		
		double bInf = sortedVals.get(0);
		for(int b=0;b<bins;b++){
			int nextInfIdx = (sortedVals.size()*(b+1))/bins;
			double nextInf;
			if(nextInfIdx==sortedVals.size())	nextInf = sortedVals.get(sortedVals.size()-1)+0.001;
			else								nextInf = sortedVals.get(nextInfIdx);
			binWidths[b] = nextInf-bInf;
			bInf = nextInf;
		}
		
		return new Histogram(sortedVals.get(0), binWidths);
	}
	
}
