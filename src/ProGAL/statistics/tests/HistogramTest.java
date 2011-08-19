package ProGAL.statistics.tests;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import ProGAL.statistics.Histogram;

public class HistogramTest {

	@Test
	public void testGetBin() {
		Histogram hist = new Histogram(5, new double[]{1,2,3,1});
		assertEquals(0, hist.getBin(-2));
		assertEquals(0, hist.getBin(5));
		assertEquals(0, hist.getBin(5.99));
		assertEquals(1, hist.getBin(6));
		assertEquals(1, hist.getBin(7.99));
		assertEquals(2, hist.getBin(8));
		assertEquals(2, hist.getBin(10.99));
		assertEquals(3, hist.getBin(11));
		assertEquals(3, hist.getBin(11.99));
		assertEquals(3, hist.getBin(20));
	}

	@Test
	public void testHistogramArray() {
		Histogram hist = new Histogram(new double[]{1,2,3,1});
		hist.addValue(0.1);//Added to bin 0
		hist.addValue(1.1);//Added to bin 1
		hist.addValue(2.1);//Added to bin 1
		hist.addValue(3.1);//Added to bin 2
		hist.addValue(4.1);//Added to bin 2
		assertArrayEquals(new int[]{1,2,2,0}, hist.histogramArray());
	}

	@Test
	public void testHistogramNormalizedArray() {

		Histogram hist = new Histogram(new double[]{1,2,3,1});
		hist.addValue(0.1);//Added to bin 0
		hist.addValue(1.1);//Added to bin 1
		hist.addValue(2.1);//Added to bin 1
		hist.addValue(3.1);//Added to bin 2
		hist.addValue(4.1);//Added to bin 2
		assertEquals(0.2*1, hist.histogramNormalizedArray()[0], 0.000001);
		assertEquals(0.2*2, hist.histogramNormalizedArray()[1], 0.000001);
		assertEquals(0.2*2, hist.histogramNormalizedArray()[2], 0.000001);
		assertEquals(0.2*0, hist.histogramNormalizedArray()[3], 0.000001);
	}

	@Test
	public void testBinCenters() {
		Histogram hist = new Histogram(5,new double[]{1,2,3,1});
		hist.addValue(5+0.1);//Added to bin 0
		hist.addValue(5+1.1);//Added to bin 1
		hist.addValue(5+2.1);//Added to bin 1
		hist.addValue(5+3.1);//Added to bin 2
		hist.addValue(5+4.1);//Added to bin 2
		assertEquals( 5.5, hist.binCenters()[0], 0.000001);
		assertEquals( 7.0, hist.binCenters()[1], 0.000001);
		assertEquals( 9.5, hist.binCenters()[2], 0.000001);
		assertEquals(11.5, hist.binCenters()[3], 0.000001);
	}

	@Test
	public void testCreateBalancedHistogram() {
		List<Double> values = new ArrayList<Double>();
		values.add(1.0);
		values.add(2.0);
		values.add(3.0);
		values.add(4.0);
		values.add(5.0);
		values.add(6.0);
		values.add(7.0);
		values.add(8.0);
		Histogram hist = Histogram.createBalancedHistogram(values, 4);
		assertEquals(2.0, hist.binWidths()[0], 0.001);
		assertEquals(2.0, hist.binWidths()[1], 0.001);
		assertEquals(2.0, hist.binWidths()[2], 0.001);
		assertEquals(1.0, hist.binWidths()[3], 0.001);

		hist = Histogram.createBalancedHistogram(values, 4);
		for(Double v: values) hist.addValue(v);
		int[] histArr = hist.histogramArray();
		for(int b1=0;b1<histArr.length;b1++){
			for(int b2=0;b2<b1;b2++){
				assertTrue(Math.abs(histArr[b1]-histArr[b2])<2);//At most 1 difference between occupancy
			}
		}
	}

}
