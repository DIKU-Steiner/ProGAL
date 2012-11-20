package ProGAL;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;

public abstract class Benchmark {
	private final ThreadMXBean bean;
	
	protected Benchmark(){
		if(!ManagementFactory.getThreadMXBean().isCurrentThreadCpuTimeSupported()){
			System.err.println("Warning: Thread timing not supported. Using System.nanoTime instead.");
			bean = null;
		}else{
			bean = ManagementFactory.getThreadMXBean();
		}
	}
	
	/** Returns the threads current user-time in ms. */
	public final double getUserTime(){
		return (bean!=null?bean.getCurrentThreadUserTime():System.nanoTime())/1000000.0;
	}

	/** Returns the threads current CPU-time in ms. */
	public final double getCPUTime(){
		return (bean!=null?bean.getCurrentThreadCpuTime():System.nanoTime())/1000000.0;
	}
	
	public abstract void runBenchmark();
}
