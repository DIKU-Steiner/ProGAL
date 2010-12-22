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
	public final long getUserTime(){
		return bean!=null?bean.getCurrentThreadUserTime():System.nanoTime();
	}

	/** Returns the threads current CPU-time in ms. */
	public final long getCPUTime(){
		return bean!=null?bean.getCurrentThreadCpuTime():System.nanoTime();
	}
	
	public abstract void runBenchmark();
}
