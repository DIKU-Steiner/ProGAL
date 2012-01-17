package ProGAL.dataStructures;

import java.util.Random;

public class SorterQuick extends Sorter {
	private Set set;
	private SortTool tool;

	public void Sort(Set set, SortTool tool, boolean descending) {
		this.set = set;
		this.tool = tool;
		set.randomPermutation();
		int size = set.getSize();
		if (size > 1) partition(0, size-1);
	}
	
	private void partition(int left, int right) {
		int mid  = set.partition(tool, left, right);
		if (left < mid-1) partition(left, mid-1);
		if (mid+1 < right) partition(mid+1,right);
	}

	public static void main(String[] args) {
		Random randGen = new Random();
		Integer [] iArray = new Integer[10];
		for (int i = 0; i < 10; ++i) iArray[i] = new Integer(randGen.nextInt(99));
		for (int i = 0; i < 10; ++i) System.out.print(iArray[i] + " ");
		System.out.println();
		Set set = new Set(iArray);
		Sorter sort = new SorterQuick();
		sort.Sort(set,new SortToolInteger(), true);
		for (int i = 0; i < 10; ++i) System.out.print((Integer)set.get(i) + " ");
		System.out.println();
	}
}

