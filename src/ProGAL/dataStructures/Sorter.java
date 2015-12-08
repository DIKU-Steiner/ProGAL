package ProGAL.dataStructures;

public abstract class Sorter {
	public void Sort(Set set, SortTool tool) { Sort(set, tool, false); }
	
	public abstract void Sort(Set set, SortTool tool, boolean descending);

}

