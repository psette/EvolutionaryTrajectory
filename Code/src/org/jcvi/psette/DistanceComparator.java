package org.jcvi.psette;


import java.util.Comparator;

public class DistanceComparator implements Comparator<ETstrain> {
	@Override
	public int compare(ETstrain a, ETstrain b) {
		return a.distancefrom < b.distancefrom ? -1
				: a.distancefrom == b.distancefrom ? 0 : 1;
	}

	static public void sort(ETstrain[] ads,
			Comparator<? super ETstrain> DistanceComparator) {
	}
}
