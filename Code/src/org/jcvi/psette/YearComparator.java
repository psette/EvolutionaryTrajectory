package org.jcvi.psette;

import java.util.Comparator;

public class YearComparator implements Comparator<ETstrain> {
	@Override
	public int compare(ETstrain a, ETstrain b) {
		return a.year < b.year ? -1
				: a.year == b.year ? 0 : 1;
	}
	public void sort(ETstrain[] sads,
            Comparator<? super ETstrain> YearComparator){
	}
}


