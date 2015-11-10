package org.jcvi.psette;

import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;
import java.util.Arrays;
import java.util.regex.Pattern;

public class TrajectoryCost {
	public static int ycord;
	public double distancefrom;
	public static double lowm, totalX, totalY, confl, confu;
	public String seq, strain_name, accesion_number, geographic_location, host,
			sequence;
	public static String[] headerparts;
	public int year;
	public static int n, max;
	private static int changes, yearpenalty;
	static double[] penalty, mutationrates;
	static String[][] hosts;
	static StringBuilder mostprobable = new StringBuilder();
	static StringBuilder mutations = new StringBuilder();

	static void intializer(ETstrain[] ads, macrodata[] macroarray,
			String[] header, String[] sequence, int[] startPos, int[] endPos,
			double smg, int bigcounter, int u, StringBuilder[] psa_results,
			StringBuilder[] high, StringBuilder[] medium)
			throws FileNotFoundException, UnsupportedEncodingException {
		mostprobable = new StringBuilder(" ");
		mostprobable
				.append("The following are the most probable common ancestors: ");
		penalty = new double[ads.length];
		mutationrates = new double[ads.length];
		totalY = 0;
		totalX = 0;
		max = 99999;
		for (int i = 0; i < ads.length; i++) {
			sequence[i] = sequence[i].substring(
					Math.max(startPos[0], startPos[i]),
					Math.min(endPos[0], endPos[i]));
			headerparts = header[i].split(Pattern.quote("|"));
			ads[i].accesion_number = headerparts[1];
			ads[i].year = Integer.parseInt(headerparts[2]
					.substring(headerparts[2].length() - 4));
			ads[i].strain_name = headerparts[0];
			ads[i].host = headerparts[4];
			ads[i].geographic_location = headerparts[3];
			ads[i].seq = sequence[i];
			ads[i].xcord = ads[0].year - ads[i].year;
			for (int l = 0; l < sequence[i].length(); l++) {
				if (ads[i].seq.charAt(l) != ads[0].seq.charAt(l)) {
					ads[i].ycord = ads[i].ycord + 1;
				}
			}

		}

		lowm = 99999;
		for (int i = 0; i < ads.length; i++) {
			float m = (float) (ads[i].ycord / ads[i].xcord);
			if (m < lowm && m > 0) {
				lowm = m;

			}
		}
		for (int i = 0; i < ads.length; i++) {
			ads[i].distancefrom = ads[i].ycord / ads[i].xcord - lowm;
		}
		comparing(ads, macroarray, bigcounter, u, psa_results, high, medium,
				smg);
	}

	public static void comparing(ETstrain[] ads, macrodata[] macroarray,
			int bigcounter, int u, StringBuilder[] psa_results,
			StringBuilder[] high, StringBuilder[] medium, double smg)
			throws FileNotFoundException, UnsupportedEncodingException {
		DistanceComparator dc = new DistanceComparator();
		Arrays.sort(ads, dc);
		YearComparator yc = new YearComparator();
		double confu = 0;
		double confl = 99999;
		for (int z = 0; z < ads.length; z++) {
			changes = 0;
			yearpenalty = 0;
			ETstrain[] sads = new ETstrain[z + 1];
			for (int q = 0; q <= z; q++) {
				sads[q] = new ETstrain();
				sads[q].year = ads[q].year;
				sads[q].seq = ads[q].seq;
				sads[q].xcord = ads[q].xcord;
				sads[q].ycord = ads[q].ycord;
				if (sads[q].seq.length() < max) {
					max = sads[q].seq.length();
				}

			}

			Arrays.sort(sads, yc);
			changes = 0;
			for (int a = 0; a < max - 1; a++) {
				for (int i = 0; i < z - 2; i++) {

					if (sads[i].seq.charAt(a) != sads[i + 1].seq.charAt(a)
							&& sads[i].seq.charAt(a) != '-'
							&& sads[i + 1].seq.charAt(a) != '-') {
						changes = changes + 1;
						yearpenalty = Math.abs(sads[i].year - sads[i + 1].year);

					}
				}
			}

			penalty[z] = changes / (z + 1);
			penalty[z] = penalty[z] / Math.pow((z + 1), .5);

			changes = 0;
		}
		Arrays.sort(ads, yc);
		double lowp = 999999999;

		for (int z = 5; z < ads.length - 1; z++) {
			if (penalty[z] < lowp) {
				lowp = penalty[z];
				n = z;
			}

		}
		for (int i = 0; i < ads.length; i++) {
			if ((ads[i].ycord / ads[i].xcord <= (ads[n].ycord / ads[n].xcord))) {
				totalY = totalY + ads[i].ycord;
				totalX = totalX + ads[i].xcord;
			}
		}
		System.out.println(totalY / totalX);
		for (int z = 0; z < n; z++) {
			if ((ads[z].ycord / ads[z].xcord) > confu) {
				confu = ads[z].ycord / ads[z].xcord;
			}
			if ((ads[z].ycord / ads[z].xcord) < confl) {
				confl = ads[z].ycord / ads[z].xcord;
			}
		}

		for (int a = 0; a < n; a++) {
			mostprobable.append(ads[a].accesion_number + " ");
		}

		Interp(n, ads, macroarray, bigcounter, u, psa_results, high, medium,
				smg);
	}

	public static void Interp(int n, ETstrain[] ads, macrodata[] macroarray,
			int bigcounter, int u, StringBuilder[] psa_results,
			StringBuilder[] high, StringBuilder[] medium, double smg)
			throws FileNotFoundException, UnsupportedEncodingException {
		StringBuilder[] analyze = new StringBuilder[n];
		ETstrain[] winners = new ETstrain[n];
		int[] count = new int[n];
		int a = 0;
		int b = 0;
		int secondarycount = 0;
		int highcount = 0;
		for (int i = 0; i < n; i++) {
			count[i] = 0;
			winners[i] = new ETstrain();
			winners[i] = ads[i];

		}
		for (int i = 0; i <= n - 1; i++) {
			analyze[i] = new StringBuilder();
			if (analyze[0].indexOf(winners[i].host) == -1) {
				secondarycount = secondarycount + 1;
				analyze[0].append(winners[i].host + " ");
				for (int r = 0; r < n; r++) {
					if (winners[i].host.equals(winners[r].host)) {
						count[i] = count[i] + 1;
					}
				}
			}
			if (count[i] > highcount) {
				highcount = count[i];
				a = i;
			}
		}
		for (int i = 0; i <= n - 1; i++) {
			if (analyze[0].indexOf(winners[i].geographic_location) == -1) {
				secondarycount = secondarycount + 1;
				analyze[0].append(winners[i].geographic_location + " ");
				for (int r = 0; r < n; r++) {
					if (winners[i].geographic_location
							.equals(winners[r].geographic_location)) {
						count[i] = count[i] + 1;
					}
				}
			}
			if (count[i] > highcount) {
				highcount = count[i];
				b = i;
			}
		}
		macroarray[bigcounter] = new macrodata();
		macroarray[bigcounter].host = winners[a].host;
		macroarray[bigcounter].hostcounter = count[a];
		macroarray[bigcounter].temphostcount = count[a];
		macroarray[bigcounter].temphoststor = winners[a].host;
		macroarray[bigcounter].country = winners[b].geographic_location;
		macroarray[bigcounter].countrycounter = count[b];
		macroarray[bigcounter].tempcountrycount = count[b];
		macroarray[bigcounter].tempcountrystor = winners[b].geographic_location;
		StringBuilder[] mutations = new StringBuilder[winners[0].seq.length()
				* winners.length];
		mutations[0]= new StringBuilder();
		mutations[0].append("Mutation\tLocation\tYear\tFrom\tTo");
		int o = 0;
		for (int l = 0; l < winners[0].seq.length()-1; l++) {
			for (int i = 0; i < winners.length-1; i++) {
				if (winners[i].seq.charAt(l) != winners[i + 1].seq.charAt(l)) {
					mutations[o]= new StringBuilder();
					mutations[o].append(l + "\t" + winners[i].year + "\t"
							+ winners[i].seq.charAt(l) + "\t"
							+ (winners[i + 1].seq.charAt(l)));
					o=o+1;

				}
			}
		}
		if (u == EvolTraj.sl - 1) {
			finaldata(macroarray, bigcounter, psa_results, high, medium,
					winners, mutations,o);
		} else {
			EvolTraj.print(psa_results, high, medium, macroarray, a, b,
					mostprobable, winners, mutations,o);
		}

	}


	public static void finaldata(macrodata[] macroarray, int bigcounter,
			StringBuilder[] psa_results, StringBuilder[] high,
			StringBuilder[] medium, ETstrain[] winners,
			StringBuilder[] mutations,int o) throws FileNotFoundException,
			UnsupportedEncodingException {

		StringBuilder[] templist = new StringBuilder[1];

		int highcount = 0;
		int a = 0;
		int b = 0;
		templist[0] = new StringBuilder();
		for (int i = 1; i <= bigcounter; i++) {
			if (templist[0].indexOf(macroarray[i].host) == -1) {
				templist[0].append(macroarray[i].host + "");
				for (int r = 1; r <= bigcounter; r++) {
					if (macroarray[i].host.equals(macroarray[r].host)) {
						macroarray[i].temphostcount = macroarray[i].temphostcount
								+ macroarray[r].hostcounter;
					}
				}
			}

			if (macroarray[i].temphostcount > highcount) {
				highcount = macroarray[i].temphostcount;
				a = i;
			}

		}

		for (int i = 1; i <= bigcounter; i++) {
			if (templist[0].indexOf(macroarray[i].country) == -1) {
				templist[0].append(macroarray[i].country + "");
				for (int r = 1; r <= bigcounter; r++) {
					if (macroarray[i].country.equals(macroarray[i].country)) {
						macroarray[i].tempcountrycount = macroarray[i].tempcountrycount
								+ macroarray[r].countrycounter;
					}
				}
			}
			if (macroarray[i].tempcountrycount > highcount) {
				highcount = macroarray[i].tempcountrycount;
				b = i;
			}

		}
		EvolTraj.print(psa_results, high, medium, macroarray, a, b,mostprobable, winners,mutations,o);

	}
}
