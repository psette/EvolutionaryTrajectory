package org.jcvi.psette;

import java.io.*;
import java.util.*;
import java.util.regex.Pattern;

/*Program written by: 
 Pietro Sette
 Supervised by
 Dr. Douglas Greer and Dr. Richard Schueremann*/
public class EvolTraj {
	public static String[] host, strain_name, header, sequence, headerparts,
			geographic_location, accesion_number;
	static int[] year, endPos, startPos;
	static double[] ancestors;
	public static String ref, query, fan, whichstrain;
	static int r, xcord, lxcord, counter, len, fxcord, sc, ssc, sl, bigcounter,
			u;
	static double ycord;
	public static float m, lowm;
	static double ls_y, ls_x, d, mg, fmg, smg, tmg, tx, tbx, tby, ty, confl,
			confu, sconfl, sconfu, tconfl, tconfu, std;
	private BufferedReader in;

	public EvolTraj(String filename) {
		readSequenceFromFile(filename);
	}

	void readSequenceFromFile(String file) { // Reads through the user specified
												// file
		List<String> description = new ArrayList<String>();
		List<String> seq = new ArrayList<String>();
		try {
			StringBuffer buffer_reader = new StringBuffer();
			in = new BufferedReader(new FileReader(file));
			String line = in.readLine();
			if (line == null)
				throw new IOException(file + " is an empty file");
			if (line.charAt(0) != '>')
				throw new IOException("First line of " + file
						+ " should start with '>'");
			else
				description.add(line);
			for (line = in.readLine().trim(); line != null; line = in
					.readLine()) {
				if (line.length() > 0 && line.charAt(0) == '>') {
					seq.add(buffer_reader.toString());
					buffer_reader = new StringBuffer();
					description.add(line);
				} else
					buffer_reader.append(line.trim());
			}
			if (buffer_reader.length() != 0)
				seq.add(buffer_reader.toString());
		} catch (IOException e) {
			System.out.println("Error when reading " + file);
			e.printStackTrace();
		}
		header = new String[description.size()];
		sequence = new String[seq.size()];
		for (int i = 0; i < seq.size(); i++) {
			header[i] = description.get(i);
			if (header.length < 4)
				;
			sequence[i] = seq.get(i);
		}
	}

	public String getHeader() {
		return header[0];
	}

	public String getSequence() {
		return sequence[0];
	}

	public String getHeader(int i) {
		return header[i];
	}

	public String getSequence(int i) {
		return sequence[i];
	}

	public int size() {
		return sequence.length;
	}

	public static void main(String[] args) throws Exception {
		String fn = "";
		if (args.length > 0)
			fn = args[0];
		else {
			System.out.print("Enter the name of the BLAST output directory:");

			File segments = new File((new BufferedReader(new InputStreamReader(
					System.in))).readLine());
			final long startTime = System.currentTimeMillis(); // Program timer

			File[] segment = segments.listFiles();
			macrodata[] macroarray = new macrodata[segment.length];

			for (u = 0; u < segment.length; u++) { // Finds all the files of type fasta in the directory
				if (segment[u].getName().endsWith(".fasta")) {
					bigcounter = bigcounter + 1;
					fn = segment[u].getPath();
					sl = segment.length;
					StringBuilder[] high = new StringBuilder[segment.length];
					StringBuilder[] psa_results = new StringBuilder[segment.length];
					StringBuilder[] medium = new StringBuilder[segment.length];

					high[u] = new StringBuilder("");
					psa_results[u] = new StringBuilder("");
					medium[u] = new StringBuilder("");
					whichstrain = segment[u].getName().substring(0,
							segment[u].getName().length() - 6);
					System.out.println("\n\t\t\t"
							+ segment[u].getName().substring(0,
									segment[u].getName().length() - 6));
					EvolTraj result = new EvolTraj(fn);
					ancestors = new double[result.size()];
					ETstrain[] ads = new ETstrain[result.size()];
					for (int p = 0; p < ads.length; p++) {
						ads[p] = new ETstrain();
					}
					arrayfiller(result, psa_results, high, medium);			// creates the arrays
					psalgo(0, result, psa_results, 1, 0, year, accesion_number,	
							startPos, endPos); // Calls the  first algorithm which creates a slope estimate
					psalgo2(high, medium, result, year, accesion_number, //Calls the second alogrithm which updates the slope estimate
							startPos, endPos, mg);
					TrajectoryCost.intializer(ads, macroarray, header,
							sequence, startPos, endPos, smg, bigcounter, u,
							psa_results, high, medium);	//Initializes the arrays for Trajectory Cost class
				}
			}

			final long endTime = System.currentTimeMillis();		//Timer
			System.out.println("\nTotal execution time:\n "
					+ (endTime - startTime) + " miliseconds\n\n\t\t"); // Program
																		// written
																		// by:
																		// \n\t\t\tPietro
																		// Sette\n\t\t\tSupervised
																		// by\n\t
																		// Dr.
																		// Douglas
																		// Greer
																		// and
																		// Dr.
																		// Richard
																		// Schueremann");
		}
	}

	public static void arrayfiller(EvolTraj result,
			StringBuilder[] psa_results, StringBuilder[] high,
			StringBuilder[] medium) {
		tx = 0;
		ty = 0;
		tby = 0;
		tbx = 0;
		tmg = 0;
		psa_results[u]
				.append("The following are most probable common ancestors: \n");
		high[u].append("The following are less probable, yet still possible ancestors: \n");
		medium[u]
				.append("The following are even less probable, yet still possible ancestors:\n");
		year = new int[result.size()]; // This makes all the data about the
		startPos = new int[result.size()];
		endPos = new int[result.size()];
		strain_name = new String[result.size()];
		host = new String[result.size()];
		accesion_number = new String[result.size()];
		geographic_location = new String[result.size()];
		for (int i = 0; i < result.size(); i++) {
			headerparts = header[i].split(Pattern.quote("|"));
			year[i] = Integer.parseInt(headerparts[2].substring(headerparts[2]
					.length() - 4));

			accesion_number[i] = headerparts[1];
			strain_name[i] = headerparts[0];
			host[i] = headerparts[4];
			geographic_location[i] = headerparts[3];
			len = result.getSequence(i).length();
			for (int x = 0; x < len; x++) {
				if (result.getSequence(i).charAt(x) != '-') { // Finds the
					// beginning and
					// ending
					// position of
					// each sequence
					startPos[i] = x;
					break;
				}
			}

			for (int x = len - 1; x >= 0; x--) {
				if (result.getSequence(i).charAt(x) != '-') {
					endPos[i] = x;
					break;
				}
			}

		}
	}

	public static void psalgo(
			int r,
			EvolTraj result, // This algorithm finds the lowest slope from the
			// reference point(which is centered around the
			// origin). After it finds the point with the
			// lowest slope, it re-centers about that point.
			// The algorithm
			// re-does the calculation so that the x and y
			// coordinates are years from the point and
			// genetic differences from the point. The point
			// is at the origin. It does this as needed
			// until the rightmost point is reached.
			StringBuilder[] psa_results, int counter, int sc, int[] year,
			String[] accesion_number, int[] startPos, int[] endPos) {

		sc = sc + 1;
		std = 0;
		if (counter == 0) {
			return; // if the counter is zero, the Algorithm has found the
			// rightmost point, this is the basecase
		} else { // if it is not zero, it keeps searching for the rightmost
					// point
			counter = 0;
			lowm = 99999;

			for (int i = 0; i < result.size(); i++) {
				int start = Math.max(startPos[r], startPos[i]);
				int end = Math.min(endPos[r], endPos[i]);
				xcord = year[r] - year[i];
				if (xcord > 0) {
					ref = result.getSequence(r).substring(start, end);
					query = result.getSequence(i).substring(start, end);
					ycord = 0;
					for (int l = 0; l < ref.length(); l++) {
						if (query.charAt(l) != ref.charAt(l)) {
							ycord = ycord + 1;
						}
					}
					d = (Math.log(1 - (4 / 3 * (ycord / query.length())))) * -3
							/ 4;

					if (xcord > lxcord && sc == 1) {
						lxcord = xcord;
						fmg = ycord / xcord;
					}

					if (xcord != 0 && ycord != 0) {
						m = (float) (ycord / xcord);
						if (m < lowm && m > 0) {
							counter = 1;
							lowm = m;
							ls_x = xcord;
							ls_y = ycord;
							fan = accesion_number[i];
							r = i; // R is the next point
						}
					}

				}
			}
			if (tx == 0 || lowm < 2.1650635094 * mg) { // This is eliminating
				// outliers, so if a
				// slope is 5 times the
				// standard deviation,
				// it is discarded. This
				// equation is given by
				// 5*sqrt(3/4*1/4)

				tx = tx + ls_x;
				ty = ty + ls_y;
				ancestors[sc] = ls_y / ls_x;

				if (xcord != 0 && ycord != 0) {
					if (psa_results[u].indexOf(fan) == -1) {

						psa_results[u].append(" " + fan); // This adds the
															// result
						// of
						// the algorithm (the
						// points
						// used to re-center) to
						// a
						// list for the purpose
						// of
						// printing out.

					}
				}
				mg = ty / tx; // Finds the average slope

			}
			confl = mg - 0.09375 - 0.3675 / sc + 0.0175781 / (sc - 1); // Creates
			// a
			// confidence
			// interval
			confu = mg - 0.09375 + 0.3675 / sc + 0.0175781 / (sc - 1);
			if (counter == 0) {
				for (int x = 0; x < sc; x++) {
					std = std + Math.pow((ancestors[x] - mg), 2);

					if (x == sc - 1) {
						std = Math.sqrt(std / sc);

					}
				}
			}
			psalgo(r, result, psa_results, counter, sc, year, accesion_number,
					startPos, endPos); // the algorithm calls itself to continue

		}

	}

	public static void psalgo2(StringBuilder[] high, StringBuilder[] medium,
			EvolTraj result, int[] year, String[] accesion_number,
			int[] startPos, int[] endPos, double mg) {
		// This function is very similar to the psalgo one, however, the main difference is that it computes if the Strain is in an ancestral candiate region
		ssc = 1;
		tx = 0;
		ty = 0;
		for (int i = 0; i < result.size(); i++) {
			int start = Math.max(startPos[0], startPos[i]);
			int end = Math.min(endPos[0], endPos[i]);
			xcord = year[0] - year[i];
			ref = result.getSequence(r).substring(start, end);
			query = result.getSequence(i).substring(start, end);
			ycord = 0;
			std = 0;
			for (int l = 0; l < ref.length(); l++) {
				if (query.charAt(l) != ref.charAt(l)) {
					ycord = ycord + 1;
				}
			}
			d = (Math.log(1 - (4 / 3 * (ycord / query.length())))) * -3 / 4;

			if (high[u].indexOf(accesion_number[i]) == -1) {
				if (medium[u].indexOf(accesion_number[i]) == -1) {

					if (xcord > 0) {
						if (xcord * mg >= ycord) {
							if (xcord * mg == ycord) {
								medium[u].append(" " + accesion_number[i]);
							} else {
								high[u].append(" " + accesion_number[i]);
							}
							ancestors[ssc - 1] = ycord / xcord;
							for (int x = 0; x < ssc; x++) {
								std = std + Math.pow((ancestors[x] - smg), 2);

								if (x == ssc - 1) {
									std = Math.sqrt(std / ssc);

								}
							}
							ssc = ssc + 1;
							tx = tx + xcord;
							ty = ty + ycord;
							smg = ty / tx; // Finds the average slope
							sconfl = smg
									+ Math.pow(std, 2)
									/ 2
									- 2.9957323
									* Math.sqrt(Math.pow(std, 2) / ssc
											+ Math.pow(std, 4) / (2 * ssc - 2));
							sconfu = smg
									+ Math.pow(std, 2)
									/ 2
									+ 2.99573232
									* Math.sqrt(Math.pow(std, 2) / ssc
											+ Math.pow(std, 4) / (2 * ssc - 2));

						}
					}
				}
			}
		}
	}

	public static void print(StringBuilder[] psa_results, StringBuilder[] high,
			StringBuilder[] medium, macrodata[] macroarray, int a, int b,
			StringBuilder mostprobable, ETstrain[] winners,
			StringBuilder[] mutations, int o) throws FileNotFoundException,
			UnsupportedEncodingException {
		// // //PsAlgo results:
		// System.out.println("Mutation rate estimate: "+ mg);
		// // // PsAlgo2 Results
		// System.out.println("Mutation rate estimate: " + smg);
		for (int d = 0; d < o; d++) {
			System.out.println(mutations[d]);
		}
		PrintWriter writer = new PrintWriter("ET" + whichstrain + ".fasta",
				"UTF-8");
		for (int i = 0; i < winners.length; i++) {
			writer.println(">" + winners[i].strain_name + "|"
					+ winners[i].accesion_number + "|" + winners[i].year + "|"
					+ winners[i].geographic_location + "|" + winners[i].host);
			writer.println(winners[i].seq + "\n");
		}
		writer.close();

		if (u == EvolTraj.sl - 1) {
			System.out
					.println("\n\n\n\n*************************************SUMMARY*************************************");
			System.out.println("Host of origin: " + macroarray[a].host);
			System.out.println("Region of origin: " + macroarray[b].country);
		}
	}
}