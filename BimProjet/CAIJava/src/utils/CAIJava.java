package utils;

/**
 * Title:        Project of visualizing codon frequencies
 * Description:
 * Copyright:    Copyright (c) 2002
 * Company:      IHES, France
 * @author Andrey Zinovyev, Alessandra Carbone
 * @version 1.0
 */

import java.io.*;
import java.util.*;
import java.lang.Integer;
import java.lang.System;

import org.biojava.bio.*;
import org.biojava.bio.symbol.*;
import org.biojava.bio.seq.*;
import org.biojava.bio.seq.io.*;

public class CAIJava {

	static int taskModifier = 0;
	static String addFeature = null;
	static String idFeature = null;
	static boolean calcGCContent = false;
	static boolean fastaFormat = false;
	static int calcCodonAnalysis = 0;
	static StringBuffer AdditionalInfoFile = new StringBuffer("");
	static float GTFAverLen[] = new float[64];
	static Vector TranslationalTable;
	static String externalWValues = null;
	static String externalWValues2 = null;
	static boolean checkForStability = false;
	static boolean includeDiplets = false;
	static int classNumber = -1;
	static int numberOfDecompositions = 0;
	static PrintStream out = System.out;
	static HashMap genomes = new HashMap();
	static boolean calculateRSCU = false;
	static int onlyOneGenome = -1;
	static int minimumRefsetSize = 0;

	public static void main(String[] args) {

		/*
		 * float gcu[] = new float[64]; Random rr = new Random(); for(int
		 * i=0;i<64;i++) gcu[i] = rr.nextFloat(); String cbb =
		 * calcCodonBiasCode(gcu,"UNIVERSAL"); System.out.println(cbb);
		 */

		try {
			if (args.length == 0) {
				throw new Exception("Use: CalcCodonFreq GenBankFile [options]");
			}
			Utils.out = out;

			// String sc = "atgtcagcctaa";
			// System.out.println(cutStartAndStop(sc,"UNIVERSAL"));
			// sc = "atgtcagccacctgg";
			// System.out.println(cutStartAndStop(sc,"UNIVERSAL"));

			// float f = Utils.GCContent3(sc);

			Vector Acids = new Vector();
			for (int i = 0; i < 64; i++) {
				SymbolList sDNA = DNATools.createDNA(Utils.TripletName(i));
				SymbolList Amin = RNATools.translate(RNATools.transcribe(sDNA));
				String sss = new String(Amin.seqString());
				if ((!Acids.contains(sss)) && (!sss.equals("*")))
					Acids.addElement(sss);
			}

			TranslationalTable = new Vector();
			for (int i = 0; i < 64; i++) {
				SymbolList sDNA = DNATools.createDNA(Utils.TripletName(i));
				SymbolList Amin = RNATools.translate(RNATools.transcribe(sDNA));
				if (!Amin.seqString().equals("*")) {
					int n = Acids.indexOf(Amin.seqString());
					TranslationalTable.addElement(new Integer(n));
				} else
					TranslationalTable.addElement(new Integer(-1));
			}

			int CodonsOrder[] = orderCodonUsageByAminoacid(Acids);
			Vector CodonUsages = new Vector();

			// for(int i=0;i<Acids.size();i++){
			// System.out.println(Acids.elementAt(i));
			// }

			TranslationTable TranTable = RNATools.getGeneticCode("UNIVERSAL");
			Alphabet Codons = TranTable.getSourceAlphabet();
			try {
				// for(int i=0;i<64;i++)
				// System.out.print(Utils.TripletName(CodonsOrder[i])+"\t");
				// System.out.print("\n");
				for (int i = 0; i < 64; i++) {
					SymbolList sDNA = DNATools.createDNA(Utils
							.TripletName(CodonsOrder[i]));
					Symbol sb = ((FiniteAlphabet) Codons).getSymbol((RNATools
							.transcribe(sDNA)).toList());
					String amin = TranTable.translate(sb).getName();
					// System.out.print(amin+"\t");
				}
				// System.out.print("\n");
			} catch (Exception e) {
				e.printStackTrace();
			}

			int windowSize = 120;
			int evPosition = 21, sAnnotbeg = -1, sAnnotend = -1;
			int numberOfIterations = 0;
			int maximumGeneLength = -1;
			int minimumGeneLength = -1;
			int BayezRatio = 0;
			String fileName = "";

			for (int i = 1; i <= args.length; i++)
				if (i < args.length) {
					if (args[i].length() > 1)
						if (args[i].equals("-ew")) {
							externalWValues = args[i + 1];
							continue;
						}
					if (args[i].equals("-ew2")) {
						externalWValues2 = args[i + 1];
						continue;
					}
					if ((args[i].length() > 2)
							&& (args[i].substring(0, 3).equals("-cs"))) {
						checkForStability = true;
						if (args[i].length() > 3)
							classNumber = Integer
									.parseInt(args[i].substring(3));
						;
						continue;
					}
					if (args[i].equals("-w"))
						windowSize = Integer.parseInt(args[i + 1]);
					if (args[i].equals("-e"))
						evPosition = Integer.parseInt(args[i + 1]);
					if (args[i].equals("-f"))
						fileName = new String(args[i + 1]);
					if (args[i].equals("-i"))
						numberOfIterations = Integer.parseInt(args[i + 1]);
					if (args[i].equals("-k"))
						taskModifier = Integer.parseInt(args[i + 1]);
					if (args[i].equals("-a"))
						AdditionalInfoFile = new StringBuffer(args[i + 1]);
					if (args[i].equals("-m"))
						maximumGeneLength = Integer.parseInt(args[i + 1]);
					if (args[i].equals("-min"))
						minimumGeneLength = Integer.parseInt(args[i + 1]);
					if (args[i].equals("-minref"))
						minimumRefsetSize = Integer.parseInt(args[i + 1]);
					if (args[i].equals("-c"))
						calcCodonAnalysis = Integer.parseInt(args[i + 1]);
					if (args[i].equals("-t"))
						addFeature = new String(args[i + 1]);
					if (args[i].equals("-id"))
						idFeature = new String(args[i + 1]);
					;
					if (args[i].equals("-s"))
						fastaFormat = true;
					if (args[i].equals("-g"))
						calcGCContent = true;
					if (args[i].equals("-d"))
						includeDiplets = true;
					if (args[i].equals("-rscu"))
						calculateRSCU = true;
					if (args[i].equals("-decomp"))
						numberOfDecompositions = Integer.parseInt(args[i + 1]);
				}

			/*
			 * String ss1 = "aaaaaaaaaaaaaag"; String ss2 =
			 * "aaaagaataagaaggaag"; float wv[] = calculateCodonUsage(ss1);
			 * float ks = 0; for(int i=0;i<64;i++) ks+=wv[i];
			 * GenerateCAIWValuesTable(wv,"UNIVERSAL"); float cus[] =
			 * calculateCodonUsage(ss2); ks = 0; for(int i=0;i<64;i++)
			 * ks+=cus[i];
			 * 
			 * float cai1 = CalculateCAIValue(ss2,wv); float cai2 =
			 * CalculateCAIValue(cus,wv);
			 */

			// Prepare gene file

			HashMap GeneAddInfo = new HashMap();
			HashMap GenomeLengthes = new HashMap();
			String StringDelimiters = new String("\t");
			int GeneAddInfoCount = 0;
			out.println("AdditionalInfoFile = \"" + AdditionalInfoFile + "\"");
			if (!AdditionalInfoFile.toString().equals("")) {
				out.println("Trying to read...");
				LineNumberReader lr = new LineNumberReader(new FileReader(
						AdditionalInfoFile.toString()));
				String s;
				s = lr.readLine();
				lr.close();
				lr = new LineNumberReader(new FileReader(
						AdditionalInfoFile.toString()));
				StringTokenizer st = new StringTokenizer(s, StringDelimiters);
				int j = 0;
				while ((st.hasMoreTokens())) {
					st.nextToken();
					j++;
				}
				GeneAddInfoCount = j;

				// lr.reset();
				while (((s = lr.readLine()) != null)) {
					st = new StringTokenizer(s, StringDelimiters);
					// out.println("AddInfofile: "+s.substring(1,10)+"\t"+st.countTokens());
					String[] sInfo = new String[j];
					int k = 0;
					while ((st.hasMoreTokens())) {
						String sss = st.nextToken();
						if ((k == 1) && (sss.length() == 1))
							sss = "0" + sss;
						if (k < GeneAddInfoCount)
							sInfo[k] = sss.trim().replace(' ', '_');
						else
							out.println("AddInfofile: out of bound " + s);
						k++;
					}
					if (GeneAddInfo.containsKey(sInfo[0])) {
						String sm[] = (String[]) GeneAddInfo.get(sInfo[0]);
						for (int i = 1; i < sm.length; i++)
							sm[i] = sm[i] + ";" + sInfo[i];
					} else
						GeneAddInfo.put(sInfo[0], sInfo);
				}
			}

			// for(Iterator i=GeneAddInfo.keySet().iterator();i.hasNext();){
			// String sm[] = (String[])GeneAddInfo.get((String)(i.next()));
			// for (int ii = 0; ii < sm.length; ii++) out.print(sm[ii]+" ");
			// out.println();
			// }
			//

			// read external WValues
			float extWValues[] = new float[64];
			if (externalWValues != null) {
				LineNumberReader lr = new LineNumberReader(new FileReader(
						externalWValues));
				int k = 0;
				while (k < 64) {
					String s = lr.readLine();
					StringTokenizer st = new StringTokenizer(s, "\t ,");
					// while ((st.hasMoreTokens())) {
						//String ss = st.nextToken();
						//String tn = Utils.TripletName(k, "t", "c", "a", "g");
						//int i = Utils.TripletNumber(tn);
						//extWValues[i] = Float.parseFloat(ss);
						String tn = st.nextToken();
						int i = Utils.TripletNumber(tn);
						String ss = st.nextToken();
						extWValues[i] = Float.parseFloat(ss);
						k++;
					//}
				}
			}
			System.out.println("extWValues : ");
			Utils.printDistribution(extWValues, "t", "c", "g", "a");

			// read external WValues 2
			float extWValues2[] = new float[64];
			if (externalWValues2 != null) {
				LineNumberReader lr = new LineNumberReader(new FileReader(
						externalWValues2));
				int k = 0;
				while (k < 64) {
					String s = lr.readLine();
					StringTokenizer st = new StringTokenizer(s, "\t ,");
					while ((st.hasMoreTokens())) {
						String ss = st.nextToken();
						String tn = Utils.TripletName(k, "t", "c", "g", "a");
						int i = Utils.TripletNumber(tn);
						extWValues2[i] = Float.parseFloat(ss);
						k++;
					}
				}
			}
			Utils.printDistribution(extWValues2, "t", "c", "g", "a");

			// Random class file
			int randomclassnumber = 0;
			String randomclass[][] = new String[1000][4];
			out.println("Looking for randomclass.txt...");
			File ff = new File("randomclass.txt");
			if (ff.exists()) {
				LineNumberReader lr = new LineNumberReader(new FileReader(
						"randomclass.txt"));
				String sr = null;
				while ((sr = lr.readLine()) != null) {
					StringTokenizer st = new StringTokenizer(sr, "\t");
					randomclass[randomclassnumber][0] = st.nextToken();
					randomclass[randomclassnumber][1] = st.nextToken();
					randomclass[randomclassnumber][2] = st.nextToken();
					randomclass[randomclassnumber][3] = st.nextToken();
					randomclassnumber++;
				}
				out.println(randomclassnumber
						+ " genes are found in randomclass.txt");
			}

			HashMap GeneNoteTable = new HashMap();

			PrintStream fout;
			if (fileName.equals(""))
				fout = out;
			else
				fout = new PrintStream(new FileOutputStream(fileName));

			File GenBankFile = new File(args[0]);
			out.println("Loading sequence...");
			BufferedReader eReader = new BufferedReader(new InputStreamReader(
					new FileInputStream(GenBankFile)));
			SequenceIterator seqI = null;
			if (!fastaFormat)
				seqI = SeqIOTools.readGenbank(eReader);
			else
				seqI = SeqIOTools.readFastaDNA(eReader);
			out.println("Loaded...");

			int NumberOfPoints = 0;

			int sAnnBegM = sAnnotbeg;
			int sAnnEndM = sAnnotend;
			int TF[];
			float GTF[] = new float[64];

			// Calculating number of points and CAIW tables
			Vector ListOfGenes = new Vector();
			Vector CAIValues = new Vector();

			int numgenomes = 0;
			if (!fastaFormat)
				while (seqI.hasNext()) {
					out.println("Getting seq...");
					Sequence seq = seqI.nextSequence();
					out.println("Got " + seq.getName());
					numgenomes++;
					if (numgenomes > 1)
						onlyOneGenome = -1;
					else
						onlyOneGenome = seq.length();
					Vector genes = new Vector();
					int numCDSFeat = 0;
					for (Iterator i = seq.features(); i.hasNext();) {
						StrandedFeature Ft = (StrandedFeature) i.next();
						String sseq = Ft.getSymbols().seqString();
						if (((Ft.getType().equals("CDS")))
								&& (sseq.length() > 10)) {
							numCDSFeat++;
							try {
								String ss = new String(Ft.getSymbols()
										.seqString());
								if (maximumGeneLength != -1)
									if (ss.length() > maximumGeneLength)
										continue;
								if (minimumGeneLength != -1)
									if (ss.length() < minimumGeneLength)
										continue;

								TF = Utils.calcTripletFreq(
										cutStartAndStop(ss, "UNIVERSAL"), 1,
										cutStartAndStop(ss, "UNIVERSAL")
												.length());
								// k=CalcCodonUsageForSeq(Ft.getSymbols(),fr0,fr1,fr2);
								// ss = new String(Ft.getSymbols().seqString());
								// if (ss.length()/3.0!=(int)(ss.length()/3.0))
								// out.println("Mistake!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
								// out.println(Ft.getAnnotation().getProperty(new
								// String("gene")).toString()+"; "+ss.length()+" symbols");
								// String check = new
								// String(Ft.getAnnotation().getProperty("gene").toString());
								// GeneticCodes GCode = new GeneticCodes();
								// out.println(GCode.translate(GCode.transcribe(Ft.getSymbols())).seqString());
								// out.println(ss);
								for (int j = 0; j < 64; j++)
									GTF[j] += TF[j];
								NumberOfPoints++;
								ListOfGenes.addElement(Ft);
								float CAIV[] = new float[numberOfIterations];
								CAIValues.addElement(CAIV);

								String geneName = getGeneName(Ft);
								out.printf("-----------------------> %s\n",geneName);
								String geneID = getGeneID(Ft);

								if (Ft.getAnnotation().containsProperty("gene")) {
									genes.add(geneID);
								} else {
									genes.add(geneID);
								}

							} catch (Exception e) {
								out.print("\n" + e + " : ");
								out.println(Ft);
							}
						}

						if ((Ft.getType().equals("gene"))) {

							String geneName = getGeneName(Ft);
							String geneID = getGeneID(Ft);

							if ((Ft.getAnnotation().containsProperty("note"))
									&& (Ft.getAnnotation()
											.containsProperty("gene"))) {
								String geneNote = Ft.getAnnotation()
										.getProperty("note").toString();
								GeneNoteTable.put(geneName, geneNote);
							} else {
								GeneNoteTable.put(geneName, geneName);
							}

							if (Ft.getAnnotation().containsProperty("gene")) {
								genes.add(geneID);
							} else {
								genes.add(geneID);
							}

						}
					}
					out.println(numCDSFeat + " CDS features found");
					genomes.put(seq.getName(), genes);
					GenomeLengthes
							.put(seq.getName(), new Integer(seq.length()));
				}

			out.println("Gene Note table: ");
			Set se = GeneNoteTable.keySet();
			for (Iterator i = se.iterator(); i.hasNext();) {
				String ss = (String) i.next();
				out.println(ss + "\t" + GeneNoteTable.get(ss));
			}

			if (fastaFormat)
				while (seqI.hasNext()) {
					Sequence seq = seqI.nextSequence();
					out.println("Got" + seq.getName());
					String ss = seq.seqString();
					if (ss.indexOf("N") == -1) {
						if (maximumGeneLength != -1)
							if (ss.length() > maximumGeneLength)
								continue;
						if (minimumGeneLength != -1)
							if (ss.length() < minimumGeneLength)
								continue;

						TF = Utils.calcTripletFreq(
								cutStartAndStop(ss, "UNIVERSAL"), 1,
								cutStartAndStop(ss, "UNIVERSAL").length());
						for (int j = 0; j < 64; j++)
							GTF[j] += TF[j];
						NumberOfPoints++;
						ListOfGenes.addElement(seq);
						float CAIV[] = new float[numberOfIterations];
						CAIValues.addElement(CAIV);
					}
				}

			int sum = 0;
			for (int i = 0; i < 64; i++)
				sum += GTF[i];
			for (int i = 0; i < 64; i++)
				GTF[i] = GTF[i] / sum;

			out.print("CodonUsage" + " sphere ");
			for (int ii = 0; ii < 64; ii++)
				out.print(GTF[ii] + " ");
			out.println("0.4 255 0 0");

			float GTFC[] = new float[64];
			System.arraycopy(GTF, 0, GTFC, 0, 64);
			CodonUsages.addElement(GTFC);

			Utils.printDistribution(GTF);
			float GTFAver[] = new float[64];
			for (int i = 0; i < 64; i++)
				GTFAver[i] = GTF[i];
			GenerateCAIWValuesTable(GTF, "UNIVERSAL");
			GenerateCAIWValuesTableMethodAverage(GTFAver, "UNIVERSAL");
			GenerateCAIWValuesSuppressingRareCodons(GTFAver, "UNIVERSAL",
					ListOfGenes);
			out.println();
			// Utils.printDistribution(GTF);
			// Utils.printDistribution(GTF,"t","c","g","a");
			out.println();
			// Utils.printDistribution(GTFAver);

			out.println("Number of points = " + NumberOfPoints);

			// --- Calculating all iterations
			Random r = new Random();
			float FirstCAI[] = new float[ListOfGenes.size()];
			float first_p = 100.0f;
			float cais[] = new float[ListOfGenes.size()];
			for (int i = 0; i < numberOfIterations; i++) {
				if (i == 0) {
					for (int j = 0; j < ListOfGenes.size(); j++) {
						if (!fastaFormat) {
							if (!checkForStability)
								cais[j] = CalculateCAIValue(
										((StrandedFeature) ListOfGenes.elementAt(j))
												.getSymbols().seqString(), GTF);
							else {
								cais[j] = r.nextFloat();
								first_p = 1f;
								if (classNumber != -1) {
									String name = "**";
									String seq = "";
									try {
										name = ((StrandedFeature) ListOfGenes
												.elementAt(j)).getAnnotation()
												.getProperty("gene").toString();
										seq = ((StrandedFeature) ListOfGenes
												.elementAt(j)).getSymbols()
												.seqString();
									} catch (Exception e) {
									}
									;
									cais[j] = 0;
									for (int iii = 0; iii < randomclassnumber; iii++)
										if (randomclass[iii][classNumber]
												.equals(name + "_"
														+ seq.length())) {
											cais[j] = 1;
											out.println(Integer.toString(j + 1)
													+ "\t" + name);
										}
								}
							}
						} else {
							if (!checkForStability)
								cais[j] = CalculateCAIValue(
										((Sequence) ListOfGenes.elementAt(j))
												.seqString(),
										GTF);
							else {
								cais[j] = r.nextFloat();
								first_p = 1f;
							}
						}
					}
					// for (int ii = 0; ii < cais.length; ii++)
					// out.print(cais[ii]+" ");
					// out.println();
					for (int hh = 0; hh < ListOfGenes.size(); hh++)
						FirstCAI[hh] = cais[hh];
				}

				int prop = (int) (first_p / Math.pow(2, i + 1));
				if (prop < 1)
					prop = 1;

				out.println("Iteration " + i);
				float CodUs[] = new float[64];
				MakeCAIIteration(ListOfGenes, cais, prop, CodUs);

				out.print("CodonUsageIt" + i + " sphere ");
				for (int ii = 0; ii < 64; ii++)
					out.print(CodUs[ii] + " ");
				out.println("0.2 0 255 0");

				CodonUsages.addElement(CodUs);

				// for (int ii = 0; ii < cais.length; ii++)
				// out.print(cais[ii]+" ");
				// out.println();

				for (int j = 0; j < ListOfGenes.size(); j++) {
					float mas[] = (float[]) CAIValues.elementAt(j);
					mas[i] = cais[j];
				}

			}// ---

			// printing codon usage
			PrintStream cuout = new PrintStream(new FileOutputStream(
					"codonusage"));
			for (int j = 0; j < 61; j++) {
				SymbolList sDNA = DNATools.createDNA(Utils.TripletName(
						CodonsOrder[j], "a", "t", "g", "c"));
				String trAcid = (RNATools.translate(RNATools.transcribe(sDNA)))
						.seqString();
				cuout.print(trAcid + "\t"
						+ Utils.TripletName(CodonsOrder[j], "a", "t", "g", "c")
						+ "\t");
				for (int i = 0; i < CodonUsages.size(); i++) {
					float cU[] = (float[]) (CodonUsages.elementAt(i));
					cuout.print(cU[CodonsOrder[j]] + "\t");
				}
				cuout.println();
			}
			cuout.close();

			// printing W-values
			PrintStream wout = new PrintStream(new FileOutputStream("wvalues"));
			for (int j = 0; j < 61; j++) {
				SymbolList sDNA = DNATools.createDNA(Utils
						.TripletName(CodonsOrder[j]));
				String trAcid = (RNATools.translate(RNATools.transcribe(sDNA)))
						.seqString();
				wout.print(trAcid + "\t" + Utils.TripletName(CodonsOrder[j])
						+ "\t");
				for (int i = 0; i < CodonUsages.size(); i++) {
					float cU[] = (float[]) (CodonUsages.elementAt(i));
					float cp[] = new float[64];
					for (int k = 0; k < 64; k++)
						cp[k] = cU[k];
					GenerateCAIWValuesTable(cp, "UNIVERSAL");
					wout.print(cp[CodonsOrder[j]] + "\t");
				}
				wout.println();
			}
			wout.close();

			// table for final W-values

			// cuout = new PrintStream(new FileOutputStream("wvalues"));
			int ncu = CodonUsages.size();
			float GTFFinal[] = new float[64];
			// float GTFFinalForClassicCAI[] = new float[64];
			float CodonUsageFinal[] = (float[]) (CodonUsages.elementAt(ncu - 1));
			for (int i = 0; i < 64; i++) {
				GTFFinal[i] = CodonUsageFinal[i];
			}
			GenerateCAIWValuesTable(GTFFinal, "UNIVERSAL");
			Utils.printDistribution(GTFFinal);
			out.println();
			Utils.printDistribution(GTFFinal, "t", "c", "g", "a");
			if (externalWValues != null) {
				out.println("External W-table");
				Utils.printDistribution(extWValues, "t", "c", "g", "a");
			}
			if (externalWValues2 != null) {
				out.println("External W-table");
				Utils.printDistribution(extWValues2, "t", "c", "g", "a");
			}

			float CodonUsageInitial[] = (float[]) (CodonUsages.elementAt(0));
			float WTableInitial[] = new float[64];
			for (int i = 0; i < 64; i++) {
				WTableInitial[i] = CodonUsageInitial[i];
			}
			GenerateCAIWValuesTable(WTableInitial, "UNIVERSAL");

			// printing all files
			int dpos = args[0].indexOf(".");
			String basename = args[0].substring(0, dpos);
			Utils.printDistributionToFile(basename + ".cui", CodonUsageInitial,
					"t", "c", "g", "a");
			Utils.printDistributionToFile(basename + ".cuf", CodonUsageFinal,
					"t", "c", "g", "a");
			Utils.printDistributionToFile(basename + ".wti", WTableInitial,
					"t", "c", "g", "a");
			Utils.printDistributionToFile(basename + ".wtf", GTFFinal, "t",
					"c", "g", "a");

			// main table
			int addf = 0;
			if (GeneAddInfoCount != 0)
				addf += GeneAddInfoCount - 1;
			if (addFeature != null)
				addf++;
			if (idFeature != null)
				addf++;
			if (calcGCContent)
				addf += 4;
			if (externalWValues != null)
				addf++;
			if (externalWValues2 != null)
				addf++;
			if (includeDiplets)
				addf += 32;
			if (!fastaFormat)
				addf += 5;
			addf += numberOfDecompositions * 2;
			fout.println((94 + (numberOfIterations + addf)) + " "
					+ (int) NumberOfPoints);
			for (int i = 0; i < 64; i++)
				fout.println(Utils.TripletName(i) + " FLOAT");
			for (int i = 0; i < 20; i++)
				fout.println("AA_" + (String) Acids.elementAt(i) + " FLOAT");
			if (includeDiplets) {
				for (int i = 0; i < 4; i++)
					fout.println(Utils.BaseName(i) + " FLOAT");
				fout.println("a1 FLOAT");
				fout.println("c1 FLOAT");
				fout.println("g1 FLOAT");
				fout.println("t1 FLOAT");
				fout.println("a2 FLOAT");
				fout.println("c2 FLOAT");
				fout.println("g2 FLOAT");
				fout.println("t2 FLOAT");
				fout.println("a3 FLOAT");
				fout.println("c3 FLOAT");
				fout.println("g3 FLOAT");
				fout.println("t3 FLOAT");
				for (int i = 0; i < 16; i++)
					fout.println(Utils.DipletName(i) + " FLOAT");
			}
			fout.println("CB_CODE" + " STRING");
			fout.println("GENOME" + " STRING");
			fout.println("GENENAME" + " STRING");
			fout.println("GENE_ID" + " STRING");
			if (!fastaFormat) {
				fout.println("COMP" + " FLOAT");
				fout.println("POSITION" + " FLOAT");
				fout.println("START" + " FLOAT");
				fout.println("END" + " FLOAT");
				fout.println("STRAND" + " FLOAT");
			}
			if (numberOfDecompositions > 0) {
				for (int kk = 0; kk < numberOfDecompositions; kk++) {
					fout.println("DCMPCOEF" + (kk + 1) + " FLOAT");
					fout.println("DCMPCAI" + (kk + 1) + " FLOAT");
				}
			}
			if (addFeature != null)
				fout.println("ADDF" + " STRING");
			if (idFeature != null)
				fout.println(idFeature.toUpperCase() + " STRING");
			if (calcGCContent) {
				fout.println("GC_CONT" + " FLOAT");
				fout.println("GC_CONT3" + " FLOAT");
				fout.println("GC_SKEW" + " FLOAT");
				fout.println("AT_SKEW" + " FLOAT");
			}
			if (externalWValues != null)
				fout.println("CAIEXT" + " FLOAT");
			if (externalWValues2 != null)
				fout.println("CAIEXT2" + " FLOAT");
			fout.println("CAI" + " FLOAT");
			fout.println("CAICLASS" + " FLOAT");
			fout.println("LN_CAICLASS" + " FLOAT");
			fout.println("CAI1" + " FLOAT");
			fout.println("CAI2" + " FLOAT");
			fout.println("GENELENGTH" + " FLOAT");
			for (int i = 0; i < numberOfIterations; i++)
				fout.println("CAI_IT" + (i + 1) + " FLOAT");
			for (int i = 1; i < GeneAddInfoCount; i++)
				fout.println("ADD" + (i) + " STRING");

			/*
			 * eReader = new BufferedReader( new InputStreamReader(new
			 * FileInputStream(GenBankFile)));
			 * out.println("Loading sequence..."); seqI =
			 * SeqIOTools.readGenbank(eReader); out.println("Loaded...");
			 */

			// Here fill table for all genomes
			// while(seqI.hasNext()) {
			// Sequence seq = seqI.nextSequence();

			// for(Iterator i=seq.features(); i.hasNext(); )

			Annotatable Ft = null;

			Date SO = null;
			Date EO = null;

			StringBuffer sOut = new StringBuffer();

			float correlations[][] = new float[8][ListOfGenes.size()];
			// 0 - CAICLASS
			// 1 - STRAND
			// 2 - GC_CONT
			// 3 - GC_CONT3
			// 4 - GC_SKEW
			// 5 - AT_SKEW
			// 6 - GENELENGTH
			// 7 - CAIEXT
			String infos[][] = new String[2][ListOfGenes.size()];

			for (int i = 0; i < ListOfGenes.size(); i++) {
				// Date SO1 = new Date();
				// Feature Ft = (Feature) i.next();
				Ft = (Annotatable) ListOfGenes.elementAt(i);
				// if ((Ft.getType().equals("CDS")))
				{
					try {
						// k=CalcCodonUsageForSeq(Ft.getSymbols(),fr0,fr1,fr2);
						String ss = null;
						if (!fastaFormat)
							ss = ((StrandedFeature) Ft).getSymbols()
									.seqString();
						else
							ss = ((Sequence) Ft).seqString();
						// SO = new Date();
						ss = ss.substring(0, ss.length());
						float caiv = FirstCAI[i]; // CalculateCAIValue(ss,GTF);
						float caivclassic = CalculateCAIValue(ss, GTFFinal);
						float caiv1 = CalculateCAIValue(ss, GTFAver);
						float caiv2 = CalculateCAIValue(ss, GTFAverLen);
						// EO = new Date();
						// out.println("CAI calc: "+(EO.getTime()-SO.getTime()));
						sOut = new StringBuffer("");

						String sss = cutStartAndStop(ss, "UNIVERSAL");
						TF = Utils.calcTripletFreq(sss, 1, sss.length());
						if (calculateRSCU) {
							float tmpF[] = new float[64];
							for (int kk = 0; kk < 64; kk++)
								tmpF[kk] = TF[kk];
							recalculateCodonsForRSCU(tmpF, "UNIVERSAL");
							for (int j = 0; j <= 63; j++)
								sOut.append((float) ((float) tmpF[j]) + " ");
						} else
							for (int j = 0; j <= 63; j++)
								sOut.append((float) ((float) TF[j] / (sss
										.length() / 3)) + " ");

						float GeneCodonUsage[] = new float[64];
						for (int jj = 0; jj < 64; jj++)
							GeneCodonUsage[jj] = (float) TF[jj]
									/ (sss.length() / 3);

						// SO = new Date();
						TF = calcAminoacidFreq(sss, Acids, 1, sss.length());
						// EO = new Date();
						// out.println("Calc Freq: "+(EO.getTime()-SO.getTime()));

						for (int j = 0; j < 20; j++)
							sOut.append((float) ((float) TF[j] / (sss.length() / 3))
									+ " ");

						if (includeDiplets) {
							float baseF[] = Utils.calcBaseFreqF(ss, 1,
									ss.length());
							float diplF[] = Utils.calcDipletFreqFOverlapping(
									ss, 1, ss.length());
							float aphase[] = Utils.calcNinPhase(GeneCodonUsage,
									"a");
							float cphase[] = Utils.calcNinPhase(GeneCodonUsage,
									"c");
							float gphase[] = Utils.calcNinPhase(GeneCodonUsage,
									"g");
							float tphase[] = Utils.calcNinPhase(GeneCodonUsage,
									"t");

							for (int k = 0; k < 4; k++)
								sOut.append(baseF[k] + " ");

							sOut.append(aphase[0] + " ");
							sOut.append(cphase[0] + " ");
							sOut.append(gphase[0] + " ");
							sOut.append(tphase[0] + " ");
							sOut.append(aphase[1] + " ");
							sOut.append(cphase[1] + " ");
							sOut.append(gphase[1] + " ");
							sOut.append(tphase[1] + " ");
							sOut.append(aphase[2] + " ");
							sOut.append(cphase[2] + " ");
							sOut.append(gphase[2] + " ");
							sOut.append(tphase[2] + " ");

							for (int k = 0; k < 16; k++)
								sOut.append(diplF[k] + " ");
						}

						String cb = calcCodonBiasCode(GeneCodonUsage,
								"UNIVERSAL");
						sOut.append("_" + cb + "_");

						// for (int j=0; j<=63; j++) sOut.append(TF[j]+" ");

						String geneName = getGeneName(Ft);
						String geneID = getGeneID(Ft);
						String geneNote = "";

						// if (Ft.getAnnotation().containsProperty("gene"))
						// geneName =
						// Ft.getAnnotation().getProperty("gene").toString();
						if (Ft.getAnnotation().containsProperty("note"))
							geneNote = Ft.getAnnotation().getProperty("note")
									.toString();

						String genomeName = getGenomeName(geneName + "_"
								+ ss.length());
						sOut.append("\"" + genomeName + "\" ");

						sOut.append(" \"" + geneName + "\" ");
						infos[0][i] = geneName;
						sOut.append(" \"" + geneName + "_" + ss.length()
								+ "\" ");

						if (!fastaFormat) {
							int cmp = ((StrandedFeature) Ft).getStrand()
									.getValue();
							int ps = Utils
									.getGenePosition((StrandedFeature) Ft);
							sOut.append(" " + cmp + " ");
							sOut.append(" " + ps + " ");
							sOut.append(" "
									+ ((StrandedFeature) Ft).getLocation()
											.getMin() + " ");
							sOut.append(" "
									+ ((StrandedFeature) Ft).getLocation()
											.getMax() + " ");
							Integer gl = (Integer) GenomeLengthes
									.get(genomeName);
							if ((gl != null) || (onlyOneGenome != -1)) {
								int str = cmp;
								if (gl != null)
									if (ps > gl.intValue() / 2)
										str = -cmp;
									else if (ps > onlyOneGenome / 2)
										str = -cmp;
								sOut.append(" " + str + " ");
								correlations[1][i] = str;
							} else
								sOut.append(" 0 ");
						}

						if (addFeature != null) {
							String af = new String("");
							// if (Ft.getAnnotation().containsProperty("gene"))
							if (Ft.getAnnotation().containsProperty(addFeature))
								af = new String(Ft.getAnnotation()
										.getProperty(addFeature).toString());
							af = af.replace(' ', '_');
							if (af.length() > 100)
								af = af.substring(0, 100);
							af = " \"" + af + "\" ";
							sOut.append(af);
							infos[1][i] = af;
						}
						if (idFeature != null) {
							String af = new String("");
							if (Ft.getAnnotation().containsProperty(idFeature))
								af = new String(Ft.getAnnotation()
										.getProperty(idFeature).toString());
							af = af.replace(' ', '_');
							if (af.length() > 100)
								af = af.substring(0, 100);
							af = " \"" + af + "\" ";
							sOut.append(af);
						}

						if (calcGCContent) {
							// SO = new Date();
							float gc = Utils.GCContent(ss);
							float gc3 = Utils.GCContent3(ss);
							if (gc3 > 1)
								out.println("gc3>1: " + ss);
							float gcsk = Utils.GCSkew(ss);
							float atsk = Utils.ATSkew(ss);

							correlations[2][i] = gc;
							correlations[3][i] = gc3;
							correlations[4][i] = gcsk;
							correlations[5][i] = atsk;

							// EO = new Date();
							// out.println("GC Cont: "+(EO.getTime()-SO.getTime()));

							sOut.append(" " + gc + " " + gc3 + " " + gcsk + " "
									+ atsk + " ");
						}
						if (externalWValues != null) {
							float caivext = CalculateCAIValue(ss, extWValues);
							sOut.append(caivext + " ");
							correlations[7][i] = caivext;
							System.out.println(caivext);
						}

						if (externalWValues2 != null) {
							float caivext = CalculateCAIValue(ss, extWValues2);
							sOut.append(caivext + " ");
						}

						correlations[0][i] = caivclassic;

						sOut.append(caiv + " ");
						sOut.append(caivclassic + " ");
						sOut.append(Math.log(caivclassic) + " ");
						sOut.append(caiv1 + " ");
						sOut.append(caiv2 + " ");
						sOut.append(ss.length() + " ");
						correlations[6][i] = ss.length();

						for (int ii = 0; ii < numberOfIterations; ii++) {
							float mas[] = (float[]) CAIValues.elementAt(i);
							sOut.append(mas[ii] + " ");
						}
						if (Ft.getAnnotation().containsProperty("gene"))
							geneName = new String(Ft.getAnnotation()
									.getProperty(new String("gene")).toString());
						// if (Ft.getAnnotation().containsProperty("note"))
						// {
						// geneNote = new
						// String(Ft.getAnnotation().getProperty(new
						// String("note")).toString());
						// out.println(geneName+" "+geneNote);
						// }

						if (GeneNoteTable.containsKey(geneName))
							geneNote = (String) (GeneNoteTable.get(geneName));
						else
							geneNote = new String("");

						if ((GeneAddInfo.containsKey(geneName))
								|| (GeneAddInfo.containsKey(geneNote))) {
							String sm[];
							if (GeneAddInfo.containsKey(geneName))
								sm = (String[]) (GeneAddInfo.get(geneName));
							else {
								sm = (String[]) (GeneAddInfo.get(geneNote));
								// out.println("!!!"+geneNote);
							}
							// GeneAddInfo.get()
							for (int ii = 1; ii < GeneAddInfoCount; ii++) {
								// out.println(sm[ii].substring(1,2)+"//"+sm[ii].substring(0,1));
								if (sm[ii] == null) {
									sOut.append("\" \"");
									continue;
								}
								if ((sm[ii].length() != 0)
										&& (sm[ii].substring(0, 1).equals("\"")))
									sOut.append(sm[ii] + " ");
								else
									sOut.append("\"" + sm[ii] + "\" ");
							}
						} else
							for (int ii = 1; ii < GeneAddInfoCount; ii++)
								sOut.append("\"\" ");
						fout.println(sOut.toString());
						// Date EO1 = new Date();
						// out.println("Total("+ss.length()+"): "+(EO1.getTime()-SO1.getTime()));

					} catch (Exception e) {
						/**
						 * out.print("\n"+e+" : "); out.println(Ft);
						 **/
						// out.println(sOut.toString());
						out.print("\n" + e + " : ");
						e.printStackTrace();
					}
				}
			}
			fout.close();

			// Calc correlations and t-value
			out.println("\n// 0 - CAICLASS \n// 1 - STRAND\n// 2 - GC_CONT\n// 3 - GC_CONT3\n// 4 - GC_SKEW\n// 5 - AT_SKEW\n// 6 - GENELENGTH\n// 7 - CAIEXT\n");
			out.println("---");
			out.println(calcTValue(correlations[0], correlations[1]));
			out.println(calcCorrelationCoeff(correlations[0], correlations[2]));
			out.println(calcCorrelationCoeff(correlations[0], correlations[3]));
			out.println(calcCorrelationCoeff(correlations[0], correlations[4]));
			out.println(calcCorrelationCoeff(correlations[0], correlations[5]));
			out.println(calcCorrelationCoeff(correlations[0], correlations[6]));
			out.println(calcCorrelationCoeff(correlations[0], correlations[7]));
			// ----------------

			wout = new PrintStream(new FileOutputStream("refset"));
			wout.println("GENENAME\tADDF\tCAIEXT\tCAI\tGENELENGTH\tGC\tGC3");
			int nums[] = SortCais(correlations[0]);
			int n = nums.length / 100;
			for (int i = 0; i < n; i++) {
				int k = nums[i];
				wout.println(infos[0][k] + "\t" + infos[1][k] + "\t"
						+ correlations[7][k] + "\t" + correlations[0][k] + "\t"
						+ correlations[6][k] + "\t" + correlations[2][k] + "\t"
						+ correlations[3][k]);
			}
			wout.close();

			wout = new PrintStream(new FileOutputStream("cais.lst"));
			wout.println("GENENAME\tCAI\tADDF");
			for (int i = 0; i < nums.length; i++) {
				int k = nums[i];
				wout.println(infos[0][k] + "\t" + correlations[0][k] + "\t"
						+ infos[1][k]);
			}
			wout.close();

			// }

		} catch (Throwable t) {
			t.printStackTrace();
			System.exit(1);
		}
	}

	public static void GenerateCAIWValuesTable(float WW[], String CodeName) {
		TranslationTable TranTable = RNATools.getGeneticCode(CodeName);
		Alphabet Codons = TranTable.getSourceAlphabet();

		HashMap Mp = new HashMap();

		try {
			for (int i = 0; i < 64; i++) {
				SymbolList sDNA = DNATools.createDNA(Utils.TripletName(i));
				// out.println(sDNA.seqString());
				Symbol sb = ((FiniteAlphabet) Codons).getSymbol((RNATools
						.transcribe(sDNA)).toList());
				// out.println(sb.toString()+" => "+TranTable.translate(sb));
				float f = -1;
				String ss = TranTable.translate(sb).getName();
				if (Mp.containsKey(TranTable.translate(sb)))
					f = ((Float) Mp.get(TranTable.translate(sb))).floatValue();
				if (WW[i] > f)
					Mp.put(TranTable.translate(sb), new Float(WW[i]));
			}
			for (int i = 0; i < 64; i++) {
				SymbolList sDNA = DNATools.createDNA(Utils.TripletName(i));
				Symbol sb = ((FiniteAlphabet) Codons).getSymbol((RNATools
						.transcribe(sDNA)).toList());
				Symbol ob = TranTable.translate(sb);
				String ss = ob.getName();
				Float V = (Float) Mp.get(ob);
				if (V.floatValue() != 0)
					WW[i] /= ((Float) Mp.get(TranTable.translate(sb)))
							.floatValue();
				else
					WW[i] = 0;
			}

		} catch (Exception e) {
			out.println(e);
		}
	}

	public static void GenerateCAIWValuesTableRegardingGeneLength(float WW[],
			String CodeName, Vector V) {
		try {
			float WWW[][] = new float[V.size()][64];

			for (int i = 0; i < V.size(); i++) {
				StrandedFeature f = (StrandedFeature) (V.elementAt(i));
				String gs = f.getSymbols().seqString();
				gs = cutStartAndStop(gs, "UNIVERSAL");
				int Tc[] = Utils.calcTripletFreq(gs, 1, gs.length());
				float Tcf[] = new float[64];
				for (int j = 0; j < 64; j++)
					Tcf[j] = Tc[j] / (gs.length() / 3f);
				GenerateCAIWValuesTable(Tcf, CodeName);
				for (int j = 0; j < 64; j++)
					WWW[i][j] = Tcf[j];
			}
			for (int j = 0; j < 64; j++)
				WW[j] = 0;
			for (int i = 0; i < V.size(); i++)
				for (int j = 0; j < 64; j++)
					if (WWW[i][j] != 0)
						WW[j] += Math.log(WWW[i][j]);
			for (int j = 0; j < 64; j++)
				WW[j] /= V.size();
			for (int j = 0; j < 64; j++)
				WW[j] = (float) Math.exp(WW[j]);
		} catch (Exception e) {
			out.println(e);
		}
	}

	public static void GenerateCAIWValuesSuppressingRareCodons(float WW[],
			String CodeName, Vector V) {
		int occ[] = new int[64];
		try {
			for (int i = 0; i < V.size(); i++) {
				StrandedFeature f = null;
				String gs = null;
				if (!fastaFormat) {
					f = (StrandedFeature) (V.elementAt(i));
					gs = f.getSymbols().seqString();
				} else
					gs = ((Sequence) (V.elementAt(i))).seqString();
				gs = cutStartAndStop(gs, "UNIVERSAL");
				int Tc[] = Utils.calcTripletFreq(gs, 1, gs.length());
				for (int j = 0; j < 64; j++)
					if (Tc[j] != 0)
						occ[j]++;
			}
			GenerateCAIWValuesTable(WW, CodeName);
			for (int j = 0; j < 64; j++) {
				// WW[j]/=V.size()-occ[j]+1;
				WW[j] *= (float) occ[j] / (float) V.size();
				GTFAverLen[j] = WW[j];
			}
		} catch (Exception e) {
			out.println(e);
		}
	}

	public static void GenerateCAIWValuesTableMethodAverage(float WW[],
			String CodeName) {
		TranslationTable TranTable = RNATools.getGeneticCode(CodeName);
		Alphabet Codons = TranTable.getSourceAlphabet();

		HashMap Mp = new HashMap();
		HashMap MpN = new HashMap();

		try {
			for (int i = 0; i < 64; i++) {
				SymbolList sDNA = DNATools.createDNA(Utils.TripletName(i));
				// out.println(sDNA.seqString());
				Symbol sb = ((FiniteAlphabet) Codons).getSymbol((RNATools
						.transcribe(sDNA)).toList());
				// out.println(sb.toString()+" => "+TranTable.translate(sb));
				float f = 0;
				int nn = 0;
				if (Mp.containsKey(TranTable.translate(sb)))
					f = ((Float) Mp.get(TranTable.translate(sb))).floatValue();
				if (MpN.containsKey(TranTable.translate(sb)))
					nn = ((Integer) MpN.get(TranTable.translate(sb)))
							.intValue();
				Mp.put(TranTable.translate(sb), new Float(WW[i] + f));
				MpN.put(TranTable.translate(sb), new Integer(nn + 1));
				// out.println(sb.getName()+" "+(WW[i]+f)+" "+(nn+1));
			}
			for (int i = 0; i < 64; i++) {
				SymbolList sDNA = DNATools.createDNA(Utils.TripletName(i));
				Symbol sb = ((FiniteAlphabet) Codons).getSymbol((RNATools
						.transcribe(sDNA)).toList());
				if (((Float) Mp.get(TranTable.translate(sb))).floatValue() != 0)
					WW[i] /= ((Float) Mp.get(TranTable.translate(sb)))
							.floatValue();
				else
					WW[i] = 0;
			}

		} catch (Exception e) {
			out.println(e);
		}
	}

	public static int[] calcAminoacidFreq(String s, Vector Acids, int st, int en) {
		int[] Res = new int[20];

		try {

			String ss = s.substring(st - 1, en - 1);
			int len = ss.length();
			SymbolList Tripl;
			for (int i = 0; i < len / 3; i++) {
				String sTr = ss.substring(i * 3, i * 3 + 3);
				int n = 0;
				// out.println(Amin.seqString()+" "+sTr);
				int k = Utils.TripletNumber(sTr);
				k = ((Integer) TranslationalTable.elementAt(k)).intValue();
				if (k != -1) {
					Res[k] += 1;
				}
			}
		} catch (Exception e) {
			out.println(e);
		}
		;
		return Res;
	}

	public static int[] orderCodonUsageByAminoacid(Vector Acids) {
		int Res[] = new int[64];
		try {
			int k = 0;

			for (int j = 0; j < Acids.size(); j++) {

				String acName = new String((String) Acids.elementAt(j));

				for (int i = 0; i < 64; i++) {
					SymbolList sDNA = DNATools.createDNA(Utils.TripletName(i,
							"a", "t", "g", "c"));
					// System.out.println(sDNA.seqString());
					String trAcid = (RNATools.translate(RNATools
							.transcribe(sDNA))).seqString();
					if (trAcid.equals(acName)) {
						Res[k] = i;
						k++;
					}
				}
			}
		} catch (Exception e) {
			out.println(e + " in  orderCodonUsageByAminoacid");
		}
		return Res;
	}

	public static float CalculateCAIValue(String s, float CAIWValuesTable[])
			throws Exception {
		String ss = cutStartAndStop(s, "UINIVERSAL");
		int len = ss.length();
		float res = 0;
		int L = 0;
		for (int i = 0; i < len / 3; i++) {
			String sTr = ss.substring(i * 3, i * 3 + 3);
			int n = Utils.TripletNumber(sTr);
			if (CAIWValuesTable[n] != 0) {
				res += Math.log(CAIWValuesTable[n]);
				L++;
			} else {
				res += Math.log(0.01);
				L++;
			}
		}
		return (float) Math.exp((double) (1.0 / (float) L * res));
	}

	public static float CalculateCAIValue(float cu[], float CAIWValuesTable[]) {
		double res = 0;
		int L = 0;
		for (int i = 0; i < 64; i++) {
			if (CAIWValuesTable[i] != 0) {
				res += cu[i] * Math.log(CAIWValuesTable[i]);
			} else {
				res += cu[i] * Math.log(0.01);
			}
		}
		return (float) Math.exp((double) (res));
	}

	public static void MakeCAIIteration(Vector ListOfGenes, float cais[],
			int prop, float CodUs[]) {
		try {

			float memcais[] = new float[cais.length];
			for (int i = 0; i < cais.length; i++) {
				memcais[i] = cais[i];
			}

			int geneNumber = (int) (0.01 * prop * ListOfGenes.size());
			if (geneNumber < minimumRefsetSize) {
				geneNumber = minimumRefsetSize;
				prop = (int) (100f * geneNumber / ListOfGenes.size());
			}
			out.println("First " + geneNumber + " genes (" + prop + "%)");
			int nums[] = SortCais(cais);
			// for (int i = 0; i < nums.length; i++)
			// out.print(cais[nums[i]]+" ");
			// out.println();
			Vector BiggestGenes = new Vector();
			for (int i = 0; i < geneNumber; i++)
				BiggestGenes.addElement(ListOfGenes.elementAt(nums[i]));

			// Calculate new codon usage

			int TF[];
			float GTF[] = new float[64];

			for (int i = 0; i < BiggestGenes.size(); i++) {
				StrandedFeature Ft = null;
				String ss = null;
				if (!fastaFormat) {
					Ft = (StrandedFeature) BiggestGenes.elementAt(i);
					ss = new String(Ft.getSymbols().seqString());
				} else
					ss = ((Sequence) (BiggestGenes.elementAt(i))).seqString();
				String sss = cutStartAndStop(ss, "UNIVERSAL");
				TF = Utils.calcTripletFreq(sss, 1, sss.length());
				for (int j = 0; j < 64; j++)
					GTF[j] += TF[j];
			}

			int sum = 0;
			for (int i = 0; i < 64; i++)
				sum += GTF[i];
			for (int i = 0; i < 64; i++)
				GTF[i] = GTF[i] / sum;
			for (int i = 0; i < 64; i++)
				CodUs[i] = GTF[i];

			// Calculate new cai table

			Utils.printDistribution(GTF);
			if (taskModifier == 0)
				GenerateCAIWValuesTable(GTF, "UNIVERSAL");
			if (taskModifier == 1)
				GenerateCAIWValuesTableMethodAverage(GTF, "UNIVERSAL");
			if (taskModifier == 2)
				GenerateCAIWValuesTableRegardingGeneLength(GTF, "UNIVERSAL",
						BiggestGenes);
			if (taskModifier == 3)
				GenerateCAIWValuesSuppressingRareCodons(GTF, "UNIVERSAL",
						BiggestGenes);

			// Calculate new cais

			for (int i = 0; i < ListOfGenes.size(); i++) {
				StrandedFeature Ft = null;
				String ss = null;
				if (!fastaFormat) {
					Ft = (StrandedFeature) ListOfGenes.elementAt(i);
					ss = new String(Ft.getSymbols().seqString());
				} else
					ss = ((Sequence) ListOfGenes.elementAt(i)).seqString();
				float caiv = CalculateCAIValue(ss, GTF);
				cais[i] = caiv;
			}

		} catch (Exception e) {
			out.print("\n" + e);
		}

	}

	public static int[] SortCais(float cais[]) {
		int res[] = new int[cais.length];
		for (int i = 0; i < res.length; i++)
			res[i] = i;

		int i, j, k, inc, n = cais.length;
		float v;

		inc = 1;
		do {
			inc *= 3;
			inc++;
		} while (inc <= n);

		do {
			inc /= 3;
			for (i = inc + 1; i <= n; i++) {
				v = cais[res[i - 1]];
				j = i;
				k = res[i - 1];
				while (cais[res[j - inc - 1]] < v) {
					// cais[j]=cais[j-inc];
					res[j - 1] = res[j - inc - 1];
					j -= inc;
					if (j <= inc)
						break;
				}
				// cais[j]=v;
				res[j - 1] = k;
			}
		} while (inc > 0);

		return res;
	}

	public static String getGenomeName(String geneID) {
		String res = "";
		Set keys = genomes.keySet();
		for (Iterator i = keys.iterator(); i.hasNext();) {
			String gn = (String) i.next();
			Vector v = (Vector) (genomes.get(gn));
			if (v.indexOf(geneID) >= 0) {
				res = gn;
				break;
			}
		}
		return res;
	}

	public static void recalculateCodonsForRSCU(float WW[], String CodeName) {
		TranslationTable TranTable = RNATools.getGeneticCode(CodeName);
		Alphabet Codons = TranTable.getSourceAlphabet();

		HashMap Mp = new HashMap();
		HashMap MpN = new HashMap();

		try {
			for (int i = 0; i < 64; i++) {
				SymbolList sDNA = DNATools.createDNA(Utils.TripletName(i));
				// out.println(sDNA.seqString());
				Symbol sb = ((FiniteAlphabet) Codons).getSymbol((RNATools
						.transcribe(sDNA)).toList());
				// out.println(sb.toString()+" => "+TranTable.translate(sb));
				float f = 0;
				int nn = 0;
				if (Mp.containsKey(TranTable.translate(sb)))
					f = ((Float) Mp.get(TranTable.translate(sb))).floatValue();
				if (MpN.containsKey(TranTable.translate(sb)))
					nn = ((Integer) MpN.get(TranTable.translate(sb)))
							.intValue();
				Mp.put(TranTable.translate(sb), new Float(WW[i] + f));
				MpN.put(TranTable.translate(sb), new Integer(nn + 1));
				// out.println(sb.getName()+" "+(WW[i]+f)+" "+(nn+1));
			}
			for (int i = 0; i < 64; i++) {
				SymbolList sDNA = DNATools.createDNA(Utils.TripletName(i));
				Symbol sb = ((FiniteAlphabet) Codons).getSymbol((RNATools
						.transcribe(sDNA)).toList());
				if (((Float) Mp.get(TranTable.translate(sb))).floatValue() != 0)
					WW[i] /= ((Float) Mp.get(TranTable.translate(sb)))
							.floatValue()
							/ ((Integer) MpN.get(TranTable.translate(sb)))
									.intValue();
				else
					WW[i] = 0;
			}

		} catch (Exception e) {
			out.println(e);
		}

	}

	public static void calculateCodonUsage(Vector ListOfGenes, float[] CodUs) {
		try {
			int TF[];
			float GTF[] = new float[64];
			for (int i = 0; i < ListOfGenes.size(); i++) {
				String ss = getSeqString(ListOfGenes.elementAt(i));
				String sss = cutStartAndStop(ss, "UNIVERSAL");
				TF = Utils.calcTripletFreq(sss, 1, sss.length());
				for (int j = 0; j < 64; j++)
					GTF[j] += TF[j];
			}
			int sum = 0;
			for (int i = 0; i < 64; i++)
				sum += GTF[i];
			for (int i = 0; i < 64; i++)
				CodUs[i] = GTF[i] / sum;
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static float[] calculateCodonUsage(Object sq) {
		float res[] = new float[64];
		try {
			String ss = null;
			if (sq instanceof String)
				ss = (String) sq;
			else
				ss = getSeqString(sq);
			String sss = cutStartAndStop(ss, "UNIVERSAL");
			int TF[] = Utils.calcTripletFreq(sss, 1, sss.length());
			for (int j = 0; j < 64; j++)
				res[j] += (float) TF[j] / (sss.length() / 3.0f);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return res;
	}

	public static String getSeqString(Object sq) {
		StrandedFeature Ft = null;
		String ss = null;
		if (!fastaFormat) {
			Ft = (StrandedFeature) sq;
			ss = new String(Ft.getSymbols().seqString());
		} else
			ss = ((Sequence) sq).seqString();
		return ss;
	}

	public static String getGeneName(Object sq) {
		Annotatable Ft = (Annotatable) sq;
		String gName = "Unknown";
		if (Ft.getAnnotation().containsProperty("gene"))
			gName = Ft.getAnnotation().getProperty("gene").toString();
		else if (Ft.getAnnotation().containsProperty("standard_name"))
			gName = Ft.getAnnotation().getProperty("standard_name").toString();
		else if (Ft.getAnnotation().containsProperty(idFeature))
			gName = Ft.getAnnotation().getProperty(idFeature).toString();
		else if (((Sequence)Ft).getName().length() > 0){
			gName = ((Sequence)Ft).getName();
		}
		gName = gName.replace(' ', '_');
		return gName;
	}

	public static String getGeneID(Object sq) {
		String gn = getGeneName(sq);
		String s = getSeqString(sq);
		return gn + "_" + s.length();
	}

	public static String getGeneDescription(Object sq) {
		String res = getGeneName(sq);
		Annotatable Ft = (Annotatable) sq;
		String af = new String("");
		if (addFeature != null) {
			if (Ft.getAnnotation().containsProperty("gene"))
				if (Ft.getAnnotation().containsProperty(addFeature))
					af = new String(Ft.getAnnotation().getProperty(addFeature)
							.toString());
			af = af.replace(' ', '_');
			if (af.length() > 100)
				af = af.substring(0, 100);
		}
		return res + "\t" + af;
	}

	public static String cutStartAndStop(String s, String CodeName)
			throws Exception {
		String res = s;
		// res = s.substring(3,s.length()-3);
		TranslationTable TranTable = RNATools.getGeneticCode(CodeName);

		SymbolList sDNA = DNATools.createDNA(s.substring(0, 3));
		SymbolList Amin = RNATools.translate(RNATools.transcribe(sDNA));
		String sss = new String(Amin.seqString());
		if (sss.equals("M"))
			res = res.substring(3, res.length());

		sDNA = DNATools.createDNA(s.substring(s.length() - 3, s.length()));
		Amin = RNATools.translate(RNATools.transcribe(sDNA));
		sss = new String(Amin.seqString());
		if (sss.equals("*"))
			res = res.substring(0, res.length() - 3);

		return res;
	}

	public static float calcCorrelationCoeff(float m1[], float m2[]) {
		float res = 0;
		int N = m1.length;
		float xy = 0f, x2 = 0f, y2 = 0f, x = 0f, y = 0f;
		for (int i = 0; i < N; i++) {
			xy += m1[i] * m2[i];
			x += m1[i];
			y += m2[i];
			x2 += m1[i] * m1[i];
			y2 += m2[i] * m2[i];
		}
		res = (float) ((xy - x * y / N) / Math.sqrt((x2 - x * x / N)
				* (y2 - y * y / N)));
		return res;
	}

	public static float calcTValue(float m1[], float m2[]) {
		float res = 0;
		int N1 = 0;
		int N2 = 0;
		float xy = 0f, x2 = 0f, y2 = 0f, x = 0f, y = 0f;
		for (int i = 0; i < m1.length; i++) {
			if (m2[i] > 0) {
				x += m1[i];
				x2 += m1[i] * m1[i];
				N1++;
			} else {
				y += m1[i];
				y2 += m1[i] * m1[i];
				N2++;
			}
		}
		float disp1 = x2 / N1 - x * x / N1 / N1;
		float disp2 = y2 / N2 - y * y / N2 / N2;
		// System.out.println("Nx = "+N1+"Ny = "+N2+" x = "+x/N1+" y = "+y/N2+" dx = "+disp1+" dy = "+disp2);
		res = (float) ((x / N1 - y / N2) / Math.sqrt(disp1 + disp2));
		return res;
	}

	public static String calcCodonBiasCode(float GeneCodonUsage[],
			String codeName) {
		StringBuffer res = new StringBuffer("");

		TranslationTable TranTable = RNATools.getGeneticCode(codeName);
		Alphabet Codons = TranTable.getSourceAlphabet();

		try {

			Vector amAcids = new Vector();
			Vector amlist = new Vector();

			for (int i = 0; i < 64; i++) {
				SymbolList sDNA = DNATools.createDNA(Utils.TripletName(i));
				Symbol sb = ((FiniteAlphabet) Codons).getSymbol((RNATools
						.transcribe(sDNA)).toList());
				String amin = TranTable.translate(sb).toString();
				amlist.add(amin);
				if (amAcids.indexOf(amin) == -1)
					amAcids.add(amin);
			}

			Vector ls = new Vector();
			for (int i = 0; i < amAcids.size(); i++) {
				String am = (String) amAcids.elementAt(i);
				ls.clear();
				for (int j = 0; j < 64; j++) {
					String s = (String) amlist.elementAt(j);
					if (s.equals(am))
						ls.add(new Float(GeneCodonUsage[j]));
				}
				int imax = -1;
				float mx = -1f;
				for (int j = 0; j < ls.size(); j++) {
					float f = ((Float) ls.elementAt(j)).floatValue();
					if (f > mx) {
						mx = f;
						imax = j;
					}
				}
				res.append((imax + 1));
			}
		} catch (Exception e) {
			e.printStackTrace();
		}

		return res.toString();
	}

}
