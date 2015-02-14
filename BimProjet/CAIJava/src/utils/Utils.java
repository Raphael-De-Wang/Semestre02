package utils;

import java.io.*;
import java.util.*;
import java.text.*;

import org.biojava.bio.*;
import org.biojava.bio.symbol.*;
import org.biojava.bio.seq.*;
import org.biojava.bio.seq.io.*;

public class Utils {

	public static PrintStream out = System.out;
	public static Vector centroids;
	public static Vector deviations;
	public static Vector sizes;

	public Utils() {
	}

	public static int[] calcTripletFreq(SymbolList seq, int st, int en)
			throws Exception {
		int[] Res = new int[64];
		SymbolList LS = seq.subList(st, en);
		int len = LS.length();
		SymbolList Tripl;
		for (int i = 0; i < len / 3; i++) {
			Tripl = LS.subList(i * 3 + 1, i * 3 + 3);
			String sTr = Tripl.seqString();
			int n = TripletNumber(sTr.toString());
			Res[n] += 1;
		}
		return Res;
	}

	public static int[] calcBaseFreq(String seq, int st, int en)
			throws Exception {
		int[] Res = new int[4];
		String LS = seq.substring(st, en);
		int len = LS.length();
		SymbolList Tripl;
		for (int i = 0; i < len; i++) {
			int n = BaseNumber(LS.substring(i, i + 1));
			Res[n] += 1;
		}
		return Res;
	}

	public static float[] calcBaseFreqF(String seq, int st, int en)
			throws Exception {
		float[] Res = new float[4];
		int r[] = calcBaseFreq(seq, st, en);
		for (int i = 0; i < 4; i++)
			Res[i] = (float) r[i] / (en - st + 1);
		return Res;
	}

	public static int[] calcDipletFreqOverlapping(String seq, int st, int en)
			throws Exception {
		int[] Res = new int[16];
		String LS = seq.substring(st, en);
		int len = LS.length();
		for (int i = 0; i < len - 1; i++) {
			String dipl = LS.substring(i, i + 2);
			int n = DipletNumber(dipl);
			Res[n] += 1;
		}
		return Res;
	}

	public static int[] calcDipletFreq(String seq, int st, int en)
			throws Exception {
		int[] Res = new int[16];
		String LS = seq.substring(st, en);
		int len = LS.length();
		int k = 0;
		for (int i = 0; i < len / 2; i++) {
			String dipl = LS.substring(k, k + 2);
			int n = DipletNumber(dipl);
			Res[n] += 1;
			k += 2;
		}
		return Res;
	}

	public static int[] calcQuadroFreq(String seq, int st, int en)
			throws Exception {
		int[] Res = new int[256];
		String LS = seq.substring(st, en);
		int len = LS.length();
		int k = 0;
		for (int i = 0; i < len / 4 - 1; i++) {
			String dipl = LS.substring(k, k + 4);
			int n = QuadroNumber(dipl);
			Res[n] += 1;
			k += 4;
		}
		return Res;
	}

	public static float[] calcDipletFreqFOverlapping(String seq, int st, int en)
			throws Exception {
		float[] Res = new float[16];
		int r[] = calcDipletFreqOverlapping(seq, st, en);
		for (int i = 0; i < 16; i++)
			Res[i] = (float) r[i] / (seq.length());
		return Res;
	}

	public static float[] calcDipletFreqF(String seq, int st, int en)
			throws Exception {
		float[] Res = new float[16];
		int r[] = calcDipletFreq(seq, st, en);
		for (int i = 0; i < 16; i++)
			Res[i] = (float) r[i] / ((en - st + 1) / 2);
		return Res;
	}

	public static float[] calcQuadroFreqF(String seq, int st, int en)
			throws Exception {
		float[] Res = new float[256];
		int r[] = calcQuadroFreq(seq, st, en);
		for (int i = 0; i < 256; i++)
			Res[i] = (float) r[i] / ((en - st + 1) / 4);
		return Res;
	}

	public static int[] calcTripletFreq(String s, int st, int en)
			throws Exception {
		int[] Res = new int[64];
		// SymbolList LS = seq.subList(st,en);
		String ss = s.substring(st - 1, en);
		int len = ss.length();
		SymbolList Tripl;
		// System.out.println(ss+" "+len+" "+(int)(len/3f+0.1f));
		for (int i = 0; i < (int) (len / 3f + 0.1f); i++) {
			// Tripl = LS.subList(i*3+1,i*3+3);
			String sTr = "";
			sTr = ss.substring(i * 3, i * 3 + 3);
			// System.out.println(i+" "+sTr);
			int n = TripletNumber(sTr);
			Res[n] += 1;
		}
		return Res;
	}

	public static int[] calcTripletPairsFreq(String s, int st, int en)
			throws Exception {
		int res[] = new int[4096];
		String ss = s.substring(st - 1, en);
		int len = ss.length();
		SymbolList Tripl;
		for (int i = 0; i < (int) (len / 3f + 0.1f) - 1; i++) {
			// Tripl = LS.subList(i*3+1,i*3+3);
			String sTr1 = "", sTr2 = "";
			sTr1 = ss.substring(i * 3, i * 3 + 3);
			sTr2 = ss.substring((i + 1) * 3, (i + 1) * 3 + 3);
			// System.out.println(i+" "+sTr);
			int n = TripletPairNumber(sTr1, sTr2);
			res[n] += 1;
		}
		return res;
	}

	public static float[] calcTripletFreqF(SymbolList seq, int pos, int wSize,
			int Bayez) throws Exception {
		float[] Res = new float[64];
		int nt = (int) (wSize / 6);
		int ntb = 0;
		if (Bayez > 1)
			ntb = (int) (wSize / 6 / Bayez);
		nt = (int) (nt / 2) * 2 + 1;
		ntb = (int) (ntb / 2) * 2 + 1;
		int st = pos - 1 - (nt - 1) * 3, en = pos + 1 + (nt - 1) * 3;
		int stb = pos - 1 - (ntb - 1) * 3, enb = pos + 1 + (ntb - 1) * 3;
		int fi[] = calcTripletFreq(seq, st, en);
		int fib[] = new int[64];
		if (Bayez > 1)
			fib = calcTripletFreq(seq, stb, enb);
		for (int i = 0; i < 64; i++) {
			if (Bayez > 1)
				Res[i] = (((float) fi[i] + (float) fib[i]) / ((nt + ntb) * 2));
			else
				Res[i] = (float) fi[i] / (nt * 2);
		}
		return Res;
	}

	public static float[] calcTripletFreqF(String s, int pos, int wSize,
			int Bayez) throws Exception {
		float[] Res = new float[64];
		int nt = (int) (wSize / 3);
		int ntb = 0;
		if (Bayez > 1)
			ntb = (int) (wSize / 3 / Bayez);
		nt = (int) (nt / 2) * 2 + 1;
		ntb = (int) (ntb / 2) * 2 + 1;
		int st = pos - 1 - (nt - 1) * 3, en = pos + 1 + (nt - 1) * 3;
		int stb = pos - 1 - (ntb - 1) * 3, enb = pos + 1 + (ntb - 1) * 3;
		int fi[] = calcTripletFreq(s, st, en);
		int fib[] = new int[64];
		if (Bayez > 1)
			fib = calcTripletFreq(s, stb, enb);
		for (int i = 0; i < 64; i++) {
			if (Bayez > 1)
				Res[i] = (((float) fi[i] + (float) fib[i]) / ((nt + ntb) * 2));
			else
				Res[i] = (float) fi[i] / (nt * 2);
		}
		return Res;
	}

	/*
	 * public static float[] calcTripletFreqF(String s, int st, int en) throws
	 * Exception { int[] Res = new int[64]; //SymbolList LS =
	 * seq.subList(st,en); String ss = s.substring(st-1,en-1); int len =
	 * ss.length(); SymbolList Tripl; for(int i=0; i<len/3; i++) { //Tripl =
	 * LS.subList(i*3+1,i*3+3); String sTr = ss.substring(i*3,i*3+3); int n =
	 * TripletNumber(sTr.toString()); Res[n]+=1; } return Res; }
	 */

	public static int[] calcTripletFreqComplementary(SymbolList seq, int st,
			int en) throws Exception {
		int[] Res = new int[64];
		SymbolList LS = seq.subList(st, en);
		int len = LS.length();
		SymbolList Tripl;
		for (int i = 0; i < len / 3; i++) {
			Tripl = DNATools.complement(LS.subList(i * 3 + 1, i * 3 + 3));
			StringBuffer sTr = new StringBuffer(Tripl.seqString());
			int n = TripletNumber(sTr.reverse().toString());
			Res[n] += 1;
		}
		return Res;
	}

	public static int TripletNumber(String st) {
		char[] cA = new char[3];
		st.getChars(0, 3, cA, 0);
		return BaseNumber(new String(cA, 0, 1)) * 16
				+ BaseNumber(new String(cA, 1, 1)) * 4
				+ BaseNumber(new String(cA, 2, 1));
	}

	public static int TripletPairNumber(String st1, String st2) {
		char[] cA = new char[6];
		st1.getChars(0, 3, cA, 0);
		st2.getChars(0, 3, cA, 3);
		return BaseNumber(new String(cA, 0, 1)) * 1024
				+ BaseNumber(new String(cA, 1, 1)) * 256
				+ BaseNumber(new String(cA, 2, 1)) * 64
				+ BaseNumber(new String(cA, 3, 1)) * 16
				+ BaseNumber(new String(cA, 4, 1)) * 4
				+ BaseNumber(new String(cA, 5, 1));
	}

	public static int DipletNumber(String st) {
		char[] cA = new char[2];
		st.getChars(0, 2, cA, 0);
		return BaseNumber(new String(cA, 0, 1)) * 4
				+ BaseNumber(new String(cA, 1, 1));
	}

	public static int TripletNumber(String st, String let1, String let2,
			String let3, String let4) {
		char[] cA = new char[3];
		st.getChars(0, 3, cA, 0);
		return BaseNumber(new String(cA, 0, 1), let1, let2, let3, let4) * 16
				+ BaseNumber(new String(cA, 1, 1), let1, let2, let3, let4) * 4
				+ BaseNumber(new String(cA, 2, 1), let1, let2, let3, let4);
	}

	public static String BaseName(int n) {
		String s = new String("");
		switch (n) {
		case 0:
			s = new String("a");
			break;
		case 1:
			s = new String("c");
			break;
		case 2:
			s = new String("g");
			break;
		case 3:
			s = new String("t");
			break;
		}
		return s;
	}

	public static String BaseName(int n, String let1, String let2, String let3,
			String let4) {
		String s = new String("");
		switch (n) {
		case 0:
			s = new String(let1);
			break;
		case 1:
			s = new String(let2);
			break;
		case 2:
			s = new String(let3);
			break;
		case 3:
			s = new String(let4);
			break;
		}
		return s;
	}

	public static String TripletName(int n) {
		StringBuffer s = new StringBuffer();
		int i1 = (int) ((float) n / 16);
		int i2 = (int) ((float) (n - i1 * 16) / 4);
		int i3 = n - i1 * 16 - i2 * 4;
		s.append(BaseName(i1));
		s.append(BaseName(i2));
		s.append(BaseName(i3));
		return s.toString();
	}

	public static String QuadroName(int n) {
		StringBuffer s = new StringBuffer();
		int i1 = (int) ((float) n / 64);
		int i2 = (int) ((float) (n - i1 * 64) / 16);
		int i3 = (int) ((float) (n - i1 * 64 - i2 * 16) / 4);
		int i4 = n - i1 * 64 - i2 * 16 - i3 * 4;
		s.append(BaseName(i1));
		s.append(BaseName(i2));
		s.append(BaseName(i3));
		s.append(BaseName(i4));
		return s.toString();
	}

	public static int QuadroNumber(String st) {
		char[] cA = new char[4];
		st.getChars(0, 4, cA, 0);
		return BaseNumber(new String(cA, 0, 1)) * 64
				+ BaseNumber(new String(cA, 1, 1)) * 16
				+ BaseNumber(new String(cA, 2, 1)) * 4
				+ BaseNumber(new String(cA, 3, 1));
	}

	public static String TripletPairName(int n) {
		StringBuffer s = new StringBuffer();
		int i1 = (int) ((float) n / 1024);
		int i2 = (int) ((float) (n - i1 * 1024) / 256);
		int i3 = (int) ((float) (n - i1 * 1024 - i2 * 256) / 64);
		int i4 = (int) ((float) (n - i1 * 1024 - i2 * 256 - i3 * 64) / 16);
		int i5 = (int) ((float) (n - i1 * 1024 - i2 * 256 - i3 * 64 - i4 * 16) / 4);
		int i6 = n - i1 * 1024 - i2 * 256 - i3 * 64 - i4 * 16 - i5 * 4;
		s.append(BaseName(i1));
		s.append(BaseName(i2));
		s.append(BaseName(i3));
		s.append(BaseName(i4));
		s.append(BaseName(i5));
		s.append(BaseName(i6));
		return s.toString();
	}

	public static String DipletName(int n) {
		StringBuffer s = new StringBuffer();
		int i1 = (int) ((float) n / 4);
		int i2 = n - i1 * 4;
		s.append(BaseName(i1));
		s.append(BaseName(i2));
		return s.toString();
	}

	public static String TripletName(int n, String let1, String let2,
			String let3, String let4) {
		StringBuffer s = new StringBuffer();
		int i1 = (int) ((float) n / 16);
		int i2 = (int) ((float) (n - i1 * 16) / 4);
		int i3 = n - i1 * 16 - i2 * 4;
		s.append(BaseName(i1, let1, let2, let3, let4));
		s.append(BaseName(i2, let1, let2, let3, let4));
		s.append(BaseName(i3, let1, let2, let3, let4));
		return s.toString();
	}

	public static int BaseNumber(String cn) {
		int Res = 0;
		if (cn.equals("a"))
			Res = 0;
		if (cn.equals("c"))
			Res = 1;
		if (cn.equals("g"))
			Res = 2;
		if (cn.equals("t"))
			Res = 3;
		return Res;
	}

	public static int BaseNumber(String cn, String let1, String let2,
			String let3, String let4) {
		int Res = 0;
		if (cn.equals(let1))
			Res = 0;
		if (cn.equals(let2))
			Res = 1;
		if (cn.equals(let3))
			Res = 2;
		if (cn.equals(let4))
			Res = 3;
		return Res;
	}

	public static float GCContent(String cn) {
		float gc = 0;
		for (int i = 0; i < cn.length(); i++) {
			if ((cn.charAt(i) == 'g') || (cn.charAt(i) == 'c'))
				gc += 1;
		}
		return gc / cn.length();
	}

	public static float GCSkew(String cn) {
		float g = 0;
		float c = 0, res = 0;
		for (int i = 0; i < cn.length(); i++) {
			if (cn.charAt(i) == 'g')
				g += 1f;
			if (cn.charAt(i) == 'c')
				c += 1f;
		}
		if ((g + c) == 0)
			res = 0;
		else
			res = (g - c) / (g + c);
		return res;
	}

	public static float ATSkew(String cn) {
		float a = 0;
		float t = 0, res = 0;
		for (int i = 0; i < cn.length(); i++) {
			if (cn.charAt(i) == 'a')
				a += 1f;
			if (cn.charAt(i) == 't')
				t += 1f;
		}
		if ((a + t) == 0)
			res = 0;
		else
			res = (a - t) / (a + t);
		return res;
	}

	public static float GCContent3(String cn) {
		float gc = 0;
		for (int i = 2; i < cn.length(); i += 3) {
			if ((cn.charAt(i) == 'g')
					&& (!((cn.charAt(i - 2) == 't') && (cn.charAt(i - 1) == 'g')))) {
				gc += 1;
				continue;
			}
			if ((cn.charAt(i) == 'g')
					&& (!((cn.charAt(i - 2) == 'a') && (cn.charAt(i - 1) == 't')))) {
				gc += 1;
				continue;
			}
			// if(cn.charAt(i)=='g') gc+=1;
			if (cn.charAt(i) == 'c')
				gc += 1;
		}
		float res = gc / (cn.length() / 3.0f);
		return res;
	}

	public static int CalcIsCoding(Sequence seq, int pos) {
		int Res = 0;
		/*
		 * for(Iterator i=seq.features(); (i.hasNext())&&(Res==0); ) { Feature
		 * Ft = (Feature) i.next(); if
		 * ((Ft.getType().equals("CDS"))&&(Ft.getLocation().contains(pos)))
		 * Res=1; }
		 */
		String s = GetFeatureName(seq, pos);
		if ((s.indexOf("CDS") > -1) && (s.indexOf("intron") == -1))
			Res = 1;
		return Res;
	}

	public static int CalcIsCodingGene(Sequence seq, int pos) {
		int Res = 0;
		/*
		 * for(Iterator i=seq.features(); (i.hasNext())&&(Res==0); ) { Feature
		 * Ft = (Feature) i.next(); if
		 * ((Ft.getType().equals("CDS"))&&(Ft.getLocation().contains(pos)))
		 * Res=1; }
		 */
		String s = GetFeatureName(seq, pos);
		if (s.indexOf("gene") > -1)
			Res = 1;
		return Res;
	}

	public static int CalcIsCoding2(Sequence seq, int pos) {
		int Res = 0;
		/*
		 * for(Iterator i=seq.features(); (i.hasNext())&&(Res==0); ) { Feature
		 * Ft = (Feature) i.next(); if
		 * ((Ft.getType().equals("CDS"))&&(Ft.getLocation().contains(pos)))
		 * Res=1; }
		 */
		String s = GetFeatureName(seq, pos);
		if (((s.indexOf("CDS") > -1) || (s.indexOf("tRNA") > -1) || (s
				.indexOf("rRNA") > -1)) && (s.indexOf("intron") == -1))
			Res = 1;
		return Res;
	}

	public static int CalcIsGeneComplement(Sequence seq, int pos) {
		/*
		 * int Res = 0; for(Iterator i=seq.features(); (i.hasNext())&&(Res==0);
		 * ) { Feature Ft = (Feature) i.next(); if
		 * ((Ft.getType().equals("CDS"))&&(Ft.getLocation().contains(pos))) {
		 * Res=1; for (Iterator ai =
		 * Ft.getAnnotation().asMap().entrySet().iterator(); ai.hasNext(); ) {
		 * Map.Entry me = (Map.Entry) ai.next(); if
		 * (me.getKey().toString().equals("internal_data")) if
		 * (me.getValue().toString().indexOf("complement")>=0) Res=-1; } } }
		 * return Res;
		 */
		int Res = 0;
		for (Iterator i = seq.features(); (i.hasNext()) && (Res == 0);) {
			StrandedFeature Ft = (StrandedFeature) i.next();
			if ((Ft.getType().equals("CDS"))
					&& (Ft.getLocation().contains(pos))) {
				Res = 1;
				Res = Ft.getStrand().getValue();
			}
		}
		return Res;
	}

	public static int CalcIsGeneComplement(Feature Ft) {
		int Res = 1;
		for (Iterator ai = Ft.getAnnotation().asMap().entrySet().iterator(); ai
				.hasNext();) {
			Map.Entry me = (Map.Entry) ai.next();
			if (me.getKey().toString().equals("internal_data"))
				if (me.getValue().toString().indexOf("complement") >= 0)
					Res = -1;
		}
		return Res;
	}

	public static int getGenePosition(Feature Ft) {
		int Res = 0;
		Res = (int) (Ft.getLocation().getMin() + (Ft.getLocation().getMax() - Ft
				.getLocation().getMin()) / 2);
		return Res;
	}

	public static int NumberOfBase(StringBuffer seq, char b) {
		int c = 0;
		for (int i = 0; i < seq.length(); i++)
			if (seq.charAt(i) == b)
				c++;
		return c;
	}

	public static String GetFeatureName(Sequence seq, int pos) {
		int Res = 0;
		String sRes = "";
		for (Iterator i = seq.features(); (i.hasNext()) && (Res == 0);) {
			Feature Ft = (Feature) i.next();
			if (Ft.getLocation().contains(pos)) {
				// Res=1;
				if (!(Ft.getType().equals("source"))) {
					if (sRes.equals(""))
						sRes = new String(Ft.getType());
					else
						sRes = new String(sRes + "," + Ft.getType());
				}
			}
		}
		return sRes;
	}

	public static Vector GetFeature(Sequence seq, int pos) {
		Vector sRes = new Vector();
		for (Iterator i = seq.features(); i.hasNext();) {
			Feature Ft = (Feature) i.next();
			if (!(Ft.getType().equals("source")))
				if (Ft.getLocation().contains(pos)) {
					sRes.add(Ft);
				}
		}
		return sRes;
	}

	public static Feature GetFeatureCDS(Sequence seq, int pos) {
		Feature f = null;
		for (Iterator i = seq.features(); i.hasNext();) {
			Feature Ft = (Feature) i.next();
			if (Ft.getLocation().contains(pos))
				if (Ft.getType().equals("CDS")) {
					f = Ft;
					break;
				}
		}
		return f;
	}

	public static StrandedFeature GetStrandedFeatureCDS(Sequence seq, int pos) {
		StrandedFeature f = null;
		for (Iterator i = seq.features(); i.hasNext();) {
			StrandedFeature Ft = (StrandedFeature) i.next();
			if (Ft.getLocation().contains(pos))
				if (Ft.getType().equals("CDS")) {
					f = Ft;
					break;
				}
		}
		return f;
	}

	public static String ComplementaryTripletName(int n) {
		StringBuffer s = new StringBuffer();
		int i1 = (int) ((float) n / 16);
		int i2 = (int) ((float) (n - i1 * 16) / 4);
		int i3 = n - i1 * 16 - i2 * 4;
		s.append(ComplementaryBaseName(i1));
		s.append(ComplementaryBaseName(i2));
		s.append(ComplementaryBaseName(i3));
		return s.toString();
	}

	public static String ComplementaryBaseName(int n) {
		String s = new String("");
		switch (n) {
		case 0:
			s = new String("t");
			break;
		case 1:
			s = new String("g");
			break;
		case 2:
			s = new String("c");
			break;
		case 3:
			s = new String("a");
			break;
		}
		return s;
	}

	public static float FreqReconstruct(float fr[], String codon, int phase) {
		float f = 0;
		// out.println(codon);
		for (int i = 0; i < 64; i++)
			for (int j = 0; j < 64; j++) {
				String hexa = new String(TripletName(i) + TripletName(j));
				if (hexa.substring(phase, phase + 3).equals(codon)) {
					// out.println(hexa+" "+fr[i]+" "+fr[j]+" "+f);
					f += fr[i] * fr[j];
				}
			}
		return f;
	}

	public static float[] ReconstructPhaseDistribution(float distr[], int phase) {
		float f[] = new float[64];
		for (int i = 0; i < 64; i++)
			f[i] = FreqReconstruct(distr, TripletName(i), phase);
		return f;
	}

	public static float[] SymmetrizedDistribution(float f1[]) {
		float f[] = new float[64];
		float crf[] = ComplementaryDistribution(f1);
		for (int i = 0; i < f1.length; i++)
			f[i] = (f1[i] + crf[i]) / 2f;
		return f;
	}

	public static float[] AntiSymmetrizedDistribution(float f1[]) {
		float f[] = new float[64];
		float crf[] = ComplementaryDistribution(f1);
		for (int i = 0; i < f1.length; i++)
			f[i] = (f1[i] - crf[i]) / 2f;
		return f;
	}

	public static float DistributionDistance(float f1[], float f2[]) {
		double r = 0;
		for (int i = 0; i < 64; i++)
			r += Math.abs(f1[i] - f2[i]);// *(f1[i]-f2[i]);
		// r+=(f1[i]-f2[i])*(f1[i]-f2[i]);
		// return (float)Math.sqrt(r);
		return (float) r;
	}

	public static float DistributionDistanceEucl(float f1[], float f2[]) {
		double r = 0;
		for (int i = 0; i < 64; i++)
			r += (f1[i] - f2[i]) * (f1[i] - f2[i]);
		return (float) Math.sqrt(r);
	}

	public static float DistributionNorm(float f1[]) {
		double r = 0;
		for (int i = 0; i < 64; i++)
			r += Math.abs(f1[i]);// *f1[i];
		// r+=f1[i]*f1[i];
		// return (float)Math.sqrt(r);
		return (float) r;
	}

	public static float DistributionNormEucl(float f1[]) {
		double r = 0;
		for (int i = 0; i < 64; i++)
			r += f1[i] * f1[i];
		return (float) Math.sqrt(r);
	}

	public static float ReconstructionQuality(float distr[], float realdist[],
			int phase) {
		float f[] = ReconstructPhaseDistribution(distr, phase);
		return DistributionDistance(realdist, f) / DistributionNorm(distr);
	}

	public static float EvaluateDistributionGoodness(float f[]) {
		float p[] = ReconstructPhaseDistribution(f, 1);
		return DistributionDistance(f, p) / DistributionNorm(f) / 2;
	}

	public static float[] ComplementaryDistribution(float fr[]) {
		float f[] = new float[64];
		for (int i = 0; i < 64; i++) {
			StringBuffer s = new StringBuffer(ComplementaryTripletName(i));
			s.reverse();
			f[TripletNumber(s.toString())] = fr[i];
		}
		return f;
	}

	public static void printDistribution(float f[]) {
		DecimalFormat df = new DecimalFormat();
		df.applyPattern("#,###0.000");
		for (int i = 0; i < 8; i++) {
			for (int j = 0; j < 8; j++) {
				out.print(TripletName(i * 8 + j) + ":"
						+ df.format(f[i * 8 + j]) + "\t");
			}
			out.println();
		}
	}

	public static void printDistribution(float f[], String let1, String let2,
			String let3, String let4) {
		DecimalFormat df = new DecimalFormat();
		df.applyPattern("#,###0.000");
		for (int i = 0; i < 8; i++) {
			for (int j = 0; j < 8; j++) {
				String s = TripletName(i * 8 + j, let1, let2, let3, let4);
				int k = TripletNumber(s);
				out.print(s + ":" + df.format(f[k]) + "\t");
			}
			out.println();
		}
	}

	public static void printDistributionToFile(String fn, float f[],
			String let1, String let2, String let3, String let4)
			throws IOException {
		PrintStream cuout = new PrintStream(new FileOutputStream(fn));
		DecimalFormat df = new DecimalFormat();
		df.applyPattern("#,###0.000");
		for (int i = 0; i < 8; i++) {
			for (int j = 0; j < 8; j++) {
				String s = TripletName(i * 8 + j, let1, let2, let3, let4);
				int k = TripletNumber(s);
				cuout.print(s + ":" + df.format(f[k]) + "\t");
			}
			cuout.println();
		}
		cuout.close();
	}

	public static float ComplementaryAsymmetry(float fr0[], float fr1[],
			float fr2[], float cfr0[], float cfr1[], float cfr2[]) {
		float r = 0;
		float d00 = DistributionDistance(fr0, cfr0);
		float d01 = DistributionDistance(fr0, cfr1);
		float d02 = DistributionDistance(fr0, cfr2);
		float d10 = DistributionDistance(fr1, cfr0);
		float d11 = DistributionDistance(fr1, cfr1);
		float d12 = DistributionDistance(fr1, cfr2);
		float d20 = DistributionDistance(fr2, cfr0);
		float d21 = DistributionDistance(fr2, cfr1);
		float d22 = DistributionDistance(fr2, cfr2);
		float r1 = Math.min(d00, Math.min(d01, d02));
		float r2 = Math.min(d10, Math.min(d11, d12));
		float r3 = Math.min(d20, Math.min(d21, d22));
		r = (r1 + r2 + r3) / 3;
		return r;
	}

	public static float[] CalcJunkDistribution(float fr0[], float fr1[],
			float fr2[], float cfr0[], float cfr1[], float cfr2[]) {
		float res[] = new float[64];
		for (int i = 0; i < 64; i++)
			res[i] = (fr0[i] + fr1[i] + fr2[i] + cfr0[i] + cfr1[i] + cfr2[i]) / 6;
		return res;
	}

	public static String StringSpaces(int n) {
		StringBuffer ss = new StringBuffer("");
		for (int i = 0; i < n; i++)
			ss.append(" ");
		return ss.toString();
	}

	public static String getJunkString(Sequence seq) {
		StringBuffer ss = new StringBuffer(seq.seqString());
		for (Iterator i = seq.features(); i.hasNext();) {
			Feature Ft = (Feature) i.next();
			int i1 = Ft.getLocation().getMin();
			int i2 = Ft.getLocation().getMax() + 1;
			if (!(Ft.getType().equals("source")))
				if (Ft.getType().equals("gene")) {
					ss.replace(i1 - 1, i2 - 1, StringSpaces(i2 - i1));
				}
		}
		StringBuffer s1 = new StringBuffer("");
		for (int i = 0; i < ss.length(); i++)
			if (ss.charAt(i) != ' ')
				s1.append(ss.charAt(i));
		return s1.toString();
	}

	public static String getJunkStringSpaced(Sequence seq) {
		StringBuffer ss = new StringBuffer(seq.seqString());
		for (Iterator i = seq.features(); i.hasNext();) {
			Feature Ft = (Feature) i.next();
			int i1 = Ft.getLocation().getMin();
			int i2 = Ft.getLocation().getMax() + 1;
			if (!(Ft.getType().equals("source")))
				if (Ft.getType().equals("gene")) {
					ss.replace(i1 - 1, i2 - 1, StringSpaces(i2 - i1));
				}
		}
		return ss.toString();
	}

	public static String FormatSeq(StringBuffer s) {
		StringBuffer ss = new StringBuffer();
		int n = (int) ((float) s.length() / 60);
		for (int i = 0; i < n; i++) {
			ss.append("    " + (i * 60 + 1) + " ");
			for (int j = 0; j < 6; j++)
				ss.append(s.substring(i * 60 + j * 10, i * 60 + j * 10 + 10)
						+ " ");
			ss.append("\n");
		}
		return ss.toString();
	}

	public static byte[] nucleotideBinaryCode(String s) {
		byte res[] = new byte[4];
		if (s.equals("a")) {
			res[0] = 1;
			res[1] = 0;
			res[2] = 0;
			res[3] = 0;
		}
		if (s.equals("c")) {
			res[0] = 0;
			res[1] = 1;
			res[2] = 0;
			res[3] = 0;
		}
		if (s.equals("g")) {
			res[0] = 0;
			res[1] = 0;
			res[2] = 1;
			res[3] = 0;
		}
		if (s.equals("t")) {
			res[0] = 0;
			res[1] = 0;
			res[2] = 0;
			res[3] = 1;
		}
		return res;
	}

	public static byte[] stringBinaryCode(String s) {
		byte res[] = new byte[4 * s.length()];
		for (int i = 0; i < s.length(); i++) {
			byte r[] = nucleotideBinaryCode(s.substring(i, i + 1));
			res[i * 4] = r[0];
			res[i * 4 + 1] = r[1];
			res[i * 4 + 2] = r[2];
			res[i * 4 + 3] = r[3];
		}
		return res;
	}

	public static String generateRandomString(int len, float pa, float pc,
			float pg, float pt, Random r) {
		StringBuffer res = new StringBuffer();
		for (int i = 0; i < len; i++) {
			float rr = r.nextFloat();
			if (rr < pa)
				res.append("a");
			if ((rr >= pa) && (rr < pa + pc))
				res.append("c");
			if ((rr >= pa + pc) && (rr < pa + pc + pg))
				res.append("g");
			if ((rr > pa + pc + pg) && (rr <= 1))
				res.append("t");
		}
		return res.toString();
	}

	public static String mutateStringChange(String s, int chn, float pa,
			float pc, float pg, float pt, Random r) {
		StringBuffer sb = new StringBuffer(s);
		for (int i = 0; i < chn; i++) {
			char c = chooseLetter(100, (int) (pa * 100), (int) (pc * 100),
					(int) (pg * 100), (int) (pt * 100));
			double d = r.nextDouble();
			int ps = (int) (s.length() * d);
			if ((ps >= 0) && (ps < s.length())) {
				sb.setCharAt(ps, c);
			}
		}
		return sb.toString();
	}

	private static char chooseLetter(int len, int ac, int cc, int gc, int tc) {
		Random r = new Random();
		char res = 'a';
		double d = r.nextDouble();
		if (d < ac / len)
			res = 'a';
		if ((d >= ac / len) && (d < (ac + cc) / len))
			res = 'c';
		if ((d >= (ac + cc) / len) && (d < (ac + cc + gc) / len))
			res = 'g';
		if ((d >= (ac + cc + gc) / len) && (d < (ac + cc + gc + tc) / len))
			res = 't';
		return res;
	}

	public static String complementString(String s) {
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < s.length(); i++) {
			char c = s.charAt(i);
			if (c == 'a')
				c = 't';
			else if (c == 't')
				c = 'a';
			else if (c == 'g')
				c = 'c';
			else if (c == 'c')
				c = 'g';
			else if (c == 'A')
				c = 'T';
			else if (c == 'T')
				c = 'A';
			else if (c == 'G')
				c = 'C';
			else if (c == 'C')
				c = 'G';
			sb.append(c);
		}
		return sb.toString();
	}

	public static String reverseString(String s) {
		StringBuffer sb = new StringBuffer();
		for (int i = s.length() - 1; i >= 0; i--) {
			char c = s.charAt(i);
			sb.append(c);
		}
		return sb.toString();
	}

	public static float[] calcGCinPhase(float CodonUs[]) {
		float ret[] = new float[3];
		for (int i = 0; i < 64; i++) {
			String cn = TripletName(i);
			// System.out.println(cn+" "+CodonUs[i]);
			for (int j = 0; j < 3; j++) {
				// if((cn.charAt(j)=='a')||(cn.charAt(j)=='t'))
				// ret[j]-=CodonUs[i];
				if ((cn.charAt(j) == 'g') || (cn.charAt(j) == 'c'))
					ret[j] += CodonUs[i];
			}
		}
		return ret;
	}

	public static float[] calcNinPhase(float CodonUs[], String bp) {
		float ret[] = new float[3];
		for (int i = 0; i < 64; i++) {
			String cn = TripletName(i);
			// System.out.println(cn+" "+CodonUs[i]);
			for (int j = 0; j < 3; j++) {
				// if((cn.charAt(j)=='a')||(cn.charAt(j)=='t'))
				// ret[j]-=CodonUs[i];
				if ((cn.charAt(j) == bp.charAt(0)))
					ret[j] += CodonUs[i];
			}
		}
		return ret;
	}

	public static int calcPhase(String gene, int pos) {
		int phase = 0;
		int ps = (int) ((float) pos / 3.0f);
		phase = pos - ps * 3;
		return phase;
	}

	public static void makeClassStatisticsFiles(Vector classes, Vector points,
			Vector classList) {
		try {

			float pf[] = (float[]) points.elementAt(0);
			System.out.println("Dim = " + pf.length);

			PrintStream centr = new PrintStream(
					new FileOutputStream("cl.centr"));
			PrintStream distances = new PrintStream(new FileOutputStream(
					"cl.dist"));
			PrintStream names = new PrintStream(
					new FileOutputStream("cl.names"));
			PrintStream stds = new PrintStream(new FileOutputStream("cl.stds"));

			centroids = new Vector();
			deviations = new Vector();
			sizes = new Vector();

			if (classList == null)
				classList = new Vector();
			else {
				for (int i = 0; i < classList.size(); i++) {
					float p[] = (float[]) points.elementAt(0);
					float f[] = new float[p.length];
					centroids.add(f);
					deviations.add(new Float(0f));
					sizes.add(new Integer(1));
				}
			}

			for (int i = 0; i < classes.size(); i++) {
				String cl = (String) classes.elementAt(i);
				float p[] = (float[]) points.elementAt(i);
				int i1 = classList.indexOf(cl);
				if (i1 == -1) {
					classList.add(cl);
					float f1[] = new float[p.length];
					float f2[] = new float[p.length];
					for (int j = 0; j < p.length; j++) {
						f1[j] = p[j];
					}
					centroids.add(f1);
					deviations.add(new Float(0f));
					sizes.add(new Integer(1));
				} else {
					float f1[] = (float[]) centroids.elementAt(i1);
					for (int k = 0; k < p.length; k++)
						f1[k] += p[k];
					centroids.setElementAt(f1, i1);
					sizes.setElementAt(
							new Integer(((Integer) sizes.elementAt(i1))
									.intValue() + 1), i1);
				}
			}
			for (int i = 0; i < classList.size(); i++) {
				float p[] = (float[]) centroids.elementAt(i);
				for (int k = 0; k < p.length; k++)
					p[k] /= ((Integer) sizes.elementAt(i)).intValue();
			}
			for (int i = 0; i < classes.size(); i++) {
				String cl = (String) classes.elementAt(i);
				float p[] = (float[]) points.elementAt(i);
				int i1 = classList.indexOf(cl);
				float pp[] = (float[]) centroids.elementAt(i1);
				float d2 = 0f;
				for (int k = 0; k < p.length; k++)
					d2 += (pp[k] - p[k]) * (pp[k] - p[k]);
				float d2p = ((Float) deviations.elementAt(i1)).floatValue();
				// d2 = (float)Math.sqrt(d2);
				deviations.setElementAt(new Float(d2p + d2), i1);
				// System.out.println(Math.sqrt(d2));

				/*
				 * System.out.print(cl+"\t"); for(int
				 * j=0;j<classList.size();j++){ pp =
				 * (float[])centroids.elementAt(j); d2 = 0f; for(int
				 * k=0;k<p.length;k++) d2+=(pp[k]-p[k])*(pp[k]-p[k]);
				 * System.out.print(Math.sqrt(d2)+"\t"); } System.out.println();
				 */

			}

			for (int i = 0; i < classList.size(); i++) {
				float d2 = ((Float) deviations.elementAt(i)).floatValue();
				d2 /= ((Integer) sizes.elementAt(i)).intValue();
				deviations.setElementAt(new Float(Math.sqrt(d2)), i);
			}

			for (int i = 0; i < classList.size(); i++) {
				float p[] = (float[]) centroids.elementAt(i);
				String name = (String) classList.elementAt(i);
				float d = ((Float) deviations.elementAt(i)).floatValue();
				d *= (float) Math.sqrt(2f / p.length);
				for (int k = 0; k < p.length; k++)
					centr.print(p[k] + " ");
				centr.println();
				names.println(name);
				stds.println(d + "\t0");
				for (int k = 0; k < classList.size(); k++) {
					float f1[] = (float[]) centroids.elementAt(k);
					float dd = 0f;
					for (int kk = 0; kk < p.length; kk++)
						dd += (f1[kk] - p[kk]) * (f1[kk] - p[kk]);
					distances.print(Math.sqrt(dd) + "\t");
					// distances.print(DistributionDistanceEuclidean(f1,p)+"\t");
				}
				distances.println();
			}

			for (int i = 0; i < classList.size(); i++) {
				float pi[] = (float[]) centroids.elementAt(i);
				String clas = (String) classList.elementAt(i);
				for (int j = 0; j < classList.size(); j++)
					if (j != i) {
						float pos = 0f;
						float neg = 0f;
						int posc = 0;
						int negc = 0;
						float pj[] = (float[]) centroids.elementAt(j);
						// for(int k=0;k<pi.length;k++)
						// System.out.print(pj[k]+" ");
						// System.out.println();
						float d = 0f;
						for (int k = 0; k < pi.length; k++)
							d += (pi[k] - pj[k]) * (pi[k] - pj[k]);
						float n[] = new float[64];
						d = (float) Math.sqrt(d);
						// System.out.println(d);
						for (int k = 0; k < pi.length; k++)
							n[k] = (pj[k] - pi[k]) / d;
						for (int l = 0; l < classes.size(); l++) {
							String cl = (String) classes.elementAt(l);
							float p[] = (float[]) points.elementAt(l);
							if (cl.equals(clas)) {
								float proj = 0f;
								for (int m = 0; m < p.length; m++)
									proj += n[m] * (p[m] - pi[m]);
								if (proj > 0) {
									pos += proj;
									posc++;
								} else {
									neg += -proj;
									negc++;
								}
							}
						}
						// System.out.println("pos = "+pos+" posc = "+posc+" neg = "+neg+" negc = "+negc);
						if (posc != 0)
							pos /= posc;
						if (negc != 0)
							neg /= negc;
						stds.println(pos + "\t" + neg);
					}
			}

			centr.close();
			distances.close();
			names.close();
			stds.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static float[] calcDupletsFromSinglesQ(float f2[], float f1[]) {
		float res[] = new float[16];
		for (int i = 0; i < 16; i++) {
			String db = DipletName(i);
			String base1 = db.substring(0, 1);
			String base2 = db.substring(1, 2);
			float f = f1[BaseNumber(base1)] * f1[BaseNumber(base2)];
			if (f > 1e-10)
				res[i] = (f2[i] - f) / f;
		}
		return res;
	}

	public static float[] calcTripletsFromDipletsQ(float f3[], float f2[],
			float f1[]) {
		float res[] = new float[64];
		for (int i = 0; i < 64; i++) {
			String tr = TripletName(i);
			String base2 = tr.substring(1, 2);
			String db1 = tr.substring(0, 2);
			String db2 = tr.substring(1, 3);
			float f = f2[DipletNumber(db1)] * f2[DipletNumber(db2)]
					/ f1[BaseNumber(base2)];
			if (f > 1e-10)
				res[i] = (f3[i] - f) / f;
		}
		return res;
	}

}