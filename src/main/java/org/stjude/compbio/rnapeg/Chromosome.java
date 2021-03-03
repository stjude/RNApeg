package org.stjude.compbio.rnapeg;
// chromosome standardization.
// should probably just be replaced by code that processes
// (A) BAM header and (B) reference sequence list to juggle combinations.

import java.util.*;

public enum Chromosome {
  CHR1, CHR2, CHR3, CHR4, CHR5, CHR6, CHR7, CHR8, CHR9, CHR10,
    CHR11, CHR12, CHR13, CHR14, CHR15, CHR16, CHR17, CHR18, CHR19, CHR20,
    CHR21, CHR22, CHR23, CHR24, CHR25, CHRX, CHRY, CHRM, 
    // hot mess #2: drosophia:
    CHR2LHet, CHR2RHet, CHR3LHet, CHR3RHet, CHR2L, CHR2R,
    CHR3L, CHR3R, CHRU, CHRUextra, CHRXHet, CHRYHet;
    // drosophila end

  static HashMap<String,Chromosome> chr_map = null;
  static boolean VERBOSE = false;

  public static Chromosome valueOfByte (byte b) {
    Chromosome result;
    switch (b) {
    case 1: result = CHR1; break;
    case 2: result = CHR2; break;
    case 3: result = CHR3; break;
    case 4: result = CHR4; break;
    case 5: result = CHR5; break;
    case 6: result = CHR6; break;
    case 7: result = CHR7; break;
    case 8: result = CHR8; break;
    case 9: result = CHR9; break;
    case 10: result = CHR10; break;
    case 11: result = CHR11; break;
    case 12: result = CHR12; break;
    case 13: result = CHR13; break;
    case 14: result = CHR14; break;
    case 15: result = CHR15; break;
    case 16: result = CHR16; break;
    case 17: result = CHR17; break;
    case 18: result = CHR18; break;
    case 19: result = CHR19; break;
    case 20: result = CHR20; break;
    case 21: result = CHR21; break;
    case 22: result = CHR22; break;

    case 23: result = CHR23; break;
    case 24: result = CHR24; break;
    case 25: result = CHR25; break;
      // zebrafish

    case 'X': result = CHRX; break;
    case 'Y': result = CHRY; break;
    case 'M': result = CHRM; break;
    default: result = null;
    }
    return result;
  }

  public static Chromosome valueOfString (String key) {
    if (chr_map == null) {
      chr_map = new HashMap<String,Chromosome>();
      chr_map.put("1", CHR1);
      chr_map.put("2", CHR2);
      chr_map.put("3", CHR3);
      chr_map.put("4", CHR4);
      chr_map.put("5", CHR5);
      chr_map.put("6", CHR6);
      chr_map.put("7", CHR7);
      chr_map.put("8", CHR8);
      chr_map.put("9", CHR9);
      chr_map.put("10", CHR10);
      chr_map.put("11", CHR11);
      chr_map.put("12", CHR12);
      chr_map.put("13", CHR13);
      chr_map.put("14", CHR14);
      chr_map.put("15", CHR15);
      chr_map.put("16", CHR16);
      chr_map.put("17", CHR17);
      chr_map.put("18", CHR18);
      chr_map.put("19", CHR19);
      chr_map.put("20", CHR20);
      chr_map.put("21", CHR21);
      chr_map.put("22", CHR22);

      chr_map.put("23", CHR23);
      chr_map.put("24", CHR24);
      chr_map.put("25", CHR25);
      // hot mess:
      // - zebrafish has chr23-25
      // - human sometimes 23=X 24=Y 25=MT

      chr_map.put("x", CHRX);
      chr_map.put("y", CHRY);
      chr_map.put("m", CHRM);
      chr_map.put("mt", CHRM);

      // drosophila:
      chr_map.put("2L", CHR2L);
      chr_map.put("2R", CHR2R);
      chr_map.put("2LHet", CHR2LHet);
      chr_map.put("2RHet", CHR2RHet);
      chr_map.put("3L", CHR3L);
      chr_map.put("3R", CHR3R);
      chr_map.put("3LHet", CHR3LHet);
      chr_map.put("3RHet", CHR3RHet);
      chr_map.put("U", CHRU);
      chr_map.put("Uextra", CHRUextra);
      chr_map.put("XHet", CHRXHet);
      chr_map.put("YHet", CHRYHet);

      ArrayList<String> labels = new ArrayList<String>(chr_map.keySet());
      // separate list to avoid concurrent modification

      ArrayList<String> prefixes = new ArrayList<String>();
      prefixes.add("chr");
      // drosophila:
      prefixes.add("arm_");
      // drosophila: arm_2L
      prefixes.add("Arm");
      // drosophila: armU

      for (String label : labels) {
	Chromosome c = chr_map.get(label);
	for (String prefix : prefixes) {
	  String label2 = prefix + label;
	  if (VERBOSE) System.err.println("add label for " + label + " => " + label2 + " => " + c);  // debug
	  chr_map.put(label2, c);
	}
      }

      


    }

    Chromosome hit = chr_map.get(key);
    if (hit == null) hit = chr_map.get(key.toLowerCase());
    return hit;

  }

  public String toString() {
    String result = null;
    switch (this) {
    case CHR1: result = "chr1"; break;
    case CHR2: result = "chr2"; break;
    case CHR3: result = "chr3"; break;
    case CHR4: result = "chr4"; break;
    case CHR5: result = "chr5"; break;
    case CHR6: result = "chr6"; break;
    case CHR7: result = "chr7"; break;
    case CHR8: result = "chr8"; break;
    case CHR9: result = "chr9"; break;
    case CHR10: result = "chr10"; break;
    case CHR11: result = "chr11"; break;
    case CHR12: result = "chr12"; break;
    case CHR13: result = "chr13"; break;
    case CHR14: result = "chr14"; break;
    case CHR15: result = "chr15"; break;
    case CHR16: result = "chr16"; break;
    case CHR17: result = "chr17"; break;
    case CHR18: result = "chr18"; break;
    case CHR19: result = "chr19"; break;
    case CHR20: result = "chr20"; break;
    case CHR21: result = "chr21"; break;
    case CHR22: result = "chr22"; break;

    case CHR23: result = "chr23"; break;
    case CHR24: result = "chr24"; break;
    case CHR25: result = "chr25"; break;
      // zebrafish

    case CHRX: result = "chrX"; break;
    case CHRY: result = "chrY"; break;
    case CHRM: result = "chrM"; break;

      // drosophila:
    case CHR2L: result = "chr2L"; break;
    case CHR2R: result = "chr2R"; break;
    case CHR2LHet: result = "chr2LHet"; break;
    case CHR2RHet: result = "chr2RHet"; break;
      
    case CHR3L: result = "chr3L"; break;
    case CHR3R: result = "chr3R"; break;
    case CHR3LHet: result = "chr3LHet"; break;
    case CHR3RHet: result = "chr3RHet"; break;

    case CHRU: result = "chrU"; break;
    case CHRUextra: result = "chrUextra"; break;
    case CHRXHet: result = "chrXHet"; break;
    case CHRYHet: result = "chrYHet"; break;

    }
    return result;
  }

  public int toInt() {
    int result = 0;
    switch (this) {
    case CHR1: result = 1; break;
    case CHR2: result = 2; break;
    case CHR3: result = 3; break;
    case CHR4: result = 4; break;
    case CHR5: result = 5; break;
    case CHR6: result = 6; break;
    case CHR7: result = 7; break;
    case CHR8: result = 8; break;
    case CHR9: result = 9; break;
    case CHR10: result = 10; break;
    case CHR11: result = 11; break;
    case CHR12: result = 12; break;
    case CHR13: result = 13; break;
    case CHR14: result = 14; break;
    case CHR15: result = 15; break;
    case CHR16: result = 16; break;
    case CHR17: result = 17; break;
    case CHR18: result = 18; break;
    case CHR19: result = 19; break;
    case CHR20: result = 20; break;
    case CHR21: result = 21; break;
    case CHR22: result = 22; break;
    case CHRX: result = 23; break;
    case CHRY: result = 24; break;
    case CHRM: result = 25; break;
    }
    return result;
  }

  public static String standardize_name (String name) {
    Chromosome std = Chromosome.valueOfString(name);
    //    System.err.println("raw="+name + " std=" + std);  // debug

    return std == null ? name : std.toString();
  }

  public static void main (String[] argv) {
    ArrayList<String> tests = new ArrayList<String>();
    tests.add("1");
    tests.add("chr1");
    tests.add("Chr1");
    tests.add("CHR1");

    tests.add("x");
    tests.add("X");
    tests.add("chrX");
    tests.add("ChrX");
    tests.add("CHRX");

    tests.add("NM_12345");
    tests.add("some_custom_format");
    tests.add("Zv9_NA919");
    tests.add("23");
    tests.add("24");
    tests.add("25");

    for (String raw : tests) {
      //      Chromosome c = Chromosome.valueOfString(raw);
      //      System.err.println(raw + " => " + c);  // debug
      System.err.println(raw + " => " + Chromosome.standardize_name(raw));  // debug
    }
  }


}

