package org.stjude.compbio.rnapeg;

import htsjdk.samtools.*;
import java.io.*;
import java.util.*;
import java.net.*;
import javax.swing.*;

public class SAMUtils {
  private final static boolean DEBUG_RESTRICT_LOAD = false;
  // private final static boolean DEBUG_RESTRICT_LOAD = false;

  public static ArrayList<String> get_refname_alternates (String s) {
    ArrayList<String> results = new ArrayList<String>();
    HashSet<String> saw = new HashSet<String>();
    conditional_add(results, saw, s);
    // unique ordered list, always returning given string first

    String std = get_standardized_refname(s);
    conditional_add(results, saw, std);
    conditional_add(results, saw, "chr" + std);

    if (std.toUpperCase().equals("M")) {
      conditional_add(results, saw, "MT");
      conditional_add(results, saw, "chrMT");
    } else if (std.toUpperCase().equals("MT")) {
      conditional_add(results, saw, "M");
      conditional_add(results, saw, "chrM");
    }
    
    //    System.err.println("input="+s + " out=" + results);  // debug

    return results;
  }

  public static void main (String[] argv) {
    for (String s : get_refname_alternates("MT")) {
      System.err.println("alt="+s);  // debug
    }
  }


  private static void conditional_add (ArrayList<String> list, HashSet<String> saw, String s) {
    if (!saw.contains(s)) {
      list.add(s);
      saw.add(s);
    }
  }

  public static String get_standardized_refname (String s) {
    if (s.toLowerCase().indexOf("chr") == 0) {
      return(s.substring(3));
    } else {
      return s;
    }
  }
    
}
