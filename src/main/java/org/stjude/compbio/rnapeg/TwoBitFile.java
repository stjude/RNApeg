package org.stjude.compbio.rnapeg;
// .2bit files: http://genome.ucsc.edu/FAQ/FAQformat#format7
// MNE 2/2010

import java.io.*;
import java.util.*;

public class TwoBitFile implements ReferenceSequence {
  private TwoBitHeader tbh;

  public TwoBitFile (String filename) throws FileNotFoundException,IOException {
    tbh = new TwoBitHeader(filename);
  }

  public byte[] get_region (String sequence_name, int start_base, int length) throws IOException {
    TwoBitSequence tbs = tbh.get_sequence(sequence_name);
    byte[] result = null;
    if (tbs == null) {
      //      System.err.println("WARNING: .2bit file doesn't contain a reference sequence named " + sequence_name);  // debug
    } else {
      result = tbs.get_region(start_base, length);
    }
    return result;
  }
  
  public byte[] get_all (String sequence_name) throws IOException {
    TwoBitSequence tbs = tbh.get_sequence(sequence_name);
    byte[] result = null;
    if (tbs == null) {
      //      System.err.println("WARNING: .2bit file doesn't contain a reference sequence named " + sequence_name);  // debug
    } else {
      result = tbs.get_full_sequence();
    }
    return result;
  }

  public int get_length (String sequence_name) throws IOException {
    TwoBitSequence tbs = tbh.get_sequence(sequence_name);
    int length = ReferenceSequence.NULL_LENGTH;
    if (tbs == null) {
      //      System.err.println("WARNING: .2bit file doesn't contain a reference sequence named " + sequence_name);  // debug
      // disable message, can be used passively to test existence
    } else {
      length = tbs.get_length();
    }
    return length;
  }

  public static void main (String[] argv) {
    String fn = "c:/generatable/hg18/hg18.2bit";
    try {
      if (argv.length == 4) {
	TwoBitFile tbf = new TwoBitFile(argv[0]);
	String chr_name = argv[1];
	int start = Integer.parseInt(argv[2]);
	int len = Integer.parseInt(argv[3]);
	byte[] seq = tbf.get_region(chr_name, start, len);
	System.out.println(new String(seq));  // debug
      } else if (argv.length == 3) {
	TwoBitFile tbf = new TwoBitFile(fn);
	String chr_name = argv[0];
	int start = Integer.parseInt(argv[1]);
	int len = Integer.parseInt(argv[2]);
	byte[] seq = tbf.get_region(chr_name, start, len);
	System.out.println(new String(seq));  // debug
      } else if (true) {
	TwoBitFile tbf = new TwoBitFile(fn);
	byte[] chr = tbf.get_all("chr1");
	System.err.println("chr length = " + chr.length);  // debug
	while (true) {
	  System.gc();
	  try {
	    System.err.println("sleeping...");  // debug
	    Thread.sleep(1000 * 10);
	  } catch (InterruptedException e) {}
	}
      } else {
	System.err.println("specify chr start len");  // debug
	System.exit(1);
      }

    } catch (Exception e) {
      System.err.println("ERROR: " + e);  // debug
      e.printStackTrace();
    }
  }

  public boolean supports_sequence_list() {
    return true;
  }

  public ArrayList<String> get_sequence_names() {
    return tbh.get_sequence_names();
  }



}
