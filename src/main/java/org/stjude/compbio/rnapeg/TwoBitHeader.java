package org.stjude.compbio.rnapeg;
// .2bit files: http://genome.ucsc.edu/FAQ/FAQformat#format7
// implementation of header and sequence index section
// MNE 2/2010

import java.io.*;
import java.util.*;

public class TwoBitHeader {

  private static final int TWOBIT_SIGNATURE = 0x1A412743;
  // magic

  private File file;
  private boolean error = false;
  private int sequence_count, reserved;
  private boolean enable_caching;
  
  private ArrayList<String> sequences;
  // save list of IDs as ordered in file in case we need it
  private HashMap<String,Integer> sequence_offsets;
  // file offsets for sequence records
  private EndianReader er;
  private HashMap<String,TwoBitSequence> seqmap;

  public TwoBitHeader (File file) throws FileNotFoundException,IOException {
    this.file = file;
    setup();
  }

  public TwoBitHeader (String filename) throws FileNotFoundException,IOException {
    this.file = new File(filename);
    setup();
  }

  private RandomAccessFile open_randomaccessfile() throws IOException,FileNotFoundException {
    return new RandomAccessFile(file, "r");
  }

  private void setup() throws FileNotFoundException,IOException {
    enable_caching = true;
    sequences = new ArrayList<String>();
    sequence_offsets = new HashMap<String,Integer>();
    seqmap = new HashMap<String,TwoBitSequence>();

    RandomAccessFile raf = open_randomaccessfile();

    try {
      raf.seek(0);
      er = new BigEndianReader(raf);
      //
      //  detect file endianness
      //
      int sig = er.readInt();
      if (sig == TWOBIT_SIGNATURE) {
	System.err.println("big-endian .2bit file, UNTESTED");  // debug
      } else {
	raf.seek(0);
	er = new LittleEndianReader(raf);
	int sig2 = er.readInt();
	if (sig2 != TWOBIT_SIGNATURE) {
	  System.err.println("error: can't find .2bit signature");  // debug
	  System.err.println("sig1=" + sig + " sig2=" + sig2);  // debug
	  error = true;
	}
      }

      if (!error) {
	int version = er.readInt();
	if (version == 0) {
	  //
	  //  read sequence index
	  //
	  sequence_count = er.readInt();
	  int reserved = er.readInt();
	  byte len;
	  for (int i = 0; i < sequence_count; i++) {
	    len = raf.readByte();
	    byte[] buf = new byte[len];
	    if (raf.read(buf) == len) {
	      String s = new String(buf);
	      int offset = er.readInt();
	      sequences.add(s);
	      sequence_offsets.put(s, offset);
	      //	      System.err.println("id " + s + " at " + offset);  // debug
	    } else {
	      System.err.println("ERROR reading sequence ID!");  // debug
	      error = true;
	    }
	  }

	  long first_record_offset = er.getFilePointer();
	  long first_header_record_offset = sequence_offsets.get(sequences.get(0));

	  if (first_record_offset != first_header_record_offset) {
	    System.err.println("WARNING: first record position (i.e. after index) at " + first_record_offset + ", header says at " + first_header_record_offset + "; header corrupt??");  // debug
	  }

	} else {
	  System.err.println("ERROR: nonzero version " + version);  // debug
	  error = true;
	}
      }
      
    } catch (Exception e) {
      System.err.println("ERROR: " + e);  // debug
      e.printStackTrace();
      error = true;
    }
  }

  public void set_caching (boolean v) {
    enable_caching = v;
  }

  public TwoBitSequence get_sequence (String id) throws IOException {
    return get_sequence(id, true);
  }

  public TwoBitSequence get_sequence (String id, boolean read_fully) throws IOException {
    TwoBitSequence result = seqmap.get(id);

    RandomAccessFile raf = open_randomaccessfile();
    // feh: create a new RandomAccessFile for every sequence request.
    // Necessary because TwoBitSequence uses a BufferedReader
    // initialized from the RandomAccessFile's file descriptor, which
    // is automatically closed when the BufferedReader instance is
    // destroyed.  This leads to an invalid file handle error on
    // subsequent RandomAccessFile seek() attempts.

    if (result == null) {
      // not in cache
      Integer offset = null;
      // sequence_offsets.get(id);
      //      ArrayList<String> keys = 
      for (String key : SAMUtils.get_refname_alternates(id)) {
	offset = sequence_offsets.get(key);
	if (offset != null) {
	  //	  System.err.println("hit for " + key);  // debug
	  break;
	}
      }
      if (offset == null) {
	ChromosomeDisambiguator cd = new ChromosomeDisambiguator(sequence_offsets.keySet());
	String key = cd.find(id);
	if (key != null) offset = sequence_offsets.get(key);
      }

      if (offset == null) {
	//	System.err.println("error -- sequence ID not found in 2bit file: " + id);  // debug
	//	(new Exception()).printStackTrace();
      } else {
	er.set_randomaccessfile(raf);
	result = new TwoBitSequence(er, offset, read_fully);
	if (enable_caching) seqmap.put(id, result);
      }
    } else {
      result.set_randomaccessfile(raf);
    }

    return result;
  }

  public static void main (String[] argv) {
    String fn = "c:/generatable/hg18/hg18.2bit";
    try {
      if (argv.length == 3) {
	String chr_name = argv[0];
	int start = Integer.parseInt(argv[1]);
	int len = Integer.parseInt(argv[2]);

	TwoBitHeader tbh = new TwoBitHeader(new File(fn));
	TwoBitSequence tbs = tbh.get_sequence(chr_name);

	//      char[] seq = tbs.get_full_sequence();

	byte[] seq = tbs.get_region(start, len);

	System.err.println("sequence: " + new String(seq));  // debug
      } else {
	System.err.println("specify chr start len");  // debug
	System.exit(1);
      }

    } catch (Exception e) {
      System.err.println("ERROR: " + e);  // debug
      e.printStackTrace();
    }
  }

  public ArrayList<String> get_sequence_names() {
    return sequences;
  }

}