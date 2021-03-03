package org.stjude.compbio.rnapeg;

import java.util.*;
import java.io.*;
import htsjdk.samtools.cram.ref.ReferenceSource;

public class FASTAIndexedFAI implements ReferenceSequence {
  //
  // access a set of sequences in a FASTA file via samtools .fai index format
  //
  RandomAccessFile raf;
  HashMap<String,FAIIndexRecord> index;
  FASTAFileTools ffp = new FASTAFileTools();
  ReferenceSource rs = null;

  public FASTAIndexedFAI (String fn) throws FileNotFoundException,IOException {
    raf = new RandomAccessFile(fn, "r");
    rs = new ReferenceSource(new File(fn));
    load_index(fn + ".fai");
  }

  private void load_index (String fn) throws FileNotFoundException,IOException {
    File f = new File(fn);
    index = new HashMap<String,FAIIndexRecord>();
    if (f.exists()) {
      //      System.err.println("load index");  // debug
      BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(f)));
      String line;
      while ((line = br.readLine()) != null) {
	//	System.err.println("line="+line);  // debug
	FAIIndexRecord rec = new FAIIndexRecord(line);
	if (index.containsKey(rec.sequence_id)) {
	  System.err.println("WARNING; duplicate index entry for " + rec.sequence_id);  // debug
	}
	index.put(rec.sequence_id, rec);
      }
    } else {
      System.err.println("ERROR: no samtools FASTA index file " + fn + "; generate with \"samtools faidx file\"");  // debug
    }
  }

  public byte[] get_region (String sequence_name, int start_base, int length) throws IOException {
    // fetch a region of a reference sequence
    byte[] results = null;
    FAIIndexRecord ir = get_index_for(sequence_name);
    if (ir != null) {
      ffp.set_index(ir);
      results = ffp.get_region(raf, start_base, length);
    }
    return results;
  }

  public byte[] get_all (String sequence_name) throws IOException {
    byte[] results = null;
    FAIIndexRecord ir = get_index_for(sequence_name);
    if (ir != null) {
      ffp.set_index(ir);
      results = ffp.get_all(raf);
    }
    return results;
  }

  private FAIIndexRecord get_index_for (String id) {
    FAIIndexRecord result = null;
    for (String key : SAMUtils.get_refname_alternates(id)) {
      result = index.get(key);
      if (result != null) break;
    }

    if (result == null) {
      ChromosomeDisambiguator cd = new ChromosomeDisambiguator(index.keySet());
      String key = cd.find(id);
      if (key != null) result = index.get(key);
    }

    if (result == null) {
      System.err.println("WARNING: can't find .fai index entry for " + id);  // debug
    }
    return result;
  }


  public int get_length (String sequence_name) throws IOException {
    // sequence length
    int length = ReferenceSequence.NULL_LENGTH;
    FAIIndexRecord ir = get_index_for(sequence_name);
    if (ir != null) length = ir.sequence_length;
    return length;
  }

  public static void main (String[] argv) {
    try {
      FASTAIndexedFAI fai = new FASTAIndexedFAI("ecoli_out.padded.fasta");
      //      String id = "sequence2_10";
      String id = "NC_000913_bb";
      System.err.println("len=" + fai.get_length(id));
      byte[] all = fai.get_all(id);
      System.err.println("all="+new String(all));  // debug

      byte[] region = fai.get_region(id, 131, 10);
      System.err.println("region="+new String(region));  // debug

    } catch (Exception e) {
      System.err.println("ERROR: " + e);  // debug
    }
  }

  public boolean supports_sequence_list() {
    return true;
  }

  public ArrayList<String> get_sequence_names() {
    return new ArrayList<String>(index.keySet());
  }

  public ReferenceSource getReferenceSource() {
    return rs;
  }


  
}