package org.stjude.compbio.rnapeg;
// reference/chromosome name disambiguator.
//
// should probably just be replaced by code that processes
// (A) BAM header and (B) reference sequence list to juggle combinations.

import java.util.*;
import htsjdk.samtools.*;

public class ChromosomeDisambiguator {
  HashMap<String,String> map;
  HashMap<String,Integer> lengths;

  public ChromosomeDisambiguator (ArrayList<String> list) {
    setup(list);
  }

  public ChromosomeDisambiguator (Set<String> list) {
    setup(new ArrayList<String>(list));
  }

  public ChromosomeDisambiguator (SamReader sfr) {
    SAMFileHeader h = sfr.getFileHeader();
    SAMSequenceDictionary dict = h.getSequenceDictionary();
    ArrayList<String> refs = new ArrayList<String>();
    lengths = new HashMap<String,Integer>();
    for (SAMSequenceRecord ssr : dict.getSequences()) {
      String ref_name = ssr.getSequenceName();
      refs.add(ref_name);
      lengths.put(ref_name, ssr.getSequenceLength());
      //      System.err.println("ref="+ref_name);  // debug
    }
    setup(refs);
  }

  public ChromosomeDisambiguator (SAMFileHeader h) {
    // don't ask
    SAMSequenceDictionary dict = h.getSequenceDictionary();
    ArrayList<String> refs = new ArrayList<String>();
    lengths = new HashMap<String,Integer>();
    for (SAMSequenceRecord ssr : dict.getSequences()) {
      String ref_name = ssr.getSequenceName();
      refs.add(ref_name);
      lengths.put(ref_name, ssr.getSequenceLength());
      //      System.err.println("ref="+ref_name);  // debug
    }
    setup(refs);
  }


  private void setup (ArrayList<String> list) {
    map = new HashMap<String,String>();
    for (String raw : list) {
      String cooked = Chromosome.standardize_name(raw);
      map.put(cooked, raw);
      //      System.err.println("mapping " + cooked + " => " + raw);  // debug

      // map standard name to name appearing in raw list
    }
  }

  public String find (String raw) {
    System.err.println("raw="+raw);  // debug

    return map.get(Chromosome.standardize_name(raw));
  }

  public int get_length (String raw) {
    String cooked = find(raw);
    return cooked == null ? -1 : lengths.get(cooked);
  }

}