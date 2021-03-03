package org.stjude.compbio.rnapeg;
// parse a refGene flatfile, w/optional patching of gene symbols

import java.io.*;
import java.util.*;

public class UCSCRefGeneReader implements Iterator,Iterable<UCSCRefGene> {
  BufferedReader br;
  String next_line;
  HashMap<String,String> transcript2gene = null;

  public UCSCRefGeneReader(String file) throws FileNotFoundException,IOException {
    br = FileUtils.getUniversalBufferedReader(file);
    next_line = null;
  }

  public UCSCRefGeneReader(File file) throws FileNotFoundException,IOException {
    br = FileUtils.getUniversalBufferedReader(file);
    next_line = null;
  }
  
  public void parse (String file) throws FileNotFoundException,IOException {
    UCSCRefGene rg;
    String line;
    while ((line = br.readLine()) != null) {
      rg = new UCSCRefGene(line);
    }
  }

  // begin Iterator stub
  public boolean hasNext() {
    if (next_line == null) try {
	next_line = br.readLine();
      } catch (Exception e) {
	System.err.println("ERROR reading from refgene file");  // debug
	next_line = null;
      }
    return next_line != null;
  }

  public UCSCRefGene next() {
    UCSCRefGene rg = null;
    if (hasNext()) {
      // have a line to parse
      rg = new UCSCRefGene(next_line);
      if (transcript2gene != null) {
	String lookup = transcript2gene.get(rg.name);
	if (lookup == null) {
	  System.err.println("WTF, can't find gene map for " + rg.name);  // debug
	} else {
	  rg.name2 = lookup;
	}
      }
      next_line = null;
      // mark done
    }
    return rg;
  }

  public void set_transcript2gene (String map_file) {
    // use lookup flatfile converting transcript IDs to gene symbols
    transcript2gene = new HashMap<String,String>();
    try {
      BufferedReader br = FileUtils.getUniversalBufferedReader(map_file);
      String line;
      while ((line = br.readLine()) != null) {
	String[] f = line.split("\t");
	if (f.length == 2) {
	  transcript2gene.put(f[0], f[1]);
	} else {
	  System.err.println("transcript2gene format must be 2 fields");  // debug
	}
      }
    } catch (Exception e) {
      System.err.println("ERROR: " + e);  // debug
    }
  }

  public void remove() {
    System.err.println("remove() not implemented");  // debug
  }
  // end Iterator stub

  // begin Iterable stub
  public Iterator<UCSCRefGene> iterator() {
    return this;
  }
  // end Iterable stub

  public static void main (String[] argv) {
    try {
      //      UCSCRefGeneReader r = new UCSCRefGeneReader("refGene.txt.gz");
      UCSCRefGeneReader r = new UCSCRefGeneReader("ensGene.txt");
      r.set_transcript2gene("ensemblToGeneName.txt");
      // update name2 gene symbol annotation from lookup table

      for (UCSCRefGene rg : r) {
	System.err.println("read " + rg.name2);  // debug
      }

    } catch (Exception e) {
      System.err.println("ERROR: " +e);  // debug
    }
  }


}
