package org.stjude.compbio.rnapeg;
// FIX ME: move ENSEMBL transcript->gene lookup here, make an option

import java.util.*;

public class UCSCRefGene {
  //
  // a row in UCSC refGene table.
  //
  // tweaks for different quasi-compatible formats:
  // - AceView: parse name from transcript ID
  // - UCSC/ensGene.txt (ENSEMBL): needs lookup table to translate ENST* transcript ID to gene
  // - UCSC/knownGene_refflat.txt: ??? gene id format ??
  //   ALSO NEEDS LOOKUP FILE:
  //   http://10.4.19.34:8080/display/DTAWHS/How+to+Count+Junction+Reads+in+an+RNA-Seq+BAM#
  //   cut -f 1,5 /nfs_exports/genomes/1/Homo_sapiens/GRCh37-lite/mRNA/UCSC/kgXref.txt
  //
  // - UCSC/refGene.txt: OK
  //

  //
  // more-or-less raw fields from table:
  //
  int bin;
  public String name, chrom, strand;
  // name = refseq accession
  int txStart, txEnd;
  int cdsStart, cdsEnd;
  int exonCount;
  public int[] exonStarts, exonEnds;
  int score;
  public String name2;
  // name2 = gene symbol
  String cdsStartStat, cdsEndStat;
  int[] exonFrames;

  public UCSCRefGene() {
  }

  public UCSCRefGene(HashMap<String,String> h) {
    parse(h);
  }

  public UCSCRefGene(String line) {
    parse_string(line);
  }

  public void parse_string(String line) {
    String[] fields = line.split("\t");
    if (fields.length == 16 ||
	fields.length == 11 ||
	fields.length == 13) {
      //      System.err.println("line="+line);  // debug

      int i = 0;

      // FAIL for
      // /nfs_exports/genomes/1/Homo_sapiens/GRCh37-lite/mRNA/UCSC/knownGene_refFlat.txt 
      bin = forgiving_parse_int(fields[i++], -1);

      //      System.err.println("bin="+bin);  // debug
      name = new String(fields[i++]);
      //      System.err.println("name="+name);  // debug

      chrom = new String(fields[i++]);
      //      System.err.println("chrom="+chrom);  // debug

      strand = new String(fields[i++]);
      //      System.err.println("strand="+strand);  // debug

      txStart = Integer.parseInt(fields[i++]);
      //      System.err.println("txStart="+txStart);  // debug

      txEnd = Integer.parseInt(fields[i++]);
      //      System.err.println("txEnd="+txEnd);  // debug

      cdsStart = Integer.parseInt(fields[i++]);
      //      System.err.println("cdsStart="+txStart);  // debug

      cdsEnd = Integer.parseInt(fields[i++]);
      //      System.err.println("cdsEnd="+txEnd);  // debug

      exonCount = Integer.parseInt(fields[i++]);
      //      System.err.println("count="+exonCount);  // debug

      exonStarts = parse_int_range(fields[i++]);
      exonEnds = parse_int_range(fields[i++]);

      if (fields.length >= 13) {
	score = forgiving_parse_int(fields[i++], -1);
	name2 = new String(fields[i++]);
      }
      if (fields.length >= 16) {
	cdsStartStat = new String(fields[i++]);
	cdsEndStat = new String(fields[i++]);
	exonFrames = parse_int_range(fields[i++]);
      }

      if (fields.length == 11 && name.indexOf(".") > 0) {
	// AceView: gene name embedded in transcript name
	String[] f = name.split("\\.");
	if (f.length == 2) {
	  name2 = new String(f[0]);
	}
      }

    } else {
      System.err.println("ERROR: unknown format, expected 16 fields, has " + fields.length);  // debug
    }
  }

  public void parse(HashMap<String,String> row) {
    // FIX ME: style points for populating via Reflection!
    bin = Integer.parseInt(row.get("bin"));
    name = row.get("name");
    chrom = row.get("chrom");
    strand = row.get("strand");
    txStart = Integer.parseInt(row.get("txStart"));
    txEnd = Integer.parseInt(row.get("txEnd"));
    exonCount = Integer.parseInt(row.get("exonCount"));
    exonStarts = parse_int_range(row.get("exonStarts"));
    exonEnds = parse_int_range(row.get("exonEnds"));
    score = Integer.parseInt(row.get("score"));
    name2 = row.get("name2");
    cdsStartStat = row.get("cdsStartStat");
    cdsEndStat = row.get("cdsEndStat");
    exonFrames = parse_int_range(row.get("exonFrames"));
  }

  public static void main (String[] argv) {
    try {
      JDBCQuery ucsc = JDBCQuery.get_stjude_hg19();
      //      HashMap<String,String> hit = ucsc.query_single_row("select * from refGene where name2=\"WRN\"");
      HashMap<String,String> hit = ucsc.query_single_row("refGene", "name2", "WRN");
      System.err.println("hit="+hit);  // debug
      UCSCRefGene rg = new UCSCRefGene(hit);
      System.err.println(rg.name2);  // debug
      System.err.println(rg.exonStarts);  // debug
    } catch (Exception e) {
      System.err.println("ERROR: "+e);  // debug
    }
  }

  private int[] parse_int_range (String range) {
    String[] ranges = range.split(",");
    int[] results = new int[ranges.length];
    for (int i = 0; i < ranges.length; i++) {
      results[i] = Integer.parseInt(ranges[i]);
    }
    return results;
  }

  private int forgiving_parse_int (String v, int broken_value) {
    int result;
    if (v == null || v.length() == 0) {
      result = broken_value;
    } else {
      try {
	result = Integer.parseInt(v);
      } catch (NumberFormatException nfe) {
	result = broken_value;
      }
    }
    return result;
  }


}
