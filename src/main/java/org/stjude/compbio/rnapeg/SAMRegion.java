package org.stjude.compbio.rnapeg;

import java.text.ParseException;

public class SAMRegion {
  public String tname = null;
  public Range range;
  public String gene_name = null;

  public SAMRegion() {
    //    System.err.println("new SAMRegion");  // debug
    //    new Exception().printStackTrace();
    range = new Range();
  }

  public SAMRegion (String tname) {
    range = new Range();
    this.tname = tname;
  }

  public SAMRegion(String tname, int start, int end) {
    //    System.err.println("new SAMRegion");  // debug
    //    new Exception().printStackTrace();
    this.tname = tname;
    range = new Range(start, end);
  }

  public boolean isValid() {
    return tname != null && range.isValid();
  }

  public String toString() {
    if (isValid()) {
      return tname + ":" + range.start + "-" + range.end;
    } else if (tname != null) {
      return tname;
    } else {
      return null;
    }
  }

  public String get_ucsc_chromosome () {
    Chromosome c = Chromosome.valueOfString(tname);
    return c == null ? null : c.toString();
  }

  public int get_length() {
    return (range.end - range.start) + 1;
  }

  public void set_start (int v) {
    range.start = v;
  }

  public void set_end (int v) {
    range.end = v;
  }

  public void parse (String region) throws ParseException {
    // chr17:7571720-7590863
    boolean ok = false;
    String[] stuff = region.split(":");
    if (stuff.length == 2) {
      tname = stuff[0];
      String[] pos = stuff[1].split("-");
      if (pos.length == 2) {
	set_start(Integer.parseInt(pos[0]));
	set_end(Integer.parseInt(pos[1]));
	ok = true;
      }
    }

    if (!ok) throw new ParseException("parse error for " + region, 0);
  }

}