package org.stjude.compbio.rnapeg;

import java.io.*;

public abstract class EndianReader {

  protected RandomAccessFile raf;

  public EndianReader (RandomAccessFile raf) {
    this.raf = raf;
  }

  public abstract int readInt() throws IOException;

  public void seek (int i) throws IOException {
    raf.seek(i);
  }

  public long getFilePointer() throws IOException {
    return raf.getFilePointer();
  }

  public RandomAccessFile get_randomaccessfile() {
    return raf;
  }

  public void set_randomaccessfile(RandomAccessFile raf) {
    this.raf = raf;
  }
  

}
