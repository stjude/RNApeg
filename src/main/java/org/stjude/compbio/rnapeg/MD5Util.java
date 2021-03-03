package org.stjude.compbio.rnapeg;
import java.security.MessageDigest;

public class MD5Util {
  // maybe should just subclass MessageDigest?

  private MessageDigest md5;
  private byte[] digest;

  public MD5Util() throws java.security.NoSuchAlgorithmException {
    md5 = MessageDigest.getInstance("MD5");
    md5.reset();
    digest = null;
  }

  public void reset() {
    md5.reset();
  }

  public byte[] digest() {
    return digest = md5.digest();
  }

  public byte[] get_digest() {
    // like digest() but saves value
    if (digest == null) digest();
    return digest;
  }

  public void update (int v) {
    md5.update(string2ba(Integer.toString(v)));
  }

  public void update (byte[] v) {
    md5.update(v);
  }

  public static String digest2hex(byte[] digest) {
    int temp;
    StringBuilder hexb = new StringBuilder();
    for (int i = 0; i < digest.length; i++) {
      temp = 0xFF & digest[i];
      String s = Integer.toHexString(temp);
      if (temp <= 0x0F) s="0"+s;
      hexb.append(s);
    }
    return hexb.toString();
  }

  public String get_hex_digest() {
    if (digest == null) digest();
    return digest2hex(digest);
  }

  private byte[] string2ba (String ins) {
    char[] input = ins.toCharArray();
    byte[] output = new byte[input.length];
    for (int i=0; i < input.length; i++) {
      output[i] = (byte) input[i];
    }
    return output;
  }


  
}
