package org.stjude.compbio.rnapeg;

import java.sql.*;
import java.util.prefs.*;
import java.util.*;

public class JDBCQuery {
  private String DB_SERVER = "my_server";
  private String DB_DATABASE = "my_database";
  private String DB_USER = "my_user";
  private String DB_PASSWORD = "my_password";

  public static String UCSC_DB_SERVER = "genome-mysql.cse.ucsc.edu";
  public static String UCSC_DB_DATABASE = "hg18";
  public static String UCSC_DB_USERNAME = "genome";
  public static String UCSC_DB_PASSWORD = "";
  public static String UCSC_SNP_TABLE = "snp129";

  private ArrayList<HashMap<String,String>> results;
  private Connection connection = null;

  public static JDBCQuery get_ucsc_genome_client() {
    //
    // UCSC public server (http://genome.ucsc.edu/FAQ/FAQdownloads#download29)
    //
    JDBCQuery c = new JDBCQuery();
    c.set_db_server(UCSC_DB_SERVER);
    c.set_db_database(UCSC_DB_DATABASE);
    c.set_db_user(UCSC_DB_USERNAME);
    c.set_db_password(UCSC_DB_PASSWORD);
    return c;
  }

  public static JDBCQuery get_stjude_hg19() {
    //
    // St. Jude mirror of hg19
    //
    JDBCQuery c = new JDBCQuery();
    c.set_db_server("sjmemgb01");
    c.set_db_database("hg19");
    c.set_db_user("pallasro");
    c.set_db_password("pallasR0");
    return c;
  }

  public void set_db_server (String s) {
    DB_SERVER = s;
  }
  public void set_db_database (String s) {
    DB_DATABASE = s;
  }
  public void set_db_user (String s) {
    DB_USER = s;
  }
  public void set_db_password (String s) {
    DB_PASSWORD = s;
  }
  
  private Connection get_connection() throws SQLException {
    if (connection == null) {
      String cs = "jdbc:mysql://" + DB_SERVER + "/" + DB_DATABASE + "?user=" + DB_USER;
      if (DB_PASSWORD != null && DB_PASSWORD.length() > 0) {
	cs = cs.concat("&password=" + DB_PASSWORD);
      }
      //      System.err.println("cs="+cs);  // debug
      connection = DriverManager.getConnection(cs);
    }
    return connection;
  }

  public byte[] string2ba (String ins) {
    char[] input = ins.toCharArray();
    byte[] output = new byte[input.length];
    for (int i=0; i < input.length; i++) {
      output[i] = (byte) input[i];
    }
    return output;
  }

  public ArrayList<HashMap<String,String>> query_column (String table, String column, Collection<String> values, boolean is_numeric) throws Exception {
    // query for all rows where a column matches specified values
    ArrayList<String> mapped = new ArrayList<String>();
    if (is_numeric) {
      // verbatim
      for (String v : values) {
	mapped.add(v);
      }
    } else {
      // add quotes
      for (String v : values) {
	mapped.add('"' + v + '"');
      }
    }
    System.err.println("values size="+values.size());  // debug

    String query = "select * from " + table + " where " + column + " in (" +
      Str.join(",", mapped) + ")";
    System.err.println("query="+query);  // debug

    return query(query, false, 0);
  }

  public HashMap<String,String> query_single_row (String table, String column, String value) throws Exception {
    String query = "select * from " + table + " where " + column + "=\"" + value + "\"";
    return query_single_row(query);
  }


  public HashMap<String,String> query_single_row (String query) throws Exception {
    HashMap<String,String> result = null;
    ArrayList<HashMap<String,String>> results = query(query, true, 1);
    if (results != null && results.size() > 0) result = results.get(0);
    return result;
  }

  public ArrayList<HashMap<String,String>> query (String query) throws Exception {
    return query(query, false, 0);
  }

  public ArrayList<HashMap<String,String>> query (String query, boolean is_limited, int max_rows) throws Exception {
      results = new ArrayList<HashMap<String,String>>();
      ArrayList<String> columns = null;

      Connection c = get_connection();
      Statement st = c.createStatement();
      if (st.execute(query)) {
	ResultSet rs = st.getResultSet();

	//
	//  column setup:
	//
	ResultSetMetaData meta = rs.getMetaData();
	int cc = meta.getColumnCount();
	columns = new ArrayList<String>();
	for (int i=1; i <= cc; i++) {
	  columns.add(meta.getColumnName(i));
	}

	//
	//  fetch database results:
	//
	int row_count = 0;
	while (rs.next()) {
	  HashMap<String,String> row = new HashMap<String,String>();
	  String value;
	  ArrayList<String> raw = new ArrayList<String>();

	  if (true) {
	    for (int i=1; i <= cc; i++) {
	      value = new String(rs.getString(i));
	      raw.add(value);
	      row.put(columns.get(i - 1), value);
	    }
	  } else {
	    // old method: fetch data by column name.
	    // this catches fire, explodes, and falls over when running a 
	    // "show tables" command.
	    for (String col : columns) {
	      value = new String(rs.getString(col));
	      raw.add(value);
	      row.put(col, value);
	    }
	  }
	  results.add(row);

	  if (is_limited && ++row_count >= max_rows) break;

	}
      } else {
	System.err.println("no results!");  // debug
      }
      return results;
  }

  public int get_result_count() {
    return results.size();
  }
  
  public ArrayList<HashMap<String,String>> get_results() {
    return results;
  }

  public static void main (String[] argv) {
    //   JDBCQuery ucsc = JDBCQuery.get_ucsc_genome_client();
    System.err.println("time: " + System.currentTimeMillis());  // debug
    JDBCQuery ucsc = JDBCQuery.get_stjude_hg19();
    System.err.println("time: " + System.currentTimeMillis());  // debug
    try {
      ArrayList<HashMap<String,String>> results = ucsc.query("select * from refGene where name2=\"tp53\"");
      System.err.println("time: " + System.currentTimeMillis());  // debug
      for (HashMap<String,String> row : results) {
	System.err.println(row.get("name2"));  // debug
      }
      
      System.err.println("count=" + ucsc.get_result_count());  // debug
    } catch (Exception e) {
      System.err.println("ERROR: "+ e);  // debug
      e.printStackTrace();
    }
  }


}