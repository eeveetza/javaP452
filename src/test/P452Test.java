package test;

import main.P452;
import main.P452DigitalMaps;

import org.junit.Before;
import org.junit.Test;
import org.junit.Assert;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;


public class P452Test {
    // the results are compared to the reference MATLAB implementation of Recommendation ITU-R P.452-18
    // the test is passed when the results for transmission loss are within the tolerance from the reference
    // for different antenna heights, time percentages, frequencies, distances, and clutter profiles
    //
    //     Rev   Date        Author                          Description
    //     -------------------------------------------------------------------------------
    //     v0    19MAY24     Ivica Stevanovic, OFCOM         Implementation for ITU-R P.452-18


    TestUtil util;

    @Before
    public void setup() {

        util = new TestUtil(1e-6);
    }

    @Test
    public void test() {

        P452DigitalMaps maps = new P452DigitalMaps();

        P452 calculator = new P452();

        int sizeY = 0;

        // path to the directory where profiles are located
        String directoryPath = "src/test/validation_examples/profiles/";

        // Using File class create an object for specific directory
        File directory = new File(directoryPath);

        // Using listFiles method we get all the files of a directory
        // return type of listFiles is array
        File[] files = directory.listFiles();

        // Get name of the all files present in that path
        if (files != null) {
            for (File file : files) {
                //System.out.println(file.getName());

                List<String> lines = new ArrayList<>();

                String file_rel = directoryPath + file.getName();

                // read all the lines from the profile
                try {

                    InputStream inputStream = new FileInputStream(file_rel);
                    InputStreamReader inputStreamReader = new InputStreamReader(inputStream);
                    BufferedReader br = new BufferedReader(inputStreamReader);
                    String line;
                    while (null != (line = br.readLine())) {

                        lines.add(line);

                    }

                    sizeY = lines.size();
                    inputStream.close();

                    double[] d = new double[sizeY - 1];
                    double[] h = new double[sizeY - 1];
                    double[] g = new double[sizeY - 1];
                    int[] zone = new int[sizeY - 1];

                    // extract the vectors for distance, terrain height, clutter height and zone
                    // noting that the first line is the table header

                    for (int i = 1; i < sizeY; i++) { /* DO */

                        String[] parts = lines.get(i).trim().split(",");
                        d[i - 1] = Double.parseDouble(parts[0]);
                        h[i - 1] = Double.parseDouble(parts[1]);
                        g[i - 1] = Double.parseDouble(parts[2]) + h[i - 1];
                        zone[i - 1] = Integer.parseInt(parts[4]);

                    }

                    // at this point, the profile is read
                    // read the reference results next

                    List<String> lines_r = new ArrayList<>();
                    String directoryPath_r = "src/test/validation_examples/results/";

                    String file_rel_r = directoryPath_r + file.getName().replaceFirst("profile", "result");
                    //System.out.println(file_rel_r);

                    InputStream inputStream_r = new FileInputStream(file_rel_r);
                    InputStreamReader inputStreamReader_r = new InputStreamReader(inputStream_r);
                    BufferedReader br_r = new BufferedReader(inputStreamReader_r);
                    String line_r;
                    while (null != (line_r = br_r.readLine())) {

                        lines_r.add(line_r);
                    }

                    sizeY = lines_r.size();
                    inputStream_r.close();

                    for (int i = 1; i < sizeY; i++) { /* DO */

                        String[] parts = lines_r.get(i).trim().split(",");

                        double f = Double.parseDouble(parts[1]);
                        double p = Double.parseDouble(parts[2]);

                        double htg = Double.parseDouble(parts[3]);
                        double hrg = Double.parseDouble(parts[4]);
                        double phit_e = Double.parseDouble(parts[5]);
                        double phit_n = Double.parseDouble(parts[6]);
                        double phir_e = Double.parseDouble(parts[7]);
                        double phir_n = Double.parseDouble(parts[8]);

                        double Gt = Double.parseDouble(parts[9]);
                        double Gr = Double.parseDouble(parts[10]);
                        double pol = Double.parseDouble(parts[11]);
                        double dct = Double.parseDouble(parts[12]);
                        double dcr = Double.parseDouble(parts[13]);


                        double press = Double.parseDouble(parts[14]);
                        double temp = Double.parseDouble(parts[15]);

                        double Lb_ref = Double.parseDouble(parts[37]);

                        double result = calculator.tl_p452(maps, f, p, d, h, g, zone, htg, hrg, phit_e, phit_n, phir_e, phir_n, Gt, Gr, pol, dct, dcr, press, temp, false);
                        util.assertDoubleEquals(Lb_ref, result);

                    }


                } catch (Exception ex) {

                    throw new IllegalArgumentException("Could not load the file: '" + file_rel + "'");
                }


            }


        } else {
            System.out.println("Did not find any files");
        }


    }

}